#include <algorithm>
#include "Scores.h"
#include "SanityCheck.h"
#include "CompositionSorter.h"
#include "PseudoRandom.h"
#include "Normalizer.h"
#include "ssl.h"
#include "Reset.h"


// From Algorithm S3 of the percolator-RESET supplementary material
// s - the probability of assigning a decoy to the training set (default: s = 1/2)
int Reset::splitIntoTrainAndTest(Scores &allScores, Scores &train, Scores &test, double fractionTraining) {

    cerr << "Inside split" << endl;

    std::for_each(allScores.begin(), allScores.end(), [&](const ScoreHolder& score) {
        cerr << score.score << " " << score.isTarget() << endl;
        if (score.isTarget()) {
            train.addScoreHolder(score);
            test.addScoreHolder(score);
        } else if (PseudoRandom::lcg_uniform_rand() < fractionTraining) {
            train.addScoreHolder(score);
        } else {
            test.addScoreHolder(score);
        }
    });
    return 0;
}

// Same as above but testing with pointers to SH instead of the direct objects 
// From Algorithm S3 of the percolator-RESET supplementary material
// s - the probability of assigning a decoy to the training set (default: s = 1/2)
int Reset::splitIntoTrainAndTest(Scores &allScores, vector<ScoreHolder*> &train, vector<ScoreHolder*> &test, double fractionTraining) {

    cerr << "Inside split" << endl;

    std::for_each(allScores.begin(), allScores.end(), [&](ScoreHolder& score) {
        cerr << score.score << " " << score.isTarget() << endl;
        if (score.isTarget()) {
            train.push_back(&score);
            test.push_back(&score);
        } else if (PseudoRandom::lcg_uniform_rand() < fractionTraining) {
            train.push_back(&score);
        } else {
            test.push_back(&score);
        }
    });
    return 0;
}


int Reset::iterationOfReset(Scores &train, double selectionFDR) {
    train.calcBalancedFDR(selectionFDR);
    train.generateNegativeTrainingSet(*pSVMInput_, 1.0);
    train.generatePositiveTrainingSet(*pSVMInput_, selectionFDR, 1.0, false);

    vector_double pWeights, Outputs;
    pWeights.d = static_cast<int>(FeatureNames::getNumFeatures()) + 1;
    pWeights.vec = new double[pWeights.d];
    
    size_t numInputs = static_cast<std::size_t>(pSVMInput_->positives + pSVMInput_->negatives);
    Outputs.vec = new double[numInputs];
    Outputs.d = static_cast<int>(numInputs);
    
    for (int ix = 0; ix < pWeights.d; ix++) {
        pWeights.vec[ix] = 0;
    }
    for (int ix = 0; ix < Outputs.d; ix++) {
        Outputs.vec[ix] = 0;
    }

    // TODO: Set CPos,CNeg by X-val

    L2_SVM_MFN(*pSVMInput_, 
        options_, 
        pWeights, 
        Outputs, 
        1.0, // Cpos
        3.0 * train.getTargetDecoySizeRatio()); // CNeg
    

    for (std::size_t i = FeatureNames::getNumFeatures() + 1; i--;) {
        w_[i] = static_cast<double>(pWeights.vec[i]);
    }

    train.onlyCalcScores(w_);
    return train.calcBalancedFDR(selectionFDR);
}

int Reset::reset(Scores &psms, double selectionFDR, SanityCheck* pCheck, double fractionTraining, unsigned int decoysPerTarget) {

    CompositionSorter sorter;

    // select the representative peptides to train on
    Scores winnerPeptides(false);

    cerr << "Starting reset: psmAndPeptide" << endl;
    
    sorter.psmAndPeptide(psms, winnerPeptides, decoysPerTarget);

    cerr << "Splitting into train/test" << endl;

    // Setting up the training and test sets
    double factor = 1.0*(decoysPerTarget + 1)  - fractionTraining*decoysPerTarget;
    unsigned int numIterations(5);
    Scores train(false), test(false);
    splitIntoTrainAndTest(winnerPeptides, train, test, fractionTraining);
    train.setNullTargetWinProb(factor/(1.0 + decoysPerTarget));
    test.setNullTargetWinProb(1.0/factor);
    cerr << "Setting sizes" << endl;

    train.recalculateSizes(); test.recalculateSizes();

    pCheck->getInitDirection(train, Normalizer::getNormalizer(), w_, selectionFDR, selectionFDR);
    train.getInitDirection(selectionFDR, w_);
    train.onlyCalcScores(w_);

    // Initialize the input for the SVM

    cerr << "Setting up SVM training" << endl;
    pSVMInput_ = new AlgIn(train.size() + test.negSize() , static_cast<int>(FeatureNames::getNumFeatures()) + 1);

    cerr << "Training" << endl;
    for (unsigned int i = 0; i < numIterations; i++) {    
        unsigned int foundPositives = iterationOfReset(train, selectionFDR);
    }
    return 0;
}

