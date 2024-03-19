#include <algorithm>
#include "Scores.h"
#include "CompositionSorter.h"
#include "PseudoRandom.h"
#include "Reset.h"
#include "ssl.h"


// From Algorithm S3 of the percolator-RESET supplementary material
// s - the probability of assigning a decoy to the training set (default: s = 1/2)
int Reset::splitIntoTrainAndTest(Scores &allScores, Scores &train, Scores &test, double fractionTraining) {
    std::for_each(allScores.begin(), allScores.end(), [&](const ScoreHolder& score) {
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


int Reset::iterationOfReset(Scores &train, double selectionFDR, double threshold) {
    train.calcBalancedFDR(selectionFDR);
    train.generateNegativeTrainingSet(*pSVMInput_, 1.0);
    train.generatePositiveTrainingSet(*pSVMInput_, selectionFdr, 1.0, trainBestPositive_);

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

int Reset::reset(Scores &psms, double selectionFDR, double fractionTraining, int decoysPerTarget) {

    double factor = 1.0*(decoysPerTarget + 1)  - fractionTraining*decoysPerTarget;
    unsigned int numIterations(5);
    Scores train(false), test(false);
    splitIntoTrainAndTest(psms, train, test, fractionTraining);
    train.setNullTargetWinProb(factor/(1.0 + decoysPerTarget));
    test.setNullTargetWinProb(1.0/factor);
    pSVMInput_ = new AlgIn(train.size() + test.getNegativeSize() , static_cast<int>(FeatureNames::getNumFeatures()) + 1));

    for (unsigned int i = 0; i < numIterations; i++) {    
        unsigned int foundPositives = iterationOfReset(train, selectionFdr);
    }
    return 0;
}

