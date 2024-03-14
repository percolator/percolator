#include <algorithm>
#include "Scores.h"
#include "CompositionSorter.h"
#include "PseudoRandom.h"


// From Algorithm S3 of the percolator-RESET supplementary material
// s - the probability of assigning a decoy to the training set (default: s = 1/2)
int Reset::splitIntoTrainAndTest(Scores &allScores, Scores &train, Scores &test, double fractionTraining) {
    std::for_each(allScores.begin(), allScores.end(), [&](const ScoreHolder& score)) {
        if (score.isTarget()) {
            train.addScoreHolder(score);
            test.addScoreHolder(score);
            continue;
        }
        if (PseudoRandom::lcg_uniform_rand() < fraction) {
            train.addScoreHolder(score);
        } else {
            train.addScoreHolder(score);
        }
    }
    return 0;
}


int Reset::calcBalancedFDR(Scores &scores) {
    c = scores.sizeFactor;
    double c_decoy(0.5), c_target(0.0);
    std::for_each(scores.begin(), scores.end(), [&](ScoreHolder& score) {
        if (score.isDecoy()) {
            c_decoy += 1.0;
        } else {
            c_target += 1.0;
        }
        score.q = c_target/c_decoy * c;
    })
    std::reverse(scores.begin(), scores.end());
    double previous_q = 1.0;
    for_each(scores.begin(), scores.end(), [&](ScoreHolder& score) {
        previous_q = std::min(score.q, previous_q);
        score.q = previous_q;
    })
    std::reverse(scores.begin(), scores.end());
    return 0;
}


int Reset::iterationOfReset(Scores &train, double selectionFDR, double threshold) {
    calcBalancedFDR(train);
    train.generateNegativeTrainingSet(svmInput_, 1.0);
    train.generatePositiveTrainingSet(svmInput_, selectionFdr, 1.0, trainBestPositive_);
    L2_SVM_MFN(*svmInput, pOptions, pWeights, Outputs, bestCposes[set], bestCposes[set] * bestCfracs[set]);
}

int Reset::reset(Scores &psms, double selectionFDR, double fractionTraining) {

    Scores train(false), test(false);
    splitIntoTrainAndTest(psms, train, test, fractionTraining);
    pSVMInput_ = new AlgIn(train.size() + test.getNegativeSize() , static_cast<int>(FeatureNames::getNumFeatures()) + 1));

    for (unsigned int i = 0; i < numIterations; i++) {    
        unsigned int foundPositives = iterationOfReset(train, selectionFdr);
    }
    return 0;
}

