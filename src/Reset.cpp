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
        // cerr << score.score << " " << score.isTarget() << endl;
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
int Reset::splitIntoTrainAndTest(std::vector<ScoreHolder*> &allScores, vector<ScoreHolder*> &train, vector<ScoreHolder*> &test, double fractionTraining) {

    cerr << "Inside split" << endl;

    for(auto pScore : allScores) {
        // cerr << pScore->score << " " << pScore->isTarget() << endl;
        if (pScore->isTarget()) {
            train.push_back(pScore);
            test.push_back(pScore);
        } else if (PseudoRandom::lcg_uniform_rand() < fractionTraining) {
            train.push_back(pScore);
        } else {
            test.push_back(pScore);
        }
    }
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

int calcBalancedFDR(vector<ScoreHolder*> &scores, double nullTargetWinProb, double treshold) {
    double c_decoy(0.5), c_target(0.0), factor( nullTargetWinProb / ( 1.0 - nullTargetWinProb ) );
    for_each(scores.begin(), scores.end(), [&](ScoreHolder* pScore) {
        if (pScore->isDecoy()) {
            c_decoy += 1.0;
        } else {
            c_target += 1.0;
        }
        pScore->q = ( c_target / c_decoy ) * factor;
    });
    std::reverse(scores.begin(), scores.end());
    double previous_q = 1.0;
    for_each(scores.begin(), scores.end(), [&](ScoreHolder* pScore) {
        previous_q = std::min(pScore->q, previous_q);
        pScore->q = previous_q;
    });
    std::reverse(scores.begin(), scores.end());

    // Calcultaing the number of target PSMs with q-value less than treshold
    ScoreHolder limit(treshold, 1);
    auto pointerComp = [](const ScoreHolder* a, const ScoreHolder* b) { return *a < *b; };
    auto upper = std::lower_bound(scores.begin(), scores.end(), &limit, pointerComp);
    return upper - scores.begin();
}

void generateNegativeTrainingSet(AlgIn& data, std::vector<ScoreHolder*> scores, const double cneg) {
    std::size_t ix2 = 0;
    // Using range-based for loop with auto
    for (const auto scorePtr : scores) {  // scorePtr is automatically a ScoreHolder* from scores
        cerr << scorePtr->label << " " << scorePtr->isTarget() << " " << scorePtr->score << endl;
        if (!scorePtr->isTarget()) {  // Directly use scorePtr, which is a ScoreHolder*
            data.vals[ix2] = scorePtr->pPSM->features;  // Access members directly through the pointer
            data.Y[ix2] = -1;
            data.C[ix2++] = cneg;
        }
    }
    data.negatives = static_cast<int>(ix2);
    cerr << "Generated negative training set of size " << ix2 << endl;
}

void generatePositiveTrainingSet(AlgIn& data, const std::vector<ScoreHolder*>& scores,
                                         const double fdr, const double cpos,
                                         const bool trainBestPositive) {
    std::size_t ix2 = static_cast<std::size_t>(data.negatives);
    int p = 0;

    auto last = scores.end();
    // Range-based loop over the sorted and uniqued scores
    for (auto it = scores.begin(); it != last; ++it) {
        const auto scorePtr = *it;  // Dereferencing iterator to get the pointer to ScoreHolder
        if (scorePtr->isTarget()) {
            if (scorePtr->q <= fdr) {
                data.vals[ix2] = scorePtr->pPSM->features;
                data.Y[ix2] = 1;
                data.C[ix2++] = cpos;
                ++p;
            }
        }
    }

    data.positives = p;
    data.m = static_cast<int>(ix2);
    cerr << "Generated positive training set of size " << p << endl;
}

double calcScore(const double* feat, const std::vector<double>& w) {
    std::size_t ix = FeatureNames::getNumFeatures();
    double score = w[ix];
    for (; ix--;) {
        score += feat[ix] * w[ix];
    }
    return score;
}

int onlyCalcScores(std::vector<ScoreHolder*> scores, std::vector<double>& w) {
    std::size_t ix;
    for (auto scorePtr : scores) {  // scorePtr is automatically a ScoreHolder* from scores
      scorePtr->score = calcScore(scorePtr->pPSM->features, w);
    }
    auto reversePointerComp = [](const ScoreHolder* a, const ScoreHolder* b) { return *a > *b; }; // Pointer GreaterThen
    sort(scores.begin(), scores.end(),reversePointerComp);
    return 0;
}

int calcScores(std::vector<ScoreHolder*> scores, std::vector<double>& w, double fdr, bool skipDecoysPlusOne) {
    onlyCalcScores(scores, w);
    return calcBalancedFDR(scores, fdr, skipDecoysPlusOne);
}

// Sets init direction to the feature giving best discrimination
int setInitDirection(vector<double>& w, vector<ScoreHolder*>& scores, double nullTargetWinProb, double selectionFDR) {
    int maxIds(-1), maxDir(0), maxSign(0);

    for (int ix = 0; ix < w.size(); ix++) {
        for(auto scorePtr : scores) {
            scorePtr->score = scorePtr->pPSM->features[ix];
        }
        std::sort(scores.begin(), scores.end(), [](const ScoreHolder* a, const ScoreHolder* b) { return *a > *b; });
        for (int dir = 1; dir > -2; dir -= 2) {
            int numId = calcBalancedFDR(scores, nullTargetWinProb, selectionFDR);
            if (numId > maxIds) {
                maxIds = numId;
                maxDir = ix;
                maxSign = dir;
            }
            if(dir == 1) {
                std::reverse(scores.begin(), scores.end());
            }
        }
    }
    for (int ix = 0; ix < w.size(); ix++) {
        w[ix] = 0;
    }
    w[maxDir] = static_cast<double>(maxSign);
    return calcScores(scores, w, selectionFDR, false);
}


int Reset::iterationOfReset(vector<ScoreHolder*> &train, double nullTargetWinProb, double selectionFDR) {
    cerr << "Inside iterationOfReset" << endl;
    calcBalancedFDR(train, nullTargetWinProb, selectionFDR);
    cerr << "Generating training sets" << endl;
    generateNegativeTrainingSet(*pSVMInput_, train, 1.0);
    generatePositiveTrainingSet(*pSVMInput_, train, selectionFDR, 1.0, false);


    vector_double pWeights, Outputs;
    pWeights.d = static_cast<int>(FeatureNames::getNumFeatures()) + 1;
    pWeights.vec = new double[pWeights.d];
    cerr << "Made weight vector of length " << pWeights.d << endl;
    
    size_t numInputs = static_cast<std::size_t>(pSVMInput_->positives + pSVMInput_->negatives);
    Outputs.vec = new double[numInputs];
    Outputs.d = static_cast<int>(numInputs);
    cerr << "Made output vector of length " << Outputs.d << endl;
    
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
        3.0 * pSVMInput_->positives / max(1.0, (double) pSVMInput_->negatives)); // CNeg
    

    for (std::size_t i = FeatureNames::getNumFeatures() + 1; i--;) {
        w_[i] = static_cast<double>(pWeights.vec[i]);
    }

    onlyCalcScores(train, w_);
    return calcBalancedFDR(train, nullTargetWinProb, selectionFDR);
}

int Reset::reset(Scores &psms, double selectionFDR, SanityCheck* pCheck, double fractionTraining, unsigned int decoysPerTarget) {

    CompositionSorter sorter;

    // select the representative peptides to train on
    std::vector<ScoreHolder*> winnerPeptides;

    cerr << "Starting reset: psmAndPeptide" << endl;
    
    sorter.psmAndPeptide(psms, winnerPeptides, decoysPerTarget);

    cerr << "Splitting into train/test" << endl;

    // Setting up the training and test sets
    double factor = 1.0*(decoysPerTarget + 1)  - fractionTraining*decoysPerTarget;
    unsigned int numIterations(5);
    std::vector<ScoreHolder*> train(false), test(false);
    splitIntoTrainAndTest(winnerPeptides, train, test, fractionTraining);
    double trainNullTargetWinProb(factor/(1.0 + decoysPerTarget)), testNullTargetWinProb(1/factor);

    // setInitDirection(w_, train, trainNullTargetWinProb, selectionFDR);

    // Initialize the input for the SVM

    cerr << "Setting up SVM training for a size of " << winnerPeptides.size() << " peptides." << endl;
    // pSVMInput_ = new AlgIn(train.size() + test.negSize() , static_cast<int>(FeatureNames::getNumFeatures()) + 1);
    assert(pSVMInput_ == nullptr);
    pSVMInput_ = new AlgIn(winnerPeptides.size(), static_cast<int>(FeatureNames::getNumFeatures()) + 1);

    cerr << "Training" << endl;
    for (unsigned int i = 0; i < numIterations; i++) {    
        unsigned int foundPositives = iterationOfReset(train, trainNullTargetWinProb, selectionFDR);
    }
    return 0;
}

