#include <algorithm>
#include "Scores.h"
#include "SanityCheck.h"
#include "CompositionSorter.h"
#include "PosteriorEstimator.h"
#include "PseudoRandom.h"
#include "Normalizer.h"
#include "ssl.h"
#include "DataSet.h"
#include "Reset.h"


// From Algorithm S3 of the percolator-RESET supplementary material
// s - the probability of assigning a decoy to the training set (default: s = 1/2)
int Reset::splitIntoTrainAndTest(Scores &allScores, Scores &train, Scores &test, double fractionTraining) {

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

    for(auto& pScore : allScores) {
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

int calcBalancedFDR(vector<ScoreHolder*> &scores, double nullTargetWinProb, double treshold) {
    // This function requies the scoreholders to be sorted in score decending order
    // prior to function call. 
    double c_decoy(0.5), c_target(0.0), factor( nullTargetWinProb / ( 1.0 - nullTargetWinProb ) );
    for(auto& pScore : scores) {
        if (pScore->isDecoy()) {
            c_decoy += 1.0;
        } else {
            c_target += 1.0;
        }
        pScore->q = ( c_decoy / c_target ) * factor;
    }
    std::reverse(scores.begin(), scores.end());
    double previous_q = 1.0;
    for(auto& pScore : scores) {
        previous_q = std::min(pScore->q, previous_q);
        pScore->q = previous_q;
    }
    std::reverse(scores.begin(), scores.end());

    // Calcultaing the number of target PSMs with q-value less than treshold
    ScoreHolder limit; limit.q = treshold;
    auto pointerComp = [](const ScoreHolder* a, const ScoreHolder* b) { return a->q < b->q; };
    auto upper = std::lower_bound(scores.begin(), scores.end(), &limit, pointerComp);
    cerr << "Balanced FDR reports " << upper - scores.begin() << " of " << scores.end() - scores.begin() << " scoreholders over FDR treshold. ";
    cerr << (*upper)->q << ' ' << (*upper)->score << " ";
    cerr << (*upper-1)->q << endl;
    return upper - scores.begin();
}

void generateTrainingSet(AlgIn& data, std::vector<ScoreHolder*>& scores, const double cneg, const double cpos, const double fdr) {
    std::size_t ix2 = 0;
    // Using range-based for loop with auto
    for (const auto scorePtr : scores) {  // scorePtr is automatically a ScoreHolder* from scores
        if (!scorePtr->isTarget()) {  // Directly use scorePtr, which is a ScoreHolder*
            data.vals[ix2] = scorePtr->pPSM->features;  // Access members directly through the pointer
            data.Y[ix2] = -1;
            // data.C[ix2] = cneg;
            ix2++;
        }
    }
    data.negatives = static_cast<int>(ix2);
    int p = 0;
   // Range-based loop over the sorted and uniqued scores
    for (const auto scorePtr : scores) {  // scorePtr is automatically a ScoreHolder* from scores
        // cerr << scorePtr->label << " " << scorePtr->isTarget() << " " << scorePtr->score << " " << scorePtr->q << endl;
        if (scorePtr->isTarget()) {
            if (scorePtr->q <= fdr) {
                data.vals[ix2] = scorePtr->pPSM->features;
                data.Y[ix2] = 1;
                // data.C[ix2] = cpos;
                ix2++;
                ++p;
            }
        }
    }
    data.positives = p;
    data.m = static_cast<int>(ix2);
}


double calcScore(const double* feat, const std::vector<double>& w) {
    std::size_t ix = FeatureNames::getNumFeatures();
    double score = w[ix];
    for (; ix--;) {
        score += feat[ix] * w[ix];
    }
    return score;
}

// This function assign score to all score holders based on the weight vector w
// and sorts the score holders in decending order
int onlyCalcScores(std::vector<ScoreHolder*> scores, std::vector<double>& w) {
    std::size_t ix;
    for (auto& scorePtr : scores) {  // scorePtr is automatically a ScoreHolder* from scores
      scorePtr->score = calcScore(scorePtr->pPSM->features, w);
    }
    auto reversePointerComp = [](const ScoreHolder* a, const ScoreHolder* b) { return *a > *b; }; // Pointer GreaterThen
    sort(scores.begin(), scores.end(),reversePointerComp);
    return 0;
}

// This function assign score to all score holders based on the weight vector w
// and sorts the score holders in decending order and reports the number of scor holders over a FDR treshold.
int calcScores(std::vector<ScoreHolder*> scores, std::vector<double>& w, double fdr, bool skipDecoysPlusOne) {
    onlyCalcScores(scores, w);
    return calcBalancedFDR(scores, fdr, skipDecoysPlusOne);
}

// Sets init direction to the feature giving best discrimination
int setInitDirection(vector<double>& w, vector<ScoreHolder*>& scores, double nullTargetWinProb, double selectionFDR) {
    int maxIds(-1), maxDir(0), maxSign(0);

    for (int ix = 0; ix < w.size(); ix++) {
        for(auto& scorePtr : scores) {
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

void getScoreLabelPairs(std::vector<ScoreHolder*> scores, std::vector<pair<double, bool> >& combined) {
    combined.clear(); // Clear the combined vector first
    // Use a lambda function to transform each ScoreHolder pointer to a pair<double, bool>
    std::transform(scores.begin(), scores.end(), std::back_inserter(combined),
        [](ScoreHolder* sh) { return sh->toPair(); }); // Assuming toPair() is a member function of ScoreHolder
}


int Reset::gridSearchC(vector<ScoreHolder*> &train, const double nullTargetWinProb, const double selectionFDR) {
    std::sort(train.begin(),train.end(), [](const ScoreHolder* a, const ScoreHolder* b) { return *a > *b; });
    calcBalancedFDR(train, nullTargetWinProb, selectionFDR);
    generateTrainingSet(*pSVMInput_, train, 1.0, 1.0, selectionFDR);

    int bestResult(0);
    double bestCpos(0.), bestCfrac(0.);

    std::vector<double> cPosCandidates = {10000.,10.,1.0,0.1,0.0001};
    std::vector<double> cFracCandidates = {0.001, 0.3, 1.0, 3.0, 1000.0};
    for (auto cPos : cPosCandidates) {
        for (auto cFrac : cFracCandidates) {
            svmTrain(cPos, cFrac);
            onlyCalcScores(train, w_);
            int result = calcBalancedFDR(train, nullTargetWinProb, selectionFDR);
            if (VERB>2) 
                cerr << "GridSearch step with C+=" << cPos << ", and C-/C+=" << cFrac << " resulted in " << result << " peptides." << endl;  
            if (result>=bestResult) {
                bestResult = result;
                bestCpos = cPos;
                bestCfrac = cFrac;
            }
        }
    }
    if (VERB>1)
        cerr << "GridSearch found the optimal hyperParameters for SVM training, C+=" << bestCpos << ", and C-/C+=" << bestCfrac << "." << endl;  
    cPos_ = bestCpos;
    cFrac_ = bestCfrac;
}

void Reset::svmTrain(const double cPos, const double cFrac) {

    vector_double pWeights, Outputs;
    pWeights.d = static_cast<int>(FeatureNames::getNumFeatures()) + 1;
    pWeights.vec = new double[pWeights.d];
    
    size_t numInputs = static_cast<std::size_t>(pSVMInput_->positives + pSVMInput_->negatives);
    Outputs.vec = new double[numInputs];
    Outputs.d = static_cast<int>(numInputs);
    
    for (int ix = 0; ix < pWeights.d; ix++) {
        pWeights.vec[ix] = 0.;
    }
    for (int ix = 0; ix < Outputs.d; ix++) {
        Outputs.vec[ix] = 0.;
    }
    double cNeg = cFrac * cPos;

    L2_SVM_MFN(*pSVMInput_, 
        options_, 
        pWeights, 
        Outputs, 
        cPos, // CPos
        cNeg); // CNeg
    

    for (std::size_t i = FeatureNames::getNumFeatures() + 1; i--;) {
        w_[i] = static_cast<double>(pWeights.vec[i]);
        cerr << pWeights.vec[i] << '\t';
    }
    cerr << endl;
}

int Reset::iterationOfReset(vector<ScoreHolder*> &train, double nullTargetWinProb, double selectionFDR) {
    std::sort(train.begin(),train.end(), [](const ScoreHolder* a, const ScoreHolder* b) { return *a > *b; });
    calcBalancedFDR(train, nullTargetWinProb, selectionFDR);
    generateTrainingSet(*pSVMInput_, train, cPos_, cFrac_, selectionFDR);
    if (VERB > 1)
        cerr << "Training with a set of size " << pSVMInput_->m << " preptides, where of "  << pSVMInput_->positives << " are target, and " << pSVMInput_->negatives << " are decoys." << endl;

    svmTrain(cPos_,cFrac_);
    onlyCalcScores(train, w_);
    return calcBalancedFDR(train, nullTargetWinProb, selectionFDR);
}

int Reset::evaluateTestSet(Scores &psms, vector<ScoreHolder*> &test, double testNullTargetWinProb, double selectionFDR) {

    onlyCalcScores(test, w_); // Sorts the scores in descending order
    std::sort(test.begin(),test.end(), [](const ScoreHolder* a, const ScoreHolder* b) { return a->score > b->score; }); // Do we realy need this?
    int peptidesUnderFDR = calcBalancedFDR(test, testNullTargetWinProb, selectionFDR);
    // psms.onlyCalcScores(w_);

    std::vector<pair<double, bool> > combined;
    getScoreLabelPairs(test,combined);

    std::vector<double> peps;
    // Logistic regression on the data
    double factor( testNullTargetWinProb / ( 1.0 - testNullTargetWinProb ));
    if (VERB>2) cerr << "EstimatePEP, using factor=" << factor << endl;
    PosteriorEstimator::estimateTradPEP(combined, factor, peps,  true);
    return peptidesUnderFDR;
}

int Reset::reset(Scores &psms, Scores &outS, double selectionFDR, SanityCheck* pCheck, double fractionTraining, unsigned int decoysPerTarget, std::vector<double>& w, bool use_composition_match) {

    w_ = w;
    CompositionSorter sorter;

    // select the representative peptides to train on
    std::vector<ScoreHolder*> winnerPeptides;

    if (use_composition_match) {
        cerr << "Starting reset: psmAndPeptide" << endl;   
        sorter.psmAndPeptide(psms, winnerPeptides, decoysPerTarget);
    } else {
        cerr << "Starting reset: psmsOnly" << endl;   
        sorter.psmsOnly(psms, winnerPeptides);
    }
    std::cerr << "Splitting into train/test" << endl;

    // Setting up the training and test sets
    double factor = 1.0*(decoysPerTarget + 1)  - fractionTraining*decoysPerTarget;
    unsigned int numIterations(5);
    std::vector<ScoreHolder*> train(false), test(false);
    splitIntoTrainAndTest(winnerPeptides, train, test, fractionTraining);
    double trainNullTargetWinProb(factor/(1.0 + decoysPerTarget)), testNullTargetWinProb(1/factor);

    // setInitDirection(w_, train, trainNullTargetWinProb, selectionFDR); done before this?

    // Initialize the input for the SVM

    if (VERB>1) std::cerr << "Setting up SVM training for a size of " << winnerPeptides.size() << " peptides." << endl;
    // pSVMInput_ = new AlgIn(train.size() + test.negSize() , static_cast<int>(FeatureNames::getNumFeatures()) + 1);
    assert(pSVMInput_ == nullptr);
    pSVMInput_ = new AlgIn(winnerPeptides.size(), static_cast<int>(FeatureNames::getNumFeatures()) + 1);

    if (VERB>1) std::cerr << "Training set prepared. Starting SVM Training." << endl;

    gridSearchC(train, trainNullTargetWinProb, selectionFDR);

    for (unsigned int i = 0; i < numIterations; i++) {    
        unsigned int foundPositives = iterationOfReset(train, trainNullTargetWinProb, selectionFDR);
    }
    if (VERB>1) std::cerr << "Training Done!" << endl;

    if (VERB>1) {
        std::cerr << "Final feature weights" << endl;
        for (size_t ix=0; ix < w_.size(); ix++) {
            if (ix<FeatureNames::getNumFeatures()) {
                std::cerr << DataSet::getFeatureNames().getFeatureName(ix) << "\t" << w_[ix] << endl;
            } else {
                std::cerr << "Intercept" << "\t" << w_[ix] << endl;
            }
        }
    }
    
    int underFDR = evaluateTestSet(psms, test, testNullTargetWinProb, selectionFDR);
    cerr << "Test set evaluation finds " << underFDR << " peptides under a FDR of " << selectionFDR << endl;

    for (ScoreHolder* ptr : test) {
        if (ptr != nullptr) { // Check to ensure the pointer is not null
            outS.addScoreHolder(*ptr); // Dereference the pointer and add to scores
        }
    }
    return 0;
}

int Reset::rereset(Scores &psms, Scores &outS, double selectionFDR, SanityCheck* pCheck, double fractionTraining, unsigned int decoysPerTarget, std::vector<double>& w) {

    w_ = w;

    double factor = 1.0*(decoysPerTarget + 1)  - fractionTraining*decoysPerTarget;
    double trainNullTargetWinProb(factor/(1.0 + decoysPerTarget)), testNullTargetWinProb(1/factor);
    unsigned int numIterations(5);
    std::vector<ScoreHolder*> train(false), test(false);

    for (unsigned int i = 0; i < numIterations; i++) {    
        CompositionSorter sorter;

        // select the representative peptides to train on
        std::vector<ScoreHolder*> winnerPeptides;

        cerr << "rereset: psmAndPeptide" << endl;
        psms.onlyCalcScores(w_);

        sorter.psmAndPeptide(psms, winnerPeptides, decoysPerTarget);

        cerr << "Splitting into train/test" << endl;

        // Setting up the training and test sets
        train.clear(); test.clear();
        splitIntoTrainAndTest(winnerPeptides, train, test, fractionTraining);
    
    
        // Initialize the input for the SVM

        if(!pSVMInput_) {
            cerr << "Setting up SVM training for a size of " << winnerPeptides.size() << " peptides." << endl;
            pSVMInput_ = new AlgIn(winnerPeptides.size(), static_cast<int>(FeatureNames::getNumFeatures()) + 1);
        }

        cerr << "Training" << endl;
        unsigned int foundPositives = iterationOfReset(train, trainNullTargetWinProb, selectionFDR);
    }
    cerr << "Training Done!" << endl;
    
    int underFDR = evaluateTestSet(psms, test, testNullTargetWinProb, selectionFDR);
    cerr << "Test set evaluation finds " << underFDR << " number of peptides under FDR of " << selectionFDR << endl;

    for (ScoreHolder* ptr : test) {
        if (ptr != nullptr) { // Check to ensure the pointer is not null
            outS.addScoreHolder(*ptr); // Dereference the pointer and add to scores
        }
    }
    return 0;
}


