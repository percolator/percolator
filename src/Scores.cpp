/*******************************************************************************
 Copyright 2006-2012 Lukas Käll <lukas.kall@scilifelab.se>

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

 *******************************************************************************/

#include <cassert>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <utility>
#include <algorithm>
#include <set>
#include <map>
#include <vector>
#include <string>
#include <cmath>
#include <memory>
using namespace std;
#include "DataSet.h"
#include "Normalizer.h"
#include "SetHandler.h"
#include "Scores.h"
#include "Globals.h"
#include "PosteriorEstimator.h"
#include "ssl.h"
#include "MassHandler.h"


inline bool operator>(const ScoreHolder& one, const ScoreHolder& other) {
  return (one.score > other.score);
}

inline bool operator<(const ScoreHolder& one, const ScoreHolder& other) {
  return (one.score < other.score);
}

inline double truncateTo(double truncateMe, const char* length) {
  char truncated[64];
  char format[64];
  strcpy(format,"%.");
  strcat(format,length);
  strcat(format,"lf\n");
  sprintf(truncated, format, truncateMe);
  return atof(truncated);
}

ostream& operator<<(ostream& os, const ScoreHolder& sh) {
  if (sh.isDecoy() && !Scores::getPrintDecoysInXml()) {
    return os;
  }
  
  os << "    <psm p:psm_id=\"" << sh.pPSM->id << "\"";
  if (Scores::getPrintDecoysInXml()) {
    if (sh.isDecoy())
      os << " p:decoy=\"true\"";
    else 
      os << " p:decoy=\"false\"";
  }
  os << ">" << endl;
  
  os << "      <svm_score>"   << fixed   << sh.score   << "</svm_score>" << endl;
  os << "      <q_value>"   << scientific << sh.q   << "</q_value>" << endl;
  os << "      <pep>"          << scientific << sh.pep << "</pep>" << endl;
  
  if(Scores::getShowExpMass()) 
  {
    os << "      <exp_mass>" << fixed << setprecision (4) << sh.pPSM->expMass << "</exp_mass>" << endl;
  }   
  
  os << "      <calc_mass>" << fixed << setprecision (3) << sh.pPSM->calcMass << "</calc_mass>" << endl;
  
  if (DataSet::getCalcDoc()) os << "      <retentionTime observed=\"" 
          << PSMDescription::unnormalize(sh.pPSM->retentionTime)
          << "\" predicted=\""
          << PSMDescription::unnormalize(sh.pPSM->predictedTime) << "\"/>"
          << endl;

  if (sh.pPSM->getPeptideSequence().size() > 0) {
    string n = sh.pPSM->getFlankN();
    string c = sh.pPSM->getFlankC();
    string centpep = sh.pPSM->getPeptideSequence();
    os << "      <peptide_seq n=\"" << n << "\" c=\"" << c << "\" seq=\"" << centpep << "\"/>" << endl;
  }
  
  std::set<std::string>::const_iterator pidIt = sh.pPSM->proteinIds.begin();
  for ( ; pidIt != sh.pPSM->proteinIds.end() ; ++pidIt) {
    os << "      <protein_id>" << getRidOfUnprintablesAndUnicode(*pidIt) << "</protein_id>" << endl;
  }
  
  os << "      <p_value>" << scientific << sh.p << "</p_value>" <<endl;
  os << "    </psm>" << endl;
  return os;
}

ostream& operator<<(ostream& os, const ScoreHolderPeptide& sh) {
  if (sh.isDecoy() && !Scores::getPrintDecoysInXml()) {
    return os;
  }
  
  os << "    <peptide p:peptide_id=\"" << sh.pPSM->getPeptideSequence() << "\"";
  if (Scores::getPrintDecoysInXml()) {
    if (sh.isDecoy())
      os << " p:decoy=\"true\"";
    else 
      os << " p:decoy=\"false\"";
  }
  os << ">" << endl;
  
  os << "      <svm_score>" << fixed       << sh.score     << "</svm_score>" << endl;
  os << "      <q_value>"   << scientific  << sh.q   << "</q_value>" << endl;
  os << "      <pep>"        << scientific  << sh.pep << "</pep>" << endl;
  
  if (Scores::getShowExpMass()) {
    os << "      <exp_mass>" << fixed << setprecision (4) << sh.pPSM->expMass << "</exp_mass>" << endl;
  }
  os << "      <calc_mass>" << fixed << setprecision (3)  << sh.pPSM->calcMass << "</calc_mass>" << endl;
  
  std::set<std::string>::const_iterator pidIt = sh.pPSM->proteinIds.begin();
  for ( ; pidIt != sh.pPSM->proteinIds.end() ; ++pidIt) {
    os << "      <protein_id>" << getRidOfUnprintablesAndUnicode(*pidIt) << "</protein_id>" << endl;
  }
  
  os << "      <p_value>" << scientific << sh.p << "</p_value>" <<endl;
  os << "      <psm_ids>" << endl;
  
  // output all psms that contain the peptide
  std::vector<std::string>::const_iterator psmIt = sh.psms_list.begin();
  for ( ; psmIt != sh.psms_list.end() ; ++psmIt) {
    os << "        <psm_id>" << *psmIt << "</psm_id>" << endl;
  }
  os << "      </psm_ids>" << endl;
  os << "    </peptide>" << endl;
  
  return os;
}

bool Scores::usePi0_ = false;
bool Scores::printDecoysInXml_ = false;
bool Scores::showExpMass_ = false;
unsigned long Scores::seed_ = 1;

void Scores::merge(vector<Scores>& sv, double fdr) {
  scores_.clear();
  for (vector<Scores>::iterator a = sv.begin(); a != sv.end(); a++) {
    sort(a->begin(), a->end(), greater<ScoreHolder> ());
    if (usePi0_) a->estimatePi0();
    a->calcQ(fdr);
    a->normalizeScores(fdr);
    copy(a->begin(), a->end(), back_inserter(scores_));
  }
  sort(scores_.begin(), scores_.end(), greater<ScoreHolder> ());
  totalNumberOfDecoys_ = count_if(scores_.begin(),
      scores_.end(),
      mem_fun_ref(&ScoreHolder::isDecoy));
  totalNumberOfTargets_ = count_if(scores_.begin(),
      scores_.end(),
      mem_fun_ref(&ScoreHolder::isTarget));
  targetDecoySizeRatio_ = totalNumberOfTargets_ / max(1.0, (double)totalNumberOfDecoys_);
  if (usePi0_) estimatePi0();
  else pi0_ = 1.0;
}

void Scores::printRetentionTime(ostream& outs, double fdr) {
  std::vector<ScoreHolder>::iterator scoreIt = scores_.begin();
  for ( ; scoreIt != scores_.end(); ++scoreIt) {
    if (scoreIt->isTarget()) 
      outs << PSMDescription::unnormalize(scoreIt->pPSM->retentionTime) << "\t"
        << PSMDescription::unnormalize(doc_.estimateRT(scoreIt->pPSM->retentionFeatures))
        << "\t" << scoreIt->pPSM->peptide << endl;
  }
}

double Scores::calcScore(const double* feat) const {
  register int ix = FeatureNames::getNumFeatures();
  register double score = svmWeights_[ix];
  for (; ix--;) {
    score += feat[ix] * svmWeights_[ix];
  }
  return score;
}

/**
 * Returns the ScoreHolder object associated with a list of features
 * @param d array of features
 * @return pointer to ScoreHolder
 */
ScoreHolder* Scores::getScoreHolder(const double* d) {
  if (scoreMap_.size() == 0) {
    std::vector<ScoreHolder>::iterator scoreIt = scores_.begin();
    for ( ; scoreIt != scores_.end(); ++scoreIt) {
      scoreMap_[scoreIt->pPSM->features] = &(*scoreIt);
    }
  }
  std::map<const double*, ScoreHolder*>::iterator res = scoreMap_.find(d);
  if (res != scoreMap_.end()) {
    return res->second;
  }
  return NULL;
}

void Scores::fillFeatures(SetHandler& setHandler) {
  scores_.clear();
  setHandler.fillFeatures(scores_,1);
  setHandler.fillFeatures(scores_,-1);
  totalNumberOfTargets_ = setHandler.getSizeFromLabel(1);
  totalNumberOfDecoys_ = setHandler.getSizeFromLabel(-1);
  targetDecoySizeRatio_ = (double)totalNumberOfTargets_ / totalNumberOfDecoys_;
}

// Park–Miller random number generator
// from wikipedia
unsigned long Scores::lcg_rand() {
  //uint64_t
  seed_ = (seed_ * 279470273u) % 4294967291u;
  return seed_;
}

/**
 * Divides the PSMs from pin file into xval_fold cross-validation sets based on
 * their spectrum scan number
 * @param train vector containing the training sets of PSMs
 * @param test vector containing the test sets of PSMs
 * @param xval_fold: number of folds in train and test
 */
void Scores::createXvalSetsBySpectrum(std::vector<Scores>& train, 
    std::vector<Scores>& test, const unsigned int xval_fold) {
  // set the number of cross validation folds for train and test to xval_fold
  train.resize(xval_fold);
  test.resize(xval_fold);
  // remain keeps track of residual space available in each fold
  std::vector<int> remain(xval_fold);
  // set values for remain: initially each fold is assigned (tot number of
  // scores_ / tot number of folds)
  int fold = xval_fold, ix = scores_.size();
  while (fold--) {
    remain[fold] = ix / (fold + 1);
    ix -= remain[fold];
  }

  // store possible spectra with relative scores_
  multimap<unsigned int,ScoreHolder> spectraScores;
  // populate spectraScores
  std::vector<ScoreHolder>::iterator scoreIt = scores_.begin();
  for ( ; scoreIt != scores_.end(); ++scoreIt) {
    spectraScores.insert(pair<unsigned int,ScoreHolder>(scoreIt->pPSM->scan, *scoreIt));
  }

  // put scores into the folds; choose a fold (at random) and change it only
  // when scores from a new spectra are encountered
  // note: this works because multimap is an ordered container!
  unsigned int previousSpectrum = spectraScores.begin()->first;
  size_t randIndex = lcg_rand() % xval_fold;
  for (multimap<unsigned int, ScoreHolder>::iterator it = spectraScores.begin(); 
        it != spectraScores.end(); ++it) {
    const unsigned int curScan = (*it).first;
    const ScoreHolder sh = (*it).second;
    // if current score is from a different spectra than the one encountered in
    // the previous iteration, choose new fold
    
    if (previousSpectrum != curScan) {
      randIndex = lcg_rand() % xval_fold;
      // allow only indexes of folds that are non-full
      while (remain[randIndex] <= 0){
        randIndex = lcg_rand() % xval_fold;
      }
    }
    // insert
    for (unsigned int i = 0; i < xval_fold; i++) {
      if (i == randIndex) {
        test[i].addScoreHolder(sh);
      } else {
        train[i].addScoreHolder(sh);
      }
    }
    // update number of free position for used fold
    --remain[randIndex];
    // set previous spectrum to current one for next iteration
    previousSpectrum = curScan;
  }

  // calculate ratios of target over decoy for train and test set
  for (unsigned int i = 0; i < xval_fold; i++) {
    train[i].recalculateSizes();
    test[i].recalculateSizes();
  }
}

void Scores::recalculateSizes() {
  totalNumberOfTargets_ = 0;
  totalNumberOfDecoys_ = 0;
  std::vector<ScoreHolder>::const_iterator scoreIt = scores_.begin();
  for ( ; scoreIt != scores_.end(); ++scoreIt) {
    if (scoreIt->isTarget()) {
      ++totalNumberOfTargets_;
    } else {
      ++totalNumberOfDecoys_;
    }
  }
  targetDecoySizeRatio_ = totalNumberOfTargets_ / (double)totalNumberOfDecoys_;
}
    

void Scores::normalizeScores(double fdr) {
  // sets q=fdr to 0 and the median decoy to -1, linear transform the rest to fit
  unsigned int medianIndex = std::max(0u,totalNumberOfDecoys_/2u),decoys=0u;
  vector<ScoreHolder>::iterator it = scores_.begin();
  double q1 = it->score;
  double median = q1 + 1.0;

  for (; it != scores_.end(); ++it) {
    if (it->q < fdr)
      q1 = it->score;
    if (it->isDecoy()) {
      if(++decoys==medianIndex) {
        median = it->score;
        break;
      }
    }
  }
  //NOTE perhaps I should also check when q1 and median are both negatives
  //NOTE in such cases the normalization could give negative scores_ which would
  //     cause an assertion to fail in qvality
  if (q1 <= median || it == scores_.end()) {
    ostringstream temp;
    temp << "Error : the input data has too good separation between target "
         << "and decoy PSMs.\n" << std::endl;
    throw MyException(temp.str());
  }
   
  double diff = q1-median;
  std::vector<ScoreHolder>::iterator scoreIt = scores_.begin();
  for ( ; scoreIt != scores_.end(); ++scoreIt) {
    scoreIt->score -= q1;
    scoreIt->score /= diff;
    if (scoreIt->score <= 0 && VERB > 3) { // Why do we warn for this, it happens for most of the data
      std::cerr << "\nWARNING the score of the PSM " << scoreIt->pPSM->id << 
          " is less or equal than zero after normalization.\n" << std::endl;
    }
  }
  
}

/**
 * Calculates the SVM cost/score of each PSM and sorts them
 * @param w normal vector used for SVM cost
 * @param fdr FDR threshold specified by user (default 0.01)
 * @return number of true positives
 */
int Scores::calcScores(vector<double>& w, double fdr) {
  svmWeights_ = w;
  unsigned int ix;
  std::vector<ScoreHolder>::iterator scoreIt = scores_.begin();
  for ( ; scoreIt != scores_.end(); ++scoreIt) {
    scoreIt->score = calcScore(scoreIt->pPSM->features);
  }
  sort(scores_.begin(), scores_.end(), greater<ScoreHolder> ());
  if (VERB > 3) {
    if (scores_.size() >= 10) {
      cerr << "10 best scores and labels" << endl;
      for (ix = 0; ix < 10; ix++) {
        cerr << scores_[ix].score << " " << scores_[ix].label << endl;
      }
      cerr << "10 worst scores and labels" << endl;
      for (ix = scores_.size() - 10; ix < scores_.size(); ix++) {
        cerr << scores_[ix].score << " " << scores_[ix].label << endl;
      }
    } else {
      cerr << "Too few scores to display top and bottom PSMs (" << scores_.size() << " scores_ found)." << endl;
    }
  }
  return calcQ(fdr);
}

/**
 * Calculates the q-value for each psm in scores_: the q-value is the minimal
 * FDR of any set that includes the particular psm
 * @param fdr FDR threshold specified by user (default 0.01)
 * @return number of true positives
 */
int Scores::calcQ(double fdr) {
  assert(totalNumberOfDecoys_+totalNumberOfTargets_==size());

  int targets = 0, decoys = 0;
  double efp = 0.0, q; // estimated false positives, q value
  
  int numPos = 0;
  
  // NOTE check this
  std::vector<ScoreHolder>::iterator scoreIt = scores_.begin();
  for ( ; scoreIt != scores_.end(); ++scoreIt) {
    if (scoreIt->isTarget()) {
      targets++;
      scoreIt->p = (decoys+(double)1)/(totalNumberOfDecoys_+(double)1);
    } else {
      decoys++;
      efp = pi0_ * decoys * targetDecoySizeRatio_;
      scoreIt->p = (decoys)/(double)(totalNumberOfDecoys_);
    }
    if (targets) {
      q = efp / (double)targets;
    } else {
      q = pi0_;
    }
    if (q > pi0_) {
      q = pi0_;
    }
    scoreIt->q = q;
    if (fdr >= q) {
      numPos = targets;
    }
  }
  if (scores_.size() > 0) {
    for (int ix = scores_.size(); --ix;) {
      if (scores_[ix - 1].q > scores_[ix].q) {
        scores_[ix - 1].q = scores_[ix].q;
      }
    }
  }
  return numPos;
}

void Scores::generateNegativeTrainingSet(AlgIn& data, const double cneg) {
  unsigned int ix2 = 0;
  std::vector<ScoreHolder>::const_iterator scoreIt = scores_.begin();
  for ( ; scoreIt != scores_.end(); ++scoreIt) {
    if (scoreIt->isDecoy()) {
      data.vals[ix2] = scoreIt->pPSM->features;
      data.Y[ix2] = -1;
      data.C[ix2++] = cneg;
    }
  }
  data.negatives = ix2;
}

void Scores::generatePositiveTrainingSet(AlgIn& data, const double fdr,
    const double cpos) {
  unsigned int ix2 = data.negatives, p = 0;
  std::vector<ScoreHolder>::const_iterator scoreIt = scores_.begin();
  for ( ; scoreIt != scores_.end(); ++scoreIt) {
    if (scoreIt->isTarget()) {
      if (fdr < scoreIt->q) {
        break;
      }
      data.vals[ix2] = scoreIt->pPSM->features;
      data.Y[ix2] = 1;
      data.C[ix2++] = cpos;
      ++p;
    }
  }
  data.positives = p;
  data.m = ix2;
}

/**
 * Routine that sees to that only unique peptides are kept (used for analysis
 * on peptide-fdr rather than psm-fdr)
 */
void Scores::weedOutRedundant() {
  
  // lexicographically order the scores_ (based on peptides names,labels and scores)
  std::sort(scores_.begin(), scores_.end(), lexicOrderProb());

  /*
  * much faster and simpler version but it does not fill up psms_list     
  * which will simply iterate over the unique peptides and the removed list many times 
  * scores_.erase(std::unique(scores_.begin(), scores_.end(), mycmp), scores_.end());
  */

  //NOTE the weed out PSMs might not be cleaned at the end

  std::vector<ScoreHolder> uniquePeptideScores;
  std::string previousPeptide = "";
  int previousLabel = 0;
  // run a pointer down the scores_ list
  std::vector<ScoreHolder>::iterator scoreIt = scores_.begin();
  for ( ; scoreIt != scores_.end(); scoreIt++){
    // compare pointer's peptide with previousPeptide
    std::string currentPeptide = scoreIt->pPSM->getPeptideSequence();
    if (currentPeptide != previousPeptide || scoreIt->label != previousLabel) {
      // insert as a new score
      uniquePeptideScores.push_back(*scoreIt);
      // update previousPeptide
      previousPeptide = currentPeptide;
      previousLabel = scoreIt->label;
    }
    // append the psm_id
    uniquePeptideScores.back().psms_list.push_back(scoreIt->pPSM->id);
  }

  scores_ = uniquePeptideScores;
  sort(scores_.begin(), scores_.end(), greater<ScoreHolder> ());

  totalNumberOfDecoys_ = count_if(scores_.begin(), scores_.end(),
                                 mem_fun_ref(&ScoreHolder::isDecoy));
  totalNumberOfTargets_ = count_if(scores_.begin(), scores_.end(),
                                  mem_fun_ref(&ScoreHolder::isTarget));
  targetDecoySizeRatio_ = 
      totalNumberOfTargets_ / max(1.0, (double)totalNumberOfDecoys_);

  if (usePi0_) estimatePi0();
  else pi0_ = 1.0;
}

/**
 * Routine that sees to that only unique spectra are kept for TDC
 */
void Scores::weedOutRedundantTDC() {
  // order the scores_ (based on spectra id and scores_)
  std::sort(scores_.begin(), scores_.end(), OrderScanMassCharge());

  /*
  * much faster and simpler version but it does not fill up psms_list     
  * which will simply iterate over the unique peptides and the removed list many times 
  * scores_.erase(std::unique(scores_.begin(), scores_.end(), mycmp), scores_.end());
  */

  vector<ScoreHolder> uniquePSMs = vector<ScoreHolder>();
  unsigned previousSpectra = 0;
  double previousExpMass = 0.0;
  //int previousLabel;
  // run a pointer down the scores_ list
  vector<ScoreHolder>::iterator current = scores_.begin();
  for (;current!=scores_.end(); current++){
    // compare pointer's spectra with previous spectra
    unsigned currentSpectra = current->pPSM->scan;
    double currentExpMass = current->pPSM->expMass;
    //int currentLabel = current->label;
    if (currentSpectra != previousSpectra || previousExpMass != currentExpMass) {
     uniquePSMs.push_back(*current);
     previousSpectra = currentSpectra;
     previousExpMass = currentExpMass;
    }
  }
  scores_ = uniquePSMs;
  sort(scores_.begin(), scores_.end(), greater<ScoreHolder> ());
  totalNumberOfDecoys_ = count_if(scores_.begin(),
    scores_.end(),
    mem_fun_ref(&ScoreHolder::isDecoy));
  totalNumberOfTargets_ = count_if(scores_.begin(),
    scores_.end(),
    mem_fun_ref(&ScoreHolder::isTarget));
  targetDecoySizeRatio_ = totalNumberOfTargets_ / max(1.0, (double)totalNumberOfDecoys_);
  
  if (usePi0_) estimatePi0();
  else pi0_ = 1.0;
}

void Scores::recalculateDescriptionOfCorrect(const double fdr) {
  doc_.clear();
  std::vector<ScoreHolder>::const_iterator scoreIt = scores_.begin();
  for ( ; scoreIt != scores_.end(); ++scoreIt) {
    if (scoreIt->isTarget()) {
      //      if (fdr>scores_[ix1].pPSM->q) {
      if (0.0 >= scoreIt->q) {
        doc_.registerCorrect(*(scoreIt->pPSM));
      }
    }
  }
  doc_.trainCorrect();
  setDOCFeatures();
}

void Scores::setDOCFeatures() {
  std::vector<ScoreHolder>::const_iterator scoreIt = scores_.begin();
  for ( ; scoreIt != scores_.end(); ++scoreIt) {
    doc_.setFeatures(*(scoreIt->pPSM));
  }
}

int Scores::getInitDirection(const double fdr, vector<double>& direction) {
  int bestPositives = -1;
  int bestFeature = -1;
  bool lowBest = false;
  for (unsigned int featNo = 0; featNo < FeatureNames::getNumFeatures(); featNo++) {
    for (std::vector<ScoreHolder>::iterator scoreIt = scores_.begin(); 
         scoreIt != scores_.end(); ++scoreIt) {
      scoreIt->score = scoreIt->pPSM->features[featNo];
    }
    sort(scores_.begin(), scores_.end());
    // check once in forward direction (high scores are good) and once in backward
    for (int i = 0; i < 2; i++) {
      int positives = 0, decoys = 0;
      double efp = 0.0, q;
      std::vector<ScoreHolder>::const_iterator scoreIt = scores_.begin();
      for ( ; scoreIt != scores_.end(); ++scoreIt) {
        if (scoreIt->isTarget()) {
          positives++;
        } else {
          decoys++;
          efp = pi0_ * decoys * targetDecoySizeRatio_;
        }
        if (positives) {
          q = efp / (double)positives;
        } else {
          q = pi0_;
        }
        if (fdr <= q) {
          if (positives > bestPositives && scores_.begin()->score != scoreIt->score) {
            bestPositives = positives;
            bestFeature = featNo;
            lowBest = (i == 0);
          }
          if (i == 0) {
            reverse(scores_.begin(), scores_.end());
          }
          break;
        }
      }
    }
  }
  for (int ix = FeatureNames::getNumFeatures(); ix--;) {
    direction[ix] = 0;
  }
  direction[bestFeature] = (lowBest ? -1 : 1);
  if (VERB > 1) {
    cerr << "Selected feature number " << bestFeature + 1
        << " as initial search direction, could separate "
        << bestPositives << " positives in that direction" << endl;
  }
  return bestPositives;
}

void Scores::estimatePi0() {
  vector<pair<double, bool> > combined;
  vector<double> pvals;
  transform(scores_.begin(), scores_.end(), back_inserter(combined),
            mem_fun_ref(&ScoreHolder::toPair));
  // Estimate pi0_
  PosteriorEstimator::getPValues(combined, pvals);
  pi0_ = PosteriorEstimator::estimatePi0(pvals);
}

void Scores::calcPep() {
  vector<pair<double, bool> > combined;
  transform(scores_.begin(),
      scores_.end(),
      back_inserter(combined),
      mem_fun_ref(&ScoreHolder::toPair));
  vector<double> peps;
  // Logistic regression on the data
  PosteriorEstimator::estimatePEP(combined, pi0_, peps, true);
  for (size_t ix = 0; ix < scores_.size(); ix++) {
    (scores_[ix]).pep = peps[ix];
  }
}

unsigned Scores::getQvaluesBelowLevel(double level) {
  unsigned hits = 0;
  std::vector<ScoreHolder>::const_iterator scoreIt = scores_.begin();
  for ( ; scoreIt != scores_.end(); ++scoreIt) {
    if (scoreIt->isTarget() && scoreIt->q < level) {
      hits++;
    }
  }
  return hits;
}
