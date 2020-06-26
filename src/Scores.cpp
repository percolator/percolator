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

#include "Singleton.hpp"
#include "LayerOrderedHeap.hpp"
#include "Clock.hpp"
#include "DataSet.h"
#include "Normalizer.h"
#include "SetHandler.h"
#include "Scores.h"
#include "Globals.h"
#include "PosteriorEstimator.h"
#include "ssl.h"
#include "MassHandler.h"

#ifdef CRUX
#include "app/PercolatorAdapter.h"
#endif

inline bool operator>(const ScoreHolder& one, const ScoreHolder& other) {
  return (one.score > other.score) 
      || (one.score == other.score && one.pPSM->scan > other.pPSM->scan) 
      || (one.score == other.score && one.pPSM->scan == other.pPSM->scan && 
            one.pPSM->expMass > other.pPSM->expMass)
      || (one.score == other.score && one.pPSM->scan == other.pPSM->scan && 
            one.pPSM->expMass == other.pPSM->expMass && one.label > other.label);
}

inline bool operator<(const ScoreHolder& one, const ScoreHolder& other) {
  return (one.score < other.score) 
      || (one.score == other.score && one.pPSM->scan < other.pPSM->scan) 
      || (one.score == other.score && one.pPSM->scan == other.pPSM->scan && 
            one.pPSM->expMass < other.pPSM->expMass)
      || (one.score == other.score && one.pPSM->scan == other.pPSM->scan && 
            one.pPSM->expMass == other.pPSM->expMass && one.label < other.label);
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

void ScoreHolder::printPSM(ostream& os, bool printDecoys, bool printExpMass) {
  if (!isDecoy() || printDecoys) {
    os << "    <psm p:psm_id=\"" << pPSM->getId() << "\"";
    if (printDecoys) {
      if (isDecoy())
        os << " p:decoy=\"true\"";
      else 
        os << " p:decoy=\"false\"";
    }
    os << ">" << endl;
    
    os << "      <svm_score>" << fixed      << score << "</svm_score>" << endl;
    os << "      <q_value>"   << scientific << q     << "</q_value>" << endl;
    os << "      <pep>"       << scientific << pep   << "</pep>" << endl;
    
    if (printExpMass) {
      os << "      <exp_mass>" << fixed << setprecision (4) << pPSM->expMass << "</exp_mass>" << endl;
    }   
    
    os << "      <calc_mass>" << fixed << setprecision (3) << pPSM->calcMass << "</calc_mass>" << endl;
    
    if (DataSet::getCalcDoc()) {
      os << "      <retentionTime observed=\"" 
         << pPSM->getUnnormalizedRetentionTime()
         << "\" predicted=\""
         << PSMDescriptionDOC::unnormalize(pPSM->getPredictedRetentionTime()) << "\"/>"
         << endl;
    }

    if (pPSM->getPeptideSequence().size() > 0) {
      string n = pPSM->getFlankN();
      string c = pPSM->getFlankC();
      string centpep = pPSM->getPeptideSequence();
      os << "      <peptide_seq n=\"" << n << "\" c=\"" << c << "\" seq=\"" << centpep << "\"/>" << endl;
    }
    
    std::vector<std::string>::const_iterator pidIt = pPSM->proteinIds.begin();
    for ( ; pidIt != pPSM->proteinIds.end() ; ++pidIt) {
      os << "      <protein_id>" << getRidOfUnprintablesAndUnicode(*pidIt) << "</protein_id>" << endl;
    }
    
    os << "      <p_value>" << scientific << p << "</p_value>" <<endl;
    os << "    </psm>" << endl;
  }
}

void ScoreHolder::printPeptide(ostream& os, bool printDecoys, bool printExpMass, Scores& fullset) {
  if (!isDecoy() || printDecoys) {  
    os << "    <peptide p:peptide_id=\"" << pPSM->getPeptideSequence() << "\"";
    if (printDecoys) {
      if (isDecoy())
        os << " p:decoy=\"true\"";
      else 
        os << " p:decoy=\"false\"";
    }
    os << ">" << endl;
    
    os << "      <svm_score>" << fixed       << score     << "</svm_score>" << endl;
    os << "      <q_value>"   << scientific  << q   << "</q_value>" << endl;
    os << "      <pep>"        << scientific  << pep << "</pep>" << endl;
    
    if (printExpMass) {
      os << "      <exp_mass>" << fixed << setprecision (4) << pPSM->expMass << "</exp_mass>" << endl;
    }
    os << "      <calc_mass>" << fixed << setprecision (3)  << pPSM->calcMass << "</calc_mass>" << endl;
    
    std::vector<std::string>::const_iterator pidIt = pPSM->proteinIds.begin();
    for ( ; pidIt != pPSM->proteinIds.end() ; ++pidIt) {
      os << "      <protein_id>" << getRidOfUnprintablesAndUnicode(*pidIt) << "</protein_id>" << endl;
    }
    
    os << "      <p_value>" << scientific << p << "</p_value>" <<endl;
    os << "      <psm_ids>" << endl;
    
    // output all psms that contain the peptide
    std::vector<PSMDescription*>::const_iterator psmIt = fullset.getPsms(pPSM).begin();
    for ( ; psmIt != fullset.getPsms(pPSM).end() ; ++psmIt) {
      os << "        <psm_id>" << (*psmIt)->getId() << "</psm_id>" << endl;
    }
    os << "      </psm_ids>" << endl;
    os << "    </peptide>" << endl;
  }
}

void Scores::merge(std::vector<Scores>& sv, double fdr, bool skipNormalizeScores) {
  scores_.clear();
  for (std::vector<Scores>::iterator a = sv.begin(); a != sv.end(); a++) {
    sort(a->begin(), a->end(), greater<ScoreHolder> ());
    a->checkSeparationAndSetPi0();
    a->calcQ(fdr);
    if (!skipNormalizeScores) {
      a->normalizeScores(fdr);
    }
    copy(a->begin(), a->end(), back_inserter(scores_));
  }
  postMergeStep();
}

void Scores::postMergeStep() {
  sort(scores_.begin(), scores_.end(), greater<ScoreHolder> ());
  totalNumberOfDecoys_ = count_if(scores_.begin(),
      scores_.end(),
      mem_fun_ref(&ScoreHolder::isDecoy));
  totalNumberOfTargets_ = count_if(scores_.begin(),
      scores_.end(),
      mem_fun_ref(&ScoreHolder::isTarget));
  targetDecoySizeRatio_ = totalNumberOfTargets_ / max(1.0, (double)totalNumberOfDecoys_);
  checkSeparationAndSetPi0();
}

void Scores::printRetentionTime(ostream& outs, double fdr) {
  std::vector<ScoreHolder>::iterator scoreIt = scores_.begin();
  for ( ; scoreIt != scores_.end(); ++scoreIt) {
    if (scoreIt->isTarget()) 
      outs << scoreIt->pPSM->getUnnormalizedRetentionTime() << "\t"
        << PSMDescriptionDOC::unnormalize(doc_.estimateRT(scoreIt->pPSM->getRetentionFeatures()))
        << "\t" << scoreIt->pPSM->peptide << endl;
  }
}

double Scores::calcScore(const double* feat, const std::vector<double>& w) const {
  register int ix = FeatureNames::getNumFeatures();
  register double score = w[ix];
  for (; ix--;) {
    score += feat[ix] * w[ix];
  }
  return score;
}

void Scores::scoreAndAddPSM(ScoreHolder& sh, 
    const std::vector<double>& rawWeights, FeatureMemoryPool& featurePool) {
  const unsigned int numFeatures = FeatureNames::getNumFeatures();
  if (DataSet::getCalcDoc()) {
    size_t numRTFeatures = RTModel::totalNumRTFeatures();
    double* rtFeatures = new double[numRTFeatures]();
    DescriptionOfCorrect::calcRegressionFeature(sh.pPSM);
    for (size_t i = 0; i < numRTFeatures; ++i) {
      rtFeatures[i] = Normalizer::getNormalizer()->normalize(rtFeatures[i], numFeatures + i);
    }
    sh.pPSM->setRetentionFeatures(rtFeatures);
    doc_.setFeatures(sh.pPSM);
  }
  
  for (unsigned int j = 0; j < numFeatures; j++) {
    sh.score += sh.pPSM->features[j] * rawWeights[j];
  }
  sh.score += rawWeights[numFeatures];
  
  featurePool.deallocate(sh.pPSM->features);
  sh.pPSM->deleteRetentionFeatures();
  
  if (sh.label == 1) {
    ++totalNumberOfTargets_;
  } else if (sh.label == -1) {
    ++totalNumberOfDecoys_;
  }
  
  if (sh.label != 1 && sh.label != -1) {
    std::cerr << "Warning: the PSM " << sh.pPSM->getId()
        << " has a label not in {1,-1} and will be ignored." << std::endl;
    PSMDescription::deletePtr(sh.pPSM);
  } else {
    scores_.push_back(sh);
  }
}

void Scores::print(int label, std::ostream& os) {
#ifndef CRUX
  std::vector<ScoreHolder>::iterator scoreIt = scores_.begin();
  os << "PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds\n";
  for ( ; scoreIt != scores_.end(); ++scoreIt) {
    if (scoreIt->label == label) {
      std::ostringstream out;
      scoreIt->pPSM->printProteins(out);
      ResultHolder rh(scoreIt->score, scoreIt->q, scoreIt->pep, scoreIt->pPSM->getId(), scoreIt->pPSM->peptide, out.str());
      os << rh << std::endl;
    }
  }
#else
  PercolatorAdapter::printScores(this, label, os);
#endif
}

void Scores::populateWithPSMs(SetHandler& setHandler) {
  scores_.clear();
  setHandler.populateScoresWithPSMs(scores_, 1);
  setHandler.populateScoresWithPSMs(scores_, -1);
  totalNumberOfTargets_ = setHandler.getSizeFromLabel(1);
  totalNumberOfDecoys_ = setHandler.getSizeFromLabel(-1);
  targetDecoySizeRatio_ = (double)totalNumberOfTargets_ / totalNumberOfDecoys_;
  
  if (VERB > 1) {
    cerr << "Train/test set contains " << totalNumberOfTargets_
        << " positives and " << totalNumberOfDecoys_
        << " negatives, size ratio=" << targetDecoySizeRatio_
        << " and pi0=" << pi0_ << endl;
  }
  
  if (totalNumberOfTargets_ == 0) {
    ostringstream oss;
    oss << "Error: no target PSMs were provided.\n";
    if (NO_TERMINATE) {
      cerr << oss.str() << "No-terminate flag set: ignoring error." << std::endl;
    } else {
      throw MyException(oss.str());
    }
  }
  
  if (totalNumberOfDecoys_ == 0) {
    ostringstream oss;
    oss << "Error: no decoy PSMs were provided.\n";
    if (NO_TERMINATE) {
      cerr << oss.str() << "No-terminate flag set: ignoring error." << std::endl;
    } else {
      throw MyException(oss.str());
    }
  }
  
  // check for the minimum recommended number of positive and negative hits
  if (totalNumberOfTargets_ <= (unsigned)(FeatureNames::getNumFeatures() * 5)) {
    std::cerr << "Warning : the number of positive samples read is too small to perform a correct classification.\n" << std::endl;
  }
  if (totalNumberOfDecoys_ <= (unsigned)(FeatureNames::getNumFeatures() * 5)) {
    std::cerr << "Warning : the number of negative samples read is too small to perform a correct classification.\n" << std::endl;
  }
}

/**
 * Divides the PSMs from pin file into xval_fold cross-validation sets based on
 * their spectrum scan number
 * @param train vector containing the training sets of PSMs
 * @param test vector containing the test sets of PSMs
 * @param xval_fold: number of folds in train and test
 */
void Scores::createXvalSetsBySpectrum(std::vector<Scores>& train, 
    std::vector<Scores>& test, const unsigned int xval_fold, 
    FeatureMemoryPool& featurePool) {
  // set the number of cross validation folds for train and test to xval_fold
  train.resize(xval_fold, Scores(usePi0_));
  test.resize(xval_fold, Scores(usePi0_));
  // remain keeps track of residual space available in each fold
  std::vector<int> remain(xval_fold);
  // set values for remain: initially each fold is assigned (tot number of
  // scores_ / tot number of folds)
  int fold = xval_fold, ix = scores_.size();
  while (fold--) {
    remain[fold] = ix / (fold + 1);
    ix -= remain[fold];
  }
  
  std::sort(scores_.begin(), scores_.end(), OrderScanMassCharge());
  
  // put scores into the folds; choose a fold (at random) and change it only
  // when scores from a new spectra are encountered
  unsigned int previousSpectrum = scores_.begin()->pPSM->scan;
  size_t randIndex = PseudoRandom::lcg_rand() % xval_fold;
  for (std::vector<ScoreHolder>::iterator it = scores_.begin(); 
        it != scores_.end(); ++it) {
    const unsigned int curScan = (*it).pPSM->scan;
    const ScoreHolder sh = (*it);
    // if current score is from a different spectra than the one encountered in
    // the previous iteration, choose new fold
    
    if (previousSpectrum != curScan) {
      randIndex = PseudoRandom::lcg_rand() % xval_fold;
      // allow only indexes of folds that are non-full
      while (remain[randIndex] <= 0){
        randIndex = PseudoRandom::lcg_rand() % xval_fold;
      }
    }
    // insert
    for (unsigned int i = 0; i < xval_fold; ++i) {
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
  for (unsigned int i = 0; i < xval_fold; ++i) {
    train[i].recalculateSizes();
    test[i].recalculateSizes();
  }
  
  if (featurePool.isInitialized()) {
    boost::unordered_map<double*, double*> movedAddresses;
    size_t idx = 0;
    for (unsigned int i = 0; i < xval_fold; ++i) {
      bool isTarget = true;
      test[i].reorderFeatureRows(featurePool, isTarget, movedAddresses, idx);
      isTarget = false;
      test[i].reorderFeatureRows(featurePool, isTarget, movedAddresses, idx);
    }
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

void Scores::reorderFeatureRows(FeatureMemoryPool& featurePool, 
    bool isTarget, boost::unordered_map<double*, double*>& movedAddresses, size_t& idx) {
  size_t numFeatures = FeatureNames::getNumFeatures();
  std::vector<ScoreHolder>::const_iterator scoreIt = scores_.begin();
  for ( ; scoreIt != scores_.end(); ++scoreIt) {
    if (scoreIt->isTarget() == isTarget) {
      double* newAddress = featurePool.addressFromIdx(idx++);
      double* oldAddress = scoreIt->pPSM->features;
      while (movedAddresses.find(oldAddress) != movedAddresses.end()) {
        oldAddress = movedAddresses[oldAddress];
      }
      std::swap_ranges(oldAddress, oldAddress + numFeatures, newAddress);
      scoreIt->pPSM->features = newAddress;
      if (oldAddress != newAddress) {
        movedAddresses[newAddress] = oldAddress;
      }
    }
  }
}

void Scores::set_total_number_of_decoys() {
  if (total_number_of_decoys > 0)
    return;

  for (unsigned long i=0; i < scores_.size(); ++i)
    if (scores_[i].isDecoy())
      ++total_number_of_decoys;
}

double* Scores::get_vector_of_just_decoy_scores() {
  std::cout << "get vector of just ... " << std::endl;
  set_total_number_of_decoys();
  double* decoys = new double[total_number_of_decoys];
  unsigned long decoy_index = 0;
  for (unsigned long i=0; i < scores_.size(); ++i) 
    if (scores_[i].isDecoy()) {
      decoys[decoy_index] = scores_[i].score;
      ++decoy_index;
    }

  std::cout << "returning get vector of just ... " << std::endl;  
  return decoys;
}

// sets q=fdr to 0 and the median decoy to -1, linear transform the rest to fit
void Scores::normalizeScores(double fdr) {
  std::cout << "enter normalize scores " << std::endl;
  unsigned int medianIndex = std::max(0u,totalNumberOfDecoys_/2u),decoys=0u;
  std::vector<ScoreHolder>::iterator it = scores_.begin();
  // fixme:remove double fdrScore = it->score;
  double fdrScore = scores_[largest_index_lt_fdr-1].score; // fixme: check inclusive/exclusive
  double medianDecoyScore = fdrScore + 1.0;
  // fixme: call loh if fdr > fdr seen so far

  // nth element to find median ( array with just decoys ) with new, delete, std::copy
  double* decoy_scores = get_vector_of_just_decoy_scores();
  const unsigned position_of_median = total_number_of_decoys/2;
  std::nth_element(decoy_scores, decoy_scores + (position_of_median), decoy_scores + total_number_of_decoys);
  medianDecoyScore = decoy_scores[position_of_median];
  delete[] decoy_scores;
  std::cout << "after del decoy scores" << std::endl;

  //NOTE perhaps I should also check when fdrScore and medianDecoyScore are both 
  //  negative. In such cases the normalization could give negative scores which 
  //  would cause an assertion to fail in qvality
  
  double diff = fdrScore - medianDecoyScore;
  std::vector<ScoreHolder>::iterator scoreIt = scores_.begin();
  for ( ; scoreIt != scores_.end(); ++scoreIt) {
    scoreIt->score -= fdrScore;
    if (diff > 0.0) {
      scoreIt->score /= diff;
    }
  }
  std::cout << "exit normalize scores " << std::endl;
}

/**
 * Calculates the SVM cost/score of each PSM and sorts them
 * @param w normal vector used for SVM cost
 * @param fdr FDR threshold specified by user (default 0.01)
 * @return number of true positives
 */

int Scores::calcScoresSorted(std::vector<double>& w, double fdr, bool skipDecoysPlusOne, bool print_scores) {
  unsigned int ix;
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
      cerr << "Too few scores to display top and bottom PSMs (" << scores_.size() << " scores found)." << endl;
    }
  }
  return calcQ(fdr, skipDecoysPlusOne, print_scores);
}

double Scores::get_fdr(unsigned tps, unsigned fps) {
  double max_bias_correction = fps*(1.0-pi0_);
  return (double(fps+max_bias_correction)) / tps; //fixme: / max(1,tps)
}

int Scores::calcScoresLOH(const double fdr_threshold, const pair<double, bool> *const orig_combined_begin, pair<double, bool> *combined_begin, pair<double, bool> *combined_end, int num_tps_at_start_of_layer, int num_fps_at_start_of_layer, LayerArithmetic* la) {
  // For starters just have pi0=1
  unsigned long n = combined_end - combined_begin;
  LayerOrderedHeap<pair<double, bool>> loh(combined_begin, n, la);

  // Iterate through layers in loh
  for (int i=loh.n_layers()-1; i>=0; --i) {
    pair<double, bool>* layer_begin = loh.layer_begin(i);
    pair<double, bool>* layer_end = loh.layer_end(i);

    unsigned layer_size = layer_end - layer_begin;
    
    while ( i > 0 && *(layer_begin-1) == layer_begin[0]) {
      std::cout << "there was  a tie, moved from" << layer_begin;
      --i;
      layer_begin = loh.layer_begin(i);
      std::cout << " to " << layer_begin << std::endl;
      layer_size = layer_end - layer_begin;
    }

    // count number of targets and decoys in layer to calculate best and worst possible q values
    unsigned num_tps_in_layer = 0;
    for (auto layer_iter = layer_begin; layer_iter < layer_end; ++layer_iter)
      num_tps_in_layer += layer_iter->second;

    unsigned num_fps_in_layer = (layer_end-layer_begin) - num_tps_in_layer;

    num_tps_at_start_of_layer -= num_tps_in_layer;
    num_fps_at_start_of_layer -= num_fps_in_layer;      

    //fixme: can do something special with homogeneous layers (BOTH IN TARGET/DECOY AND IN TIED SCORES)
    //if (num_tps_in_layer == 0 || num_fps_in_layer == 0)
    //  continue;

    // optimistic fdr
    double fdr_optimistic = get_fdr(num_tps_at_start_of_layer + num_tps_in_layer, num_fps_at_start_of_layer);

    // pessemistic fdr
    double fdr_pessemistic = get_fdr(num_tps_at_start_of_layer, num_fps_at_start_of_layer + num_fps_in_layer);

    // fdr_threshold is between the highest and lowest q values, then recurse
    if (fdr_optimistic <= fdr_threshold && fdr_threshold <= fdr_pessemistic) {
      if (n > 1) {
	std::cout << "Recurse at " << layer_begin-orig_combined_begin << std::endl;
        int result = calcScoresLOH(fdr_threshold, orig_combined_begin, layer_begin, layer_end, num_tps_at_start_of_layer+num_tps_in_layer, num_fps_at_start_of_layer + num_fps_in_layer, la);
	if (result != -1)
	  return result;
      }
      else {
	return num_tps_at_start_of_layer + num_tps_in_layer+num_fps_at_start_of_layer + num_fps_in_layer;
      }
    }
  }
  return -1;
}

int Scores::calcScores(std::vector<double>& w, double fdr, bool skipDecoysPlusOne) {
  std::vector<ScoreHolder>::iterator scoreIt = scores_.begin();
  for ( ; scoreIt != scores_.end(); ++scoreIt) {
    scoreIt->score = calcScore(scoreIt->pPSM->features, w);
  }

  constexpr bool print_time = false;
  constexpr bool USE_SORT =false;
  constexpr bool USE_LOH = true;

  int loh_score, sort_score;
  std::vector<pair<double, bool> > combined;

  if (USE_LOH) {
    // else use layer-ordered heap
    Clock c;
    getScoreLabelPairs(combined);
    constexpr double ALPHA = 2;
    Singleton<LayerArithmetic, double> s_la;
    LayerArithmetic*la = &s_la.get_instance(ALPHA);

    // todo: is there a count of all decoys and all targets identified
    // in data? we may not need this O(n) pass.
    unsigned long total_num_tps = 0;

    for (unsigned long i=0; i<combined.size(); ++i, ++scoreIt) {
      total_num_tps += combined[i].second;
    }

    unsigned long total_num_fps = combined.size() - total_num_tps;

    if (skipDecoysPlusOne) 
      loh_score = calcScoresLOH(fdr, &combined[0], &combined[0], &combined[0]+combined.size(),total_num_tps,total_num_fps, la);
    else
      // fixme: should we add +1 below? compare against their results
      loh_score = calcScoresLOH(fdr, &combined[0], &combined[0], &combined[0]+combined.size(),total_num_tps,total_num_fps+1,la);

    highest_fdr_calculated = std::max(highest_fdr_calculated, fdr);
    fdr_has_been_calculated = true;
    largest_index_lt_fdr = std::max(largest_index_lt_fdr, loh_score); // fixme: check inclusive or exclusive

    double loh_time = c.tock();

    double eps=0.000001;
    unsigned long i=0;

    std::vector<ScoreHolder>::iterator scoreIt = scores_.begin();
    for (; scoreIt != scores_.end(); ++scoreIt, ++i) {
      if (i < loh_score)
	scoreIt->q = fdr-0.001;
      else
	scoreIt->q = fdr+0.001;
    }

    

    return loh_score;
  }

  if (USE_SORT) {
    Clock c;
    sort_score = calcScoresSorted(w,fdr,skipDecoysPlusOne,false);
    double sort_time = c.tock();

    return sort_score;
  }
}

void Scores::getScoreLabelPairs(std::vector<pair<double, bool> >& combined) {
  combined.clear();
  transform(scores_.begin(), scores_.end(), back_inserter(combined),
            mem_fun_ref(&ScoreHolder::toPair));
}

/**
 * Calculates the q-value for each psm in scores_: the q-value is the minimal
 * FDR of any set that includes the particular psm
 * @param fdr FDR threshold specified by user (default 0.01)
 * @return number of true positives
 */
int Scores::calcQ(double fdr, bool skipDecoysPlusOne, bool print_scores) {
  assert(totalNumberOfDecoys_+totalNumberOfTargets_==size());
  
  std::vector<pair<double, bool> > combined;
  getScoreLabelPairs(combined);
  
  std::vector<double> qvals;
  PosteriorEstimator::setNegative(true); // also get q-values for decoys
  PosteriorEstimator::getQValues(pi0_, combined, qvals, skipDecoysPlusOne);
  if (print_scores) {
    std::cout << "Q VALUES: ";
    for (int i=0; i < qvals.size(); ++i)
      std::cout << qvals[i] << ", ";
    std::cout << std::endl;    
  }

  // set q-values and count number of positives
  std::vector<double>::const_iterator qIt = qvals.begin();
  std::vector<ScoreHolder>::iterator scoreIt = scores_.begin();
  
  int numPos = 0;
  for (; qIt != qvals.end(); ++qIt, ++scoreIt) {
    scoreIt->q = *qIt;
    if (scoreIt->q < fdr && scoreIt->isTarget()) ++numPos;
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
    const double cpos, const bool trainBestPositive) {
  unsigned int ix2 = data.negatives, p = 0;
  
  std::vector<ScoreHolder>::iterator lastUniqueIt = scores_.end();
  if (trainBestPositive) {
    std::sort(scores_.begin(), scores_.end(), OrderScanLabel());
    lastUniqueIt = std::unique(scores_.begin(), scores_.end(), UniqueScanLabel());
    std::sort(scores_.begin(), lastUniqueIt, greater<ScoreHolder> ());
  }

  std::vector<ScoreHolder>::const_iterator scoreIt = scores_.begin();
  unsigned long i=0;
  for (; i<largest_index_lt_fdr; ++i) {
    if (scores_[i].isTarget()) {
        data.vals[ix2] = scores_[i].pPSM->features;
        data.Y[ix2] = 1;
        data.C[ix2++] = cpos;
        ++p;
    }
  }
  // fixme: remove once working
  // for ( ; scoreIt != lastUniqueIt; ++scoreIt) {
  //   if (scoreIt->isTarget()) {
  //     if (scoreIt->q <= fdr) {
  //       data.vals[ix2] = scoreIt->pPSM->features;
  //       data.Y[ix2] = 1;
  //       data.C[ix2++] = cpos;
  //       ++p;
  //     }
  //   }
  // }
  data.positives = p;
  data.m = ix2;
}

void Scores::weedOutRedundant() {
  std::map<std::string, unsigned int> peptideSpecCounts;
  double specCountQvalThreshold = -1.0;
  weedOutRedundant(peptideSpecCounts, specCountQvalThreshold);
}

/**
 * Routine that sees to that only unique peptides are kept (used for analysis
 * on peptide-fdr rather than psm-fdr)
 */
void Scores::weedOutRedundant(std::map<std::string, unsigned int>& peptideSpecCounts, double specCountQvalThreshold) {
  // lexicographically order the scores_ (based on peptides names,labels and scores)
  std::sort(scores_.begin(), scores_.end(), lexicOrderProb());
  
  /*
  * much simpler version but it does not fill up the peptide-PSM map:
  * scores_.erase(std::unique(scores_.begin(), scores_.end(), mycmp), scores_.end());
  */
  
  std::string previousPeptide = "";
  int previousLabel = 0;
  size_t lastWrittenIdx = 0u;
  for (size_t idx = 0u; idx < scores_.size(); ++idx){
    std::string currentPeptide = scores_.at(idx).pPSM->getPeptideSequence();
    int currentLabel = scores_.at(idx).label;
    if (currentPeptide != previousPeptide || currentLabel != previousLabel) {
      // insert as a new score
      scores_.at(lastWrittenIdx++) = scores_.at(idx);
      previousPeptide = currentPeptide;
      previousLabel = currentLabel;
    }
    // append the psm
    peptidePsmMap_[scores_.at(lastWrittenIdx - 1).pPSM].push_back(scores_.at(idx).pPSM);
    if (specCountQvalThreshold > 0.0 && scores_.at(idx).q < specCountQvalThreshold) {
      ++peptideSpecCounts[currentPeptide];
    }
  }
  scores_.resize(lastWrittenIdx);
  postMergeStep();
}

/**
 * Routine that sees to that only unique spectra are kept for TDC
 */
void Scores::weedOutRedundantTDC() {
  // order the scores (based on spectra id and score)
  std::sort(scores_.begin(), scores_.end(), OrderScanMassCharge());  
  scores_.erase(std::unique(scores_.begin(), scores_.end(), UniqueScanMassCharge()), scores_.end());
  
  /* does not actually release memory because of memory fragmentation
  double previousExpMass = 0.0;
  unsigned int previousScan = 0u;
  size_t lastWrittenIdx = 0u;
  for (size_t idx = 0u; idx < scores_.size(); ++idx){
    double currentExpMass = scores_.at(idx).pPSM->expMass;
    int currentScan = scores_.at(idx).pPSM->scan;
    if (currentExpMass != previousExpMass || currentScan != previousScan) {
      // insert as a new score
      scores_.at(lastWrittenIdx++).swap(scores_.at(idx));
      previousScan = currentScan;
      previousExpMass = currentExpMass;
    } else {
      PSMDescription::deletePtr(scores_.at(idx).pPSM);
    }
  }
  scores_.resize(lastWrittenIdx);
  */
  postMergeStep();
}

/**
 * Routine that sees to that only 1 target and 1 decoy spectra are kept for 
 * mix-max when using multiple hits per spectrum and separate searches
 */
void Scores::weedOutRedundantMixMax() {
  // order the scores (based on spectra id and score)
  std::sort(scores_.begin(), scores_.end(), OrderScanMassLabelCharge());  
  scores_.erase(std::unique(scores_.begin(), scores_.end(), UniqueScanMassLabelCharge()), scores_.end());
  
  postMergeStep();
}

void Scores::recalculateDescriptionOfCorrect(const double fdr) {
  doc_.clear();

  std::vector<ScoreHolder>::const_iterator scoreIt = scores_.begin();
  for (unsigned long i=0; i<largest_index_lt_fdr; ++i) 
    if (scores_[i].isTarget())
      doc_.registerCorrect(scores_[i].pPSM);
  
  doc_.trainCorrect();
}

void Scores::setDOCFeatures(Normalizer* pNorm) {
  std::vector<ScoreHolder>::const_iterator scoreIt = scores_.begin();
  for ( ; scoreIt != scores_.end(); ++scoreIt) {
    doc_.setFeaturesNormalized(scoreIt->pPSM, pNorm);
  }
}

int Scores::getInitDirection(const double initialSelectionFdr, std::vector<double>& direction) {
  int bestPositives = -1;
  int bestFeature = -1;
  bool lowBest = false;
  
  // for determining the initial direction, the decoys+1 in the FDR estimates 
  // is too restrictive for small datasets
  bool skipDecoysPlusOne = true; 
  
  for (unsigned int featNo = 0; featNo < FeatureNames::getNumFeatures(); featNo++) {
    for (std::vector<ScoreHolder>::iterator scoreIt = scores_.begin(); 
         scoreIt != scores_.end(); ++scoreIt) {
      scoreIt->score = scoreIt->pPSM->features[featNo];
    }
    sort(scores_.begin(), scores_.end());
    // check once in forward direction (i = 0, higher scores are better) and 
    // once in backward direction (i = 1, lower scores are better)
    for (int i = 0; i < 2; i++) {
      if (i == 1) {
        reverse(scores_.begin(), scores_.end());
      }
      int positives = calcQ(initialSelectionFdr, skipDecoysPlusOne);
      if (positives > bestPositives) {
        bestPositives = positives;
        bestFeature = featNo;
        lowBest = (i == 0);
      }
    }
  }
  for (int ix = FeatureNames::getNumFeatures(); ix--;) {
    direction[ix] = 0;
  }
  
  if (bestPositives <= 0) {
    ostringstream oss;
    oss << "Error in the input data: cannot find an initial direction with " 
        << "positive training examples. "
        << "Consider setting/raising the initial training FDR threshold (--train-initial-fdr)." << std::endl;
    if (NO_TERMINATE) {
      cerr << oss.str();
      std::cerr << "No-terminate flag set: setting initial direction to the "
                << "first feature and ignoring the error." << std::endl;
      bestFeature = 0;
    } else {
      throw MyException(oss.str() + "Terminating.\n");
    }
  }
  
  if (bestFeature >= 0) {
    direction[bestFeature] = (lowBest ? -1 : 1);
  }
  
  if (VERB > 1) {
    cerr << "Selected feature " << bestFeature + 1
        << " as initial direction. Could separate "
        << bestPositives << " training set positives with q<" << initialSelectionFdr << " in that direction." << endl;
  }
  return bestPositives;
}

void Scores::checkSeparationAndSetPi0() {
  std::vector<pair<double, bool> > combined;
  getScoreLabelPairs(combined);
  
  std::vector<double> pvals;
  PosteriorEstimator::getPValues(combined, pvals);
  
  pi0_ = 1.0;
  bool tooGoodSeparation = PosteriorEstimator::checkSeparation(pvals);
  if (tooGoodSeparation) {
    ostringstream oss;
    oss << "Error in the input data: too good separation between target "
        << "and decoy PSMs.\n";
    if (NO_TERMINATE) {
      cerr << oss.str();
      if (usePi0_) {
        std::cerr << "No-terminate flag set: setting pi0 = 1 and ignoring error." << std::endl;
      } else {
        std::cerr << "No-terminate flag set: ignoring error." << std::endl;
      }
    } else {
      throw MyException(oss.str() + "Terminating.\n");
    }
  } else if (usePi0_) {
    pi0_ = PosteriorEstimator::estimatePi0(pvals);
  }
}

void Scores::calcPep() {
  std::vector<pair<double, bool> > combined;
  getScoreLabelPairs(combined);
  
  std::vector<double> peps;
  // Logistic regression on the data
  PosteriorEstimator::estimatePEP(combined, usePi0_, pi0_, peps, true);
  for (size_t ix = 0; ix < scores_.size(); ix++) {
    scores_[ix].pep = peps[ix];
  }
}

unsigned Scores::getQvaluesBelowLevel(double level) {
  unsigned hits = 0;

  for (unsigned long i=0; i<largest_index_lt_fdr; ++i)
    if (scores_[i].isTarget())
      ++hits;
  return hits;
}
