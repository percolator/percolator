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

#include "ProteinProbEstimator.h"

const double ProteinProbEstimator::target_decoy_ratio = 1.0;
const double ProteinProbEstimator::psmThresholdMayu = 0.90;
const double ProteinProbEstimator::prior_protein = 0.5;
bool ProteinProbEstimator::calcProteinLevelProb = false;
/** Helper functions **/

template<class T> 
void bootstrap(const vector<T>& in, vector<T>& out, size_t max_size = 1000) {
  out.clear();
  double n = in.size();
  size_t num_draw = min(in.size(), max_size);
  for (size_t ix = 0; ix < num_draw; ++ix) {
    size_t draw = (size_t)((double)PseudoRandom::lcg_rand() / ((double)PseudoRandom::kRandMax + (double)1) * n);
    out.push_back(in[draw]);
  }
  // sort in desending order
  std::sort(out.begin(), out.end());
}


ProteinProbEstimator::ProteinProbEstimator(bool trivialGrouping, double absenceRatio, 
					     bool outputEmpirQVal, std::string decoyPattern) : 
	  trivialGrouping_(trivialGrouping), absenceRatio_(absenceRatio), pi0_(1.0),
	  numberDecoyProteins_(0u), numberTargetProteins_(0u), usePi0_(true), 
	  outputEmpirQVal_(outputEmpirQVal), decoyPattern_(decoyPattern), fdr_(1.0) {}

ProteinProbEstimator::~ProteinProbEstimator() {
  FreeAll(qvalues);
  FreeAll(qvaluesEmp);
  FreeAll(pvalues_);
  
  for (std::multimap<double, std::vector<std::string> >::iterator it = pepProteinMap_.begin();
        it != pepProteinMap_.end(); it++) {
	  FreeAll(it->second);
  }
      
  for (std::map<const std::string,ProteinScoreHolder*>::iterator it = proteins_.begin(); 
        it != proteins_.end(); it++) {
	  if (it->second)
	    delete it->second;
  }
}

bool ProteinProbEstimator::initialize(Scores& peptideScores) {
  setTargetandDecoysNames(peptideScores);
  return true;
}

/** use the Mayu approach to calculating FDRs **/
void ProteinProbEstimator::computeFDR() {
  if (VERB > 1)
    std::cerr << "Estimating Protein FDR ... " << std::endl;
    
  ProteinFDRestimator fastReader;
  fastReader.setDecoyPrefix(decoyPattern_);
  fastReader.setTargetDecoyRatio(target_decoy_ratio);
  fastReader.setEqualDeepthBinning(binning_equal_deepth);
  fastReader.setNumberBins(number_bins);
  fastReader.correctIdenticalSequences(targetProteins_, decoyProteins_);
  //These guys are the number of target and decoys proteins but from the subset of PSM with FDR < threshold
  std::set<std::string> numberTP;
  std::set<std::string> numberFP;
  getTPandPFfromPeptides(psmThresholdMayu,numberTP,numberFP);
    
  double fptol = fastReader.estimateFDR(numberTP,numberFP);
    
  if (fptol == -1) {
    fdr_ = 1.0;
    
    if(VERB > 1)
      std::cerr << "There was an error estimating the Protein FDR..\n" << std::endl;
  } else {	
    fdr_ = (fptol/(double)numberTP.size());
    if(fdr_ <= 0 || fdr_ >= 1.0) fdr_ = 1.0;      
    
    if(VERB > 1) {
      std::cerr << "Estimated Protein FDR at ( " << psmThresholdMayu << ") PSM FDR is : " 
      << fdr_ << " with " << fptol << " expected number of false positives proteins\n" << std::endl;
    }
  }
}

void ProteinProbEstimator::computeStatistics() {
  if (pvalues_.size() == 0) {
    estimatePValues();
  }
  
  /** compute q values from PEPs, does not require pi0 **/
  estimateQValues();
  
  if (usePi0_ && !mayufdr && outputEmpirQVal_) {
    pi0_ = estimatePi0();
    if (pi0_ <= 0.0 || pi0_ > 1.0) pi0_ = *qvalues.rbegin();
    if (VERB > 1) {
      std::cerr << "protein pi0 estimate = " << pi0_ << std::endl;
    }
  }
  
  estimateQValuesEmp();
  updateProteinProbabilities();

  if (VERB > 1) {
    if (trivialGrouping_) {
      std::cerr << "Number of protein groups identified at q-value = 0.01: ";
    } else {
      std::cerr << "Number of proteins identified at q-value = 0.01: ";
    }
    std::cerr << getQvaluesBelowLevel(0.01) << std::endl;
  }
}

void ProteinProbEstimator::printOut(const std::string &proteinFN, 
				      const std::string &proteinDecoyFN) {
  if(!proteinFN.empty() || !proteinDecoyFN.empty()) {
    if(!proteinFN.empty()) {
      ofstream proteinOut(proteinFN.data(), ios::out);
      print(proteinOut,false);
      proteinOut.close();	
    }
    if(!proteinDecoyFN.empty()) {
      ofstream proteinOut(proteinDecoyFN.data(), ios::out);
      print(proteinOut,true);
      proteinOut.close();
    }
  } else {
    print(std::cout);
  }
}

void ProteinProbEstimator::getCombinedList(
    std::vector<std::pair<double, bool> >& combined) {
  std::multimap<double, std::vector<std::string> >::const_iterator it = pepProteinMap_.begin();
  for (; it != pepProteinMap_.end(); it++) {
    double prob = it->first;
    if (trivialGrouping_) {
      bool isTarget = (countTargets(it->second) > 0);
      combined.push_back(std::make_pair(prob, isTarget));
    } else {
      std::vector<std::string> proteinList = it->second;
      for(std::vector<std::string>::const_iterator itP = proteinList.begin();
	          itP != proteinList.end(); itP++) {
        std::string proteinName = *itP;
        bool isDecoy = proteins_[proteinName]->getIsDecoy();
        combined.push_back(std::make_pair(prob, !isDecoy));
      }
    }
  }
}

void ProteinProbEstimator::estimatePValues() {
  // assuming combined sorted in best hit first order
  std::vector<std::pair<double, bool> > combined;
  getCombinedList(combined);
  
  pvalues_.clear();
  PosteriorEstimator::getPValues(combined, pvalues_);
}

void ProteinProbEstimator::getTPandPFfromPeptides(double psm_threshold, 
						  std::set<std::string> &numberTP, 
						  std::set<std::string> &numberFP) {
  /* The original paper of Mayu describes a protein as :
   * FP = if only if all its peptides with q <= threshold are decoy
   * TP = at least one of its peptides with q <= threshold is target
   * However, the way mayu estimates it on the program is like this :
   * FP = any protein that contains a decoy psm with q <= threshold
   * TP = any protein that contains a target psm with q <= threshold
   * They do not consider protein containing both decoy and target psms
   * Also, mayu estimates q as the empirical (target-decoy) q value.
   * Percolator estimates q as the empirical (target-decoy) q value and adjusted by pi0
   * Mayu extracts the list of TP and FP proteins from PSM level whereas percolator
   * extract the list of TP and FP proteins from peptide level, this avoids redundancy and
   * gives a better calibration since peptide level q values are re-adjusted in percolator.
   * This creates sometimes a difference in the number of TP and FP proteins between percolator and Mayus 
   * which causes a slight difference in the estimated protein FDR
   */
  for (std::map<std::string,ProteinScoreHolder*>::const_iterator it = proteins_.begin();
       it != proteins_.end(); it++) {
    unsigned num_target_confident = 0;
    unsigned num_decoy_confident = 0;
    std::string protname = it->first;
    std::vector<ProteinScoreHolder::Peptide*> peptides = it->second->getPeptides();
    for(std::vector<ProteinScoreHolder::Peptide*>::const_iterator itP = peptides.begin();
          itP != peptides.end(); itP++) {
      ProteinScoreHolder::Peptide *p = *itP;
      if(p->q <= psm_threshold && p->isdecoy)
        ++num_decoy_confident;
      if(p->q <= psm_threshold && !p->isdecoy)
        ++num_target_confident;
    }
    if (num_decoy_confident > 0) {
      numberFP.insert(protname);
    }
    if (num_target_confident > 0) {
      numberTP.insert(protname);
    }
  }
}

double ProteinProbEstimator::estimatePi0(const unsigned int numBoot) {
  return PosteriorEstimator::estimatePi0(pvalues_, numBoot);
}

unsigned ProteinProbEstimator::getQvaluesBelowLevel(double level) {
  std::set<int> identifiedGroupIds;
  unsigned nP = 0;
  for (std::map<const std::string, ProteinScoreHolder*>::const_iterator myP = proteins_.begin(); 
          myP != proteins_.end(); ++myP) {
    if (myP->second->getQemp() < level && !myP->second->getIsDecoy()) {
      nP++;
      identifiedGroupIds.insert(myP->second->getGroupId());
    }
  }
  
  if (trivialGrouping_) {
    return identifiedGroupIds.size();
  } else {
    return nP;
  }
}

unsigned ProteinProbEstimator::getQvaluesBelowLevelDecoy(double level) { 
  std::set<int> identifiedGroupIds;
  unsigned nP = 0;
  for (std::map<const std::string, ProteinScoreHolder*>::const_iterator myP = proteins_.begin(); 
          myP != proteins_.end(); ++myP) {
    if (myP->second->getQ() < level && myP->second->getIsDecoy()) {
      nP++;
      identifiedGroupIds.insert(myP->second->getGroupId());
    }
  }
  
  if (trivialGrouping_) {
    return identifiedGroupIds.size();
  } else {
    return nP;
  }
}


void ProteinProbEstimator::estimateQValues() {
  std::vector<double> peps; // Posterior Error Probabilities
  std::vector<size_t> idxs; // index of last element in peps vector for each protein group
  
  for (std::multimap<double,std::vector<std::string> >::const_iterator it = pepProteinMap_.begin(); 
          it != pepProteinMap_.end(); it++) {
    int ntargets = countTargets(it->second);
    int ndecoys = it->second.size() - ntargets;
    if (trivialGrouping_) {
      if (ntargets > 0) ntargets = 1;
      if (ndecoys > 0) ndecoys = 1;
    }
    
    // in case we want to count and use target and decoys proteins while estimating qvalue from PEP
    if (countDecoyQvalue_) {
      ntargets += ndecoys;
    }
    for (size_t i = 0; i < ntargets; ++i) {
      peps.push_back(it->first);
    }
    idxs.push_back(peps.size() - 1);
  }
  
  std::vector<double> qvaluesTmp;
  PosteriorEstimator::getQValuesFromPEP(peps, qvaluesTmp);
  
  qvalues.clear();
  for (std::vector<size_t>::const_iterator it = idxs.begin(); it != idxs.end(); ++it) {
    qvalues.push_back(qvaluesTmp.at(*it));
  }
}

void ProteinProbEstimator::estimateQValuesEmp() {
  std::vector<pair<double, bool> > combined;
  getCombinedList(combined);
  
  qvaluesEmp.clear();
  PosteriorEstimator::setNegative(true); // also get q-values for decoys
  PosteriorEstimator::getQValues(pi0_, combined, qvaluesEmp);
}

void ProteinProbEstimator::updateProteinProbabilities() {
  std::vector<double> peps; // posterior error probabilities, not peptides
  std::vector<std::vector<std::string> > proteinNames;
  std::transform(pepProteinMap_.begin(), pepProteinMap_.end(), std::back_inserter(peps), RetrieveKey());
  std::transform(pepProteinMap_.begin(), pepProteinMap_.end(), std::back_inserter(proteinNames), RetrieveValue());
  unsigned protIdx = 0;
  for (unsigned pepIdx = 0; pepIdx < peps.size(); pepIdx++) {
    double pep = peps[pepIdx];
    std::vector<std::string> proteinlist = proteinNames[pepIdx];
    for (unsigned j = 0; j < proteinlist.size(); j++) { 
      int protGroupId = protIdx + 1;
      size_t idx = protIdx;
      if (trivialGrouping_) {
        idx = pepIdx;
        protGroupId = pepIdx + 1;
      }
      
      std::string proteinName = proteinlist[j];
      proteins_[proteinName]->setPEP(pep);
      proteins_[proteinName]->setQ(qvalues[idx]);
      proteins_[proteinName]->setQemp(qvaluesEmp[idx]);
      proteins_[proteinName]->setP(pvalues_[idx]);
      proteins_[proteinName]->setGroupId(protGroupId);
      ++protIdx;
    }
  }

}

void ProteinProbEstimator::setTargetandDecoysNames(Scores& peptideScores) {
  unsigned int numGroups = 0;
  for (vector<ScoreHolder>::iterator psm = peptideScores.begin(); psm!= peptideScores.end(); ++psm) {
    // for each protein
    for (std::vector<std::string>::iterator protIt = psm->pPSM->proteinIds.begin(); protIt != psm->pPSM->proteinIds.end(); protIt++) {
      ProteinScoreHolder::Peptide *peptide = new ProteinScoreHolder::Peptide(
          psm->pPSM->getPeptideSequence(), psm->isDecoy(),
          psm->p, psm->pep, psm->q, psm->score);
      if (proteins_.find(*protIt) == proteins_.end()) {
	      ProteinScoreHolder *newprotein = new ProteinScoreHolder(*protIt, psm->isDecoy(), peptide, ++numGroups);
	      proteins_.insert(std::make_pair(*protIt,newprotein));
	
	      if (psm->isDecoy()) {
	        falsePosSet.insert(*protIt);
	      } else {
	        truePosSet.insert(*protIt);
	      }
      } else {
      	proteins_[*protIt]->addPeptide(peptide);
      }
    }
  }  
  numberDecoyProteins_ = falsePosSet.size();
  numberTargetProteins_ = truePosSet.size();
}

/** Used by XMLInterface to read in the proteins with its sequence and store them for the Mayu method **/
void ProteinProbEstimator::addProteinDb(bool isDecoy, std::string name, std::string sequence, double length) {
  if (isDecoy)
    decoyProteins_.insert(std::make_pair(name,std::make_pair(sequence,length)));
  else
    targetProteins_.insert(std::make_pair(name,std::make_pair(sequence,length)));
}

unsigned ProteinProbEstimator::countTargets(const std::vector<std::string> &proteinList) {
  unsigned count = 0;
  for (std::vector<std::string>::const_iterator it = proteinList.begin(); it != proteinList.end(); it++) {
    if (useDecoyPrefix) {
      if ((*it).find(decoyPattern_) == std::string::npos) {
      	count++;
      }
    } else {
      if (truePosSet.count(*it) != 0) {
      	count++;
      }
    }
  }
  return count;
}

unsigned ProteinProbEstimator::countDecoys(const std::vector<std::string> &proteinList) {
  unsigned count = 0;
  for(std::vector<std::string>::const_iterator it = proteinList.begin(); it != proteinList.end(); it++) {
    if (useDecoyPrefix) {
      if((*it).find(decoyPattern_) != std::string::npos) {
      	count++;
      }
    } else {
      if(falsePosSet.count(*it) != 0) {
      	count++;
      }
    }
  }
  return count;
}

bool ProteinProbEstimator::isDecoy(const std::string& proteinName) {
  //NOTE faster with decoyPrefix but I assume the label that identifies decoys is in decoyPattern_
   return (bool)(useDecoyPrefix ? proteinName.find(decoyPattern_) 
	    != std::string::npos : falsePosSet.count(proteinName) != 0);
}
    
bool ProteinProbEstimator::isTarget(const std::string& proteinName) {
  //NOTE faster with decoyPrefix but I assume the label that identifies decoys is in decoyPattern_
   return (bool)(useDecoyPrefix ? proteinName.find(decoyPattern_) 
		  == std::string::npos : truePosSet.count(proteinName) != 0);
}


void ProteinProbEstimator::writeOutputToXML(string xmlOutputFN, bool outputDecoys) {
  std::vector<std::pair<std::string,ProteinScoreHolder*> > myvec(proteins_.begin(), proteins_.end());
  std::sort(myvec.begin(), myvec.end(), IntCmpProb());

  ofstream os;
  os.open(xmlOutputFN.data(), ios::app);
  // append PROTEINs tag
  os << "  <proteins>" << endl;
  for (std::vector<std::pair<std::string,ProteinScoreHolder*> > ::const_iterator myP = myvec.begin(); 
	      myP != myvec.end(); myP++) {
    if ( (!outputDecoys && !myP->second->getIsDecoy()) || (outputDecoys)) {
      os << "    <protein p:protein_id=\"" << myP->second->getName() << "\"";
	    if (outputDecoys) {
	      if (myP->second->getIsDecoy()) 
	        os << " p:decoy=\"true\"";
	      else  
	        os << " p:decoy=\"false\"";
	    }
	    os << ">" << endl;
	    
	    os << "      <pep>" << scientific << myP->second->getPEP() << "</pep>" << endl;
	  
	    if (outputEmpirQVal_) {
	      os << "      <q_value_emp>" << scientific << myP->second->getQemp() << "</q_value_emp>\n";
	    }
	  
	    os << "      <q_value>" << scientific << myP->second->getQ() << "</q_value>\n";
	    
	    if (outputEmpirQVal_) {
	      os << "      <p_value>" << scientific << myP->second->getP() << "</p_value>\n";
	    }
	  
	    std::vector<ProteinScoreHolder::Peptide*> peptides = myP->second->getPeptides();
	    for (std::vector<ProteinScoreHolder::Peptide*>::const_iterator peptIt = peptides.begin(); 
	        peptIt != peptides.end(); peptIt++) {
	      if ((*peptIt)->name != "") {
	        os << "      <peptide_seq seq=\"" << (*peptIt)->name << "\"/>"<<endl;
	      }
	    }
	    os << "    </protein>" << endl;
	  }
  }
    
  os << "  </proteins>" << endl << endl;
  os.close();
}

void ProteinProbEstimator::print(ostream& myout, bool decoy) {
  
  std::vector<std::pair<std::string,ProteinScoreHolder*> > myvec(proteins_.begin(), proteins_.end());
  std::sort(myvec.begin(), myvec.end(), IntCmpProb());
  
  myout << "ProteinId\tProteinGroupId\tq-value\tposterior_error_prob\tpeptideIds" << std::endl;
  
  for (std::vector<std::pair<std::string, ProteinScoreHolder*> > ::const_iterator myP = myvec.begin(); 
	        myP != myvec.end(); myP++) {
    if( (decoy && myP->second->getIsDecoy()) || (!decoy && !myP->second->getIsDecoy())) {
      myout << myP->second->getName() << "\t" << myP->second->getGroupId() << "\t" 
            << myP->second->getQemp() << "\t" << myP->second->getPEP() << "\t";
      std::vector<ProteinScoreHolder::Peptide*> peptides = myP->second->getPeptides();
      for(std::vector<ProteinScoreHolder::Peptide*>::const_iterator peptIt = peptides.begin(); peptIt != peptides.end(); peptIt++) {
        if((*peptIt)->name != "") {
          myout << (*peptIt)->name << " ";
        }
      }
      myout << std::endl;
    }
  }
}
