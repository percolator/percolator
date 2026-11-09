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

const double ProteinProbEstimator::prior_protein = 0.5;
bool ProteinProbEstimator::calcProteinLevelProb = false;

ProteinProbEstimator::ProteinProbEstimator(bool trivialGrouping, double absenceRatio, 
					     bool outputEmpirQVal, std::string decoyPattern, 
					     double specCountQvalThreshold) : 
	  trivialGrouping_(trivialGrouping), absenceRatio_(absenceRatio), pi0_(1.0),
	  numberDecoyProteins_(0u), numberTargetProteins_(0u), usePi0_(true), 
	  outputEmpirQVal_(outputEmpirQVal), decoyPattern_(decoyPattern), fdr_(1.0),
	  peptideSpecCounts_(), specCountQvalThreshold_(specCountQvalThreshold) {}

ProteinProbEstimator::~ProteinProbEstimator() {}

bool ProteinProbEstimator::initialize(Scores& peptideScores, const Enzyme* enzyme, std::string& protEstimatorDecoyPrefix) {
  decoyPattern_ = protEstimatorDecoyPrefix;
  setTargetandDecoysNames(peptideScores);
  return true;
}

void ProteinProbEstimator::computeStatistics() {
  std::sort(proteins_.begin(), proteins_.end(), IntCmpScore());
  
  estimatePValues();
  
  /** compute q values from PEPs, does not require pi0 **/
  estimateQValues();
  
  if (usePi0_ && outputEmpirQVal_) {
    pi0_ = estimatePi0();
    if (pi0_ <= 0.0 || pi0_ > 1.0) pi0_ = proteins_.rbegin()->getQ();
    if (VERB > 1) {
      std::cerr << "protein pi0 estimate = " << pi0_ << std::endl;
    }
  }
  
  estimateQValuesEmp();

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
  if (!proteinFN.empty() || !proteinDecoyFN.empty()) {
    if (!proteinFN.empty()) {
      ofstream proteinOut(proteinFN.data(), ios::out);
      print(proteinOut,false);
      proteinOut.close();	
    }
    if (!proteinDecoyFN.empty()) {
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
  std::vector<ProteinScoreHolder>::const_iterator it = proteins_.begin();
  bool isTarget = false;
  for (; it != proteins_.end(); it++) {
    double score = it->getScore();
    isTarget = isTarget || it->isTarget();
    if (lastProteinInGroup(it)) {
      combined.push_back(std::make_pair(score, isTarget));
      isTarget = false;
    }
  }
}

void ProteinProbEstimator::estimatePValues() {
  // assuming combined sorted in best hit first order
  std::vector<std::pair<double, bool> > combined;
  getCombinedList(combined);
  
  std::vector<double> pvalues;
  PosteriorEstimator::getPValues(combined, pvalues);
  
  std::vector<ProteinScoreHolder>::iterator protIt = proteins_.begin();
  std::vector<double>::const_iterator pIt = pvalues.begin();
  bool isTarget = false;
  for (; protIt != proteins_.end(); ++protIt) {
    protIt->setP(*pIt);
    isTarget = isTarget || protIt->isTarget();
    if (lastProteinInGroup(protIt) && isTarget && pIt+1 != pvalues.end()) {
      ++pIt;
      isTarget = false;
    }
  }
}

double ProteinProbEstimator::estimatePi0(const unsigned int numBoot) {
  std::vector<double> pvalues;
  std::transform(proteins_.begin(), proteins_.end(), 
    std::back_inserter(pvalues), std::mem_fn(&ProteinScoreHolder::getP));
  
  double pi0 = 1.0;
  bool tooGoodSeparation = PosteriorEstimator::checkSeparation(pvalues);
  if (tooGoodSeparation) {
    ostringstream oss;
    oss << "Error in the input data: too good separation between target "
        << "and decoy proteins.\n";
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
    pi0 = PosteriorEstimator::estimatePi0(pvalues);
  }
  return pi0;
}

unsigned ProteinProbEstimator::getQvaluesBelowLevel(double level) {
  std::set<int> identifiedGroupIds;
  unsigned nP = 0;
  for (std::vector<ProteinScoreHolder>::const_iterator myP = proteins_.begin(); 
          myP != proteins_.end(); ++myP) {
    if (myP->getQemp() < level && myP->isTarget()) {
      nP++;
      identifiedGroupIds.insert(myP->getGroupId());
    }
  }
  
  if (trivialGrouping_) {
    return static_cast<unsigned>(identifiedGroupIds.size());
  } else {
    return nP;
  }
}



void ProteinProbEstimator::estimateQValues() {
  std::vector<double> peps; // Posterior Error Probabilities
  
  std::vector<ProteinScoreHolder>::const_iterator it = proteins_.begin();
  bool isTarget = false;
  for (; it != proteins_.end(); it++) {
    isTarget = isTarget || it->isTarget();
    if (lastProteinInGroup(it) && (countDecoyQvalue_ || isTarget)) {
      peps.push_back(it->getPEP());
      isTarget = false;
    }
  }
  
  std::vector<double> qvaluesTmp;
  PosteriorEstimator::getQValuesFromPEP(peps, qvaluesTmp);
  
  std::vector<ProteinScoreHolder>::iterator protIt = proteins_.begin();
  std::vector<double>::const_iterator qIt = qvaluesTmp.begin();
  isTarget = false;
  for (; protIt != proteins_.end(); ++protIt) {
    isTarget = isTarget || protIt->isTarget();
    protIt->setQ(*qIt);
    if (lastProteinInGroup(protIt) && (countDecoyQvalue_ || isTarget) 
                                   && qIt+1 != qvaluesTmp.end()) {
      ++qIt;
    }
  }
}

void ProteinProbEstimator::estimateQValuesEmp() {
  std::vector<pair<double, bool> > combined;
  getCombinedList(combined);
  
  std::vector<double> qvaluesEmp;
  PosteriorEstimator::setNegative(true); // also get q-values for decoys
  PosteriorEstimator::getQValues(pi0_, combined, qvaluesEmp);
  
  std::vector<ProteinScoreHolder>::iterator protIt = proteins_.begin();
  std::vector<double>::const_iterator qIt = qvaluesEmp.begin();
  for (; protIt != proteins_.end(); ++protIt) {
    protIt->setQemp(*qIt);
    if (lastProteinInGroup(protIt) && qIt+1 != qvaluesEmp.end()) {
      ++qIt;
    }
  }
}

void ProteinProbEstimator::setTargetandDecoysNames(Scores& peptideScores) {
  int numGroups = 0;
  bool decoyFound = false;
  std::vector<ScoreHolder>::iterator psm = peptideScores.begin();
  for (; psm!= peptideScores.end(); ++psm) {
    // for each protein
    PSMDescription::ProteinList proteinIds = psm->pPSM->proteins();
    PSMDescription::ProteinList::const_iterator protIt = proteinIds.begin();
    for (; protIt != proteinIds.end(); ++protIt) {
      ProteinScoreHolder::Peptide peptide(psm->pPSM->getPeptideSequence(), 
          psm->isDecoy(), 0.0, psm->pep, psm->q, psm->score);
      if (proteinToIdxMap_.find(*protIt) == proteinToIdxMap_.end()) {
	      ProteinScoreHolder newProtein(*protIt, psm->isDecoy(), peptide, ++numGroups);
	      proteinToIdxMap_[*protIt] = proteins_.size();
	      proteins_.push_back(newProtein);
	      
	      if (!useDecoyPrefix) {
	        if (psm->isDecoy()) {
	          falsePosSet_.insert(*protIt);
	          decoyFound = true;
	        } else {
	          truePosSet_.insert(*protIt);
	        }
	      } else if (isDecoy(*protIt)) {
	        decoyFound = true;
	      }
      } else {
      	proteins_.at(proteinToIdxMap_[*protIt]).addPeptide(peptide);
      }
    }
  }
  if (!useDecoyPrefix) {
    numberDecoyProteins_ = static_cast<unsigned int>(falsePosSet_.size());
    numberTargetProteins_ = static_cast<unsigned int>(truePosSet_.size());
  }
  
  if (!decoyFound) {
    std::cerr << "Warning: No decoy proteins found. "
      << "Check if the correct decoy pattern was specified "
      << "with the -P flag;\n  if the target protein is called \"protA\" and the "
      << "decoy protein \"decoy_protA\", use the option \"-P decoy_\"." 
      << std::endl << std::endl;
  }
}

void ProteinProbEstimator::addSpectralCounts(Scores& peptideScores) {
  std::vector<ScoreHolder>::iterator psm = peptideScores.begin();
  for (; psm!= peptideScores.end(); ++psm) {
    // for each protein
    PSMDescription::ProteinList proteinIds = psm->pPSM->proteins();
    PSMDescription::ProteinList::const_iterator protIt = proteinIds.begin();
    std::set<unsigned int> seenProteinIdxs;
    for (; protIt != proteinIds.end(); ++protIt) {
      if (proteinToIdxMap_.find(*protIt) != proteinToIdxMap_.end()) {
        unsigned int proteinIdx = static_cast<unsigned int>(proteinToIdxMap_[*protIt]);
        if (seenProteinIdxs.find(proteinIdx) == seenProteinIdxs.end()) {
          seenProteinIdxs.insert(proteinIdx);
        }
      }
    }
    
    bool isUnique = (seenProteinIdxs.size() == 1);
    unsigned int psmCount = peptideSpecCounts_[psm->pPSM->getPeptideSequence()];
    std::set<unsigned int>::const_iterator protIdxIt = seenProteinIdxs.begin();
    for (; protIdxIt != seenProteinIdxs.end(); ++protIdxIt) {
      proteins_[*protIdxIt].addSpecCounts(psmCount, isUnique);
    }
  }
}

unsigned ProteinProbEstimator::countTargets(
    const std::vector<std::string>& proteinList) {
  unsigned count = 0;
  std::vector<std::string>::const_iterator it = proteinList.begin();
  for (; it != proteinList.end(); it++) {
    if (useDecoyPrefix) {
      if ((*it).find(decoyPattern_) == std::string::npos) {
      	count++;
      }
    } else {
      if (truePosSet_.count(*it) != 0) {
      	count++;
      }
    }
  }
  return count;
}

unsigned ProteinProbEstimator::countDecoys(
    const std::vector<std::string>& proteinList) {
  unsigned count = 0;
  std::vector<std::string>::const_iterator it = proteinList.begin();
  for (; it != proteinList.end(); it++) {
    if (useDecoyPrefix) {
      if((*it).find(decoyPattern_) != std::string::npos) {
      	count++;
      }
    } else {
      if(falsePosSet_.count(*it) != 0) {
      	count++;
      }
    }
  }
  return count;
}

bool ProteinProbEstimator::isDecoy(const std::string& proteinName) {
  //NOTE faster with decoyPrefix but we have to assume that the label that 
  // identifies decoys is in decoyPattern_
  return (bool)(useDecoyPrefix ? proteinName.find(decoyPattern_) 
	    != std::string::npos : falsePosSet_.count(proteinName) != 0);
}
    
bool ProteinProbEstimator::isTarget(const std::string& proteinName) {
  //NOTE faster with decoyPrefix but we have to assume that the label that 
  // identifies decoys is in decoyPattern_
  return (bool)(useDecoyPrefix ? proteinName.find(decoyPattern_) 
		  == std::string::npos : truePosSet_.count(proteinName) != 0);
}


void ProteinProbEstimator::print(ostream& myout, bool decoy) {  
  myout << "ProteinId\tProteinGroupId\tq-value\tposterior_error_prob\t";
  if (specCountQvalThreshold_ > 0.0) {
    myout << "spec_count_unique\tspec_count_all\t";
  }
  myout << "peptideIds" << std::endl;
  
  for (std::vector<ProteinScoreHolder>::const_iterator myP = proteins_.begin(); 
	        myP != proteins_.end(); myP++) {
    if( (decoy && myP->isDecoy()) || (!decoy && myP->isTarget())) {
      myout << myP->getName() << "\t" << myP->getGroupId() << "\t" 
            << myP->getQemp() << "\t" << myP->getPEP() << "\t";
      if (specCountQvalThreshold_ > 0.0) {
        myout << myP->getSpecCountsUnique() << "\t" << myP->getSpecCountsAll() << "\t";
      }
      const std::vector<ProteinScoreHolder::Peptide>& peptides = myP->getPeptidesByRef();
      std::vector<ProteinScoreHolder::Peptide>::const_iterator peptIt = peptides.begin();
      for(; peptIt != peptides.end(); peptIt++) {
        if (peptIt->name != "") {
          myout << peptIt->name << " ";
        }
      }
      myout << '\n';
    }
  }
}
