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

#include "FisherInterface.h"

FisherInterface::FisherInterface(const std::string& fastaDatabase,
    double pvalueCutoff, bool reportFragmentProteins, bool reportDuplicateProteins,
    bool trivialGrouping, double absenceRatio, bool outputEmpirQval, 
    std::string& decoyPattern) :
      ProteinProbEstimator(trivialGrouping, absenceRatio, outputEmpirQval, decoyPattern),
      fastaProteinFN_(fastaDatabase), maxPeptidePval_(pvalueCutoff),
      reportFragmentProteins_(reportFragmentProteins),
      reportDuplicateProteins_(reportDuplicateProteins),
      protInferenceMethod_(BESTPEPT) {
  if (absenceRatio == 1.0) usePi0_ = false;
}

FisherInterface::~FisherInterface() {}

bool FisherInterface::initialize(Scores* fullset) {
  peptideScores_ = fullset;
  int min_peptide_length = 1000, max_peptide_length = 0;
  int max_miscleavages = 0, max_non_enzymatic_flanks = 0;
  int total_non_enzymatic_flanks = 0;
  
  ENZYME_T enzyme;
  if (Enzyme::getEnzymeType() == Enzyme::TRYPSIN) {
    enzyme = TRYPSIN;
  } else if (Enzyme::getEnzymeType() == Enzyme::TRYPSINP) {
    enzyme = TRYPSINP;
  } else if (Enzyme::getEnzymeType() == Enzyme::CHYMOTRYPSIN) {
    enzyme = CHYMOTRYPSIN;
  } else if (Enzyme::getEnzymeType() == Enzyme::ELASTASE) {
    enzyme = ELASTASE;
  } else if (Enzyme::getEnzymeType() == Enzyme::LYSN) {
    enzyme = LYSN;
  } else if (Enzyme::getEnzymeType() == Enzyme::LYSC) {
    enzyme = LYSC;
  } else if (Enzyme::getEnzymeType() == Enzyme::ARGC) {
    enzyme = ARGC;
  } else if (Enzyme::getEnzymeType() == Enzyme::ASPN) {
    enzyme = ASPN;
  } else if (Enzyme::getEnzymeType() == Enzyme::GLUC) {
    enzyme = GLUC;
  } else {
    // not supported yet: THERMOLYSIN, PROTEINASEK, PEPSIN
    std::cerr << "Warning: specified protease " << Enzyme::getStringEnzyme()
              << " currently not supported, using trypsin to identify"
              << " duplicate and fragment proteins" << std::endl;
    Enzyme::setEnzyme(Enzyme::TRYPSIN);
    enzyme = TRYPSIN;
  }
  
  for (std::vector<ScoreHolder>::iterator shIt = peptideScores_->begin(); 
           shIt != peptideScores_->end(); ++shIt) {    
    std::string peptideSequenceFlanked = shIt->pPSM->getFullPeptideSequence();
    
    peptideSequenceFlanked = PSMDescription::removePTMs(peptideSequenceFlanked);
    std::string peptideSequence = PSMDescription::removeFlanks(peptideSequenceFlanked);
    
    int peptide_length = peptideSequence.size();
    min_peptide_length = std::min(min_peptide_length, peptide_length);
    max_peptide_length = std::max(max_peptide_length, peptide_length);
    
    int miscleavages = Enzyme::countEnzymatic(peptideSequence);
    if (miscleavages > max_miscleavages && VERB > 1) {
      std::cerr << "Miscleavage detected: " << peptideSequenceFlanked << std::endl;
      max_miscleavages = std::max(max_miscleavages, miscleavages);
    }
    
    int non_enzymatic_flanks = 0;
    if (!Enzyme::isEnzymatic(peptideSequenceFlanked[0], 
                             peptideSequenceFlanked[2]) 
           && peptideSequenceFlanked.substr(0,1) != "M") { // ignore protein N-term methionine cleavage
      ++non_enzymatic_flanks;
    }
    if (!Enzyme::isEnzymatic(peptideSequenceFlanked[peptideSequenceFlanked.size() - 3],
                             peptideSequenceFlanked[peptideSequenceFlanked.size() - 1])) {
      ++non_enzymatic_flanks;
    }
    if (non_enzymatic_flanks > max_non_enzymatic_flanks && VERB > 1) {
      std::cerr << "Non enzymatic flank detected: " << peptideSequenceFlanked << std::endl;
      max_non_enzymatic_flanks = std::max(max_non_enzymatic_flanks, non_enzymatic_flanks);
    }
    total_non_enzymatic_flanks += non_enzymatic_flanks;
  }
  // if more than half of the flanks or non enzymatic, probably the protease is wrong
  if (total_non_enzymatic_flanks > peptideScores_->size()) {
    std::cerr << "Warning: more than half of the cleavage sites are non enzymatic, "
              << "please verify that the right protease was specified." << std::endl;
  }
  
  DIGEST_T digest;
  std::string digestString = "";
  if (max_non_enzymatic_flanks == 0) {
    digest = FULL_DIGEST;
    digestString = "full";
  } else if (max_non_enzymatic_flanks == 1) {
    digest = PARTIAL_DIGEST;
    digestString = "partial";
  } else {
    digest = NON_SPECIFIC_DIGEST;
    enzyme = NO_ENZYME;
    digestString = "non-specific";
  }
  if (VERB > 1) {
    std::cerr << "Protein digestion parameters for duplicate/fragment detection (detected from PSM input):\n"
              << " enzyme=" << Enzyme::getStringEnzyme() 
              << ", digestion=" << digestString
              << ", min-pept-length=" << min_peptide_length
              << ", max-pept-length=" << max_peptide_length
              << ", max-miscleavages=" << max_miscleavages << std::endl;
  }
  fisherCaller_.initConstraints(enzyme, digest, min_peptide_length, 
                                max_peptide_length, max_miscleavages);
  return true;
}

void FisherInterface::run() {
  std::map<std::string, std::string> fragment_map, duplicate_map;
  if (fastaProteinFN_ != "auto") {
    fisherCaller_.setFastaDatabase(fastaProteinFN_);
    
    if (VERB > 1) {
      std::cerr << "Detecting protein fragments/duplicates in target database" << std::endl;
    }
    bool generateDecoys = false;
    fisherCaller_.getProteinFragmentsAndDuplicates(fragment_map, duplicate_map, generateDecoys);
    
    if (VERB > 1) {
      std::cerr << "Detecting protein fragments/duplicates in decoy database" << std::endl;
    }
    generateDecoys = true;
    fisherCaller_.getProteinFragmentsAndDuplicates(fragment_map, duplicate_map, generateDecoys);
  }
  
  std::map<std::string, std::set<std::string> > groupProteinIds;
  unsigned int numGroups = 0;
  for (vector<ScoreHolder>::iterator peptideIt = peptideScores_->begin(); 
          peptideIt != peptideScores_->end(); ++peptideIt) {
    std::string lastProteinId;
    std::set<std::string> proteinsInGroup;
    bool isFirst = true, isShared = false;
    
    if (peptideIt->p > maxPeptidePval_) continue;
    
    for (std::vector<std::string>::iterator protIt = peptideIt->pPSM->proteinIds.begin(); 
            protIt != peptideIt->pPSM->proteinIds.end(); protIt++) {
      std::string proteinId = *protIt;
      
      if (fragment_map.find(proteinId) != fragment_map.end()) {
        if (reportFragmentProteins_) proteinsInGroup.insert(*protIt);
        proteinId = fragment_map[proteinId];
      } else if (duplicate_map.find(proteinId) != duplicate_map.end()) {
        if (reportDuplicateProteins_) proteinsInGroup.insert(*protIt);
        proteinId = duplicate_map[proteinId];
      } else {
        proteinsInGroup.insert(*protIt);
      }
      
      if (isFirst) {
        lastProteinId = proteinId;
        isFirst = false;
      } else if (lastProteinId != proteinId) {
        isShared = true;
        break;
      }
    }
    
    if (proteinsInGroup.size() == 1) {
      lastProteinId = *(proteinsInGroup.begin());
    }
    
    if (!isShared) {
      Protein::Peptide *peptide = new Protein::Peptide(
          peptideIt->pPSM->getPeptideSequence(), peptideIt->isDecoy(),
			    peptideIt->p, peptideIt->pep, peptideIt->q, peptideIt->score);
      if (proteins_.find(lastProteinId) == proteins_.end()) {
        if (proteinsInGroup.size() > 1) {
          groupProteinIds[lastProteinId] = proteinsInGroup;
        }
        Protein *newprotein = new Protein(lastProteinId, peptideIt->isDecoy(),
            peptide, ++numGroups);
        proteins_.insert(std::make_pair(lastProteinId,newprotein));
        if (lastProteinId.find(decoyPattern_) == std::string::npos) {
          ++numberTargetProteins_;
        } else {
          ++numberDecoyProteins_;
        }
      } else {
        proteins_[lastProteinId]->setPeptide(peptide);
        if (proteinsInGroup.size() > 1) {
          groupProteinIds[lastProteinId].insert(proteinsInGroup.begin(), 
                                                 proteinsInGroup.end());
        }
      }
    }
  }
  
  if (reportFragmentProteins_ || reportDuplicateProteins_) {
    std::map<std::string, std::set<std::string> >::iterator groupIt;
    std::map<const std::string,Protein*>::iterator representIt;
    for (groupIt = groupProteinIds.begin(); groupIt != groupProteinIds.end(); ++groupIt) {
      representIt = proteins_.find(groupIt->first);
      if (representIt != proteins_.end()) {
        std::string newName = "";
        for (std::set<std::string>::iterator proteinIt = groupIt->second.begin(); proteinIt != groupIt->second.end(); ++proteinIt) {
          newName += *proteinIt + ",";
        }
        newName = newName.substr(0, newName.size() - 1);
        representIt->second->setName(newName);
      }
    }
  }
}

void FisherInterface::computeProbabilities(const std::string& fname) {
  for (std::map<const std::string,Protein*>::iterator it = proteins_.begin(); 
        it != proteins_.end(); it++) {
    std::vector<Protein::Peptide*> peptides = it->second->getPeptides();
    switch (protInferenceMethod_) {
      case FISHER: {
        double fisher = 0.0;
        int significantPeptides = 0;
        for (std::vector<Protein::Peptide*>::const_iterator itP = peptides.begin();
              itP != peptides.end(); itP++) {
          fisher += log((*itP)->p / maxPeptidePval_);
        }
        //double proteinPvalue = boost::math::gamma_q(peptides.size(), -1.0*fisher);
        double proteinPvalue = 0.0;
        if (proteinPvalue == 0.0) proteinPvalue = DBL_MIN;
        it->second->setP(proteinPvalue);
        it->second->setScore(proteinPvalue);
        break;
      } case PEPPROD: { // MaxQuant's strategy
        double logPepProd = 0.0;
        for (std::vector<Protein::Peptide*>::const_iterator itP = peptides.begin();
              itP != peptides.end(); itP++) {
          logPepProd += log((*itP)->pep);
        }
        it->second->setScore(logPepProd);
        break;
      } case BESTPEPT: {
        double maxScore = -1000.0;
        for (std::vector<Protein::Peptide*>::const_iterator itP = peptides.begin();
              itP != peptides.end(); itP++) {
          maxScore = std::max(maxScore, (*itP)->score);
        }
        it->second->setScore(-1.0*maxScore); // lower scores are better
        break;
      }
    }
  }
  
  std::vector<std::pair<std::string,Protein*> > protIdProtPairs(proteins_.begin(), proteins_.end());
  std::sort(protIdProtPairs.begin(), protIdProtPairs.end(), IntCmpScore());
  
  if (!usePi0_) pickedProteinStrategy(protIdProtPairs);
  
  std::vector<double> peps;
  estimatePEPs(protIdProtPairs, peps);
  
  for (size_t i = 0; i < protIdProtPairs.size(); ++i) {
    pepProteinMap_.insert(std::make_pair(peps[i], std::vector<std::string>(1, protIdProtPairs[i].first) ));
  }
}

/* less efficient than pickedProteinStrategy, but more forgiving in specifying the correct decoyPattern */
void FisherInterface::pickedProteinStrategySubstring(
    std::vector<std::pair<std::string,Protein*> >& protIdProtPairs) {
  std::vector<std::pair<std::string,Protein*> > pickedProtIdProtPairs;
  std::vector<std::string> targetProts, decoyProts;
  std::vector<std::pair<std::string,Protein*> >::iterator it = protIdProtPairs.begin();
  size_t numErased = 0;
  for (; it != protIdProtPairs.end(); ++it) {
    bool isDecoy = it->second->getIsDecoy();
    bool globalErase = false;
    std::string proteinName = it->second->getName(); 
    
    std::istringstream ss(proteinName); 
    std::string proteinId;
    while (std::getline(ss, proteinId, ',')) { // split name by comma
      bool localErase = false;
      if (isDecoy) {
        for (size_t j = 0; j < targetProts.size(); ++j) {
          if (proteinId.find(targetProts[j]) != std::string::npos) {
            localErase = true;
            break;
          }
        }
      } else {
        for (size_t j = 0; j < decoyProts.size(); ++j) {
          if (decoyProts[j].find(proteinId) != std::string::npos) {
            localErase = true;
            break;
          }
        }
      }
      if (!localErase) {
        if (isDecoy) decoyProts.push_back(proteinId);
        else targetProts.push_back(proteinId);
      } else {
        globalErase = true;
      }
    }
    if (globalErase) {
      if (isDecoy) --numberDecoyProteins_;
      else --numberTargetProteins_;
      proteins_.erase(it->first); // this is where the printed proteins are stored
      numErased += 1;
    } else {
      pickedProtIdProtPairs.push_back(*it);
    }
  }
  pickedProtIdProtPairs.swap(protIdProtPairs);
  
  if (numErased == 0) {
    std::cerr << "Warning: No target-decoy protein pairs found. Check if the "
      << "correct decoy pattern was specified with the -P flag, i.e. if the "
      << "target protein is called \"protA\" and the decoy protein "
      << "\"decoy_protA\", use the option \"-p decoy_\"." << std::endl;
  }
}


void FisherInterface::pickedProteinStrategy(
    std::vector<std::pair<std::string,Protein*> >& protIdProtPairs) {  
  std::vector<std::pair<std::string,Protein*> > pickedProtIdProtPairs;
  std::set<std::string> targetProts, decoyProts;
  std::vector<std::pair<std::string,Protein*> >::iterator it = protIdProtPairs.begin();
  size_t numErased = 0;
  for (; it != protIdProtPairs.end(); ++it) {
    bool isDecoy = it->second->getIsDecoy();
    bool globalErase = false;
    std::string proteinName = it->second->getName(); 
    
    std::istringstream ss(proteinName); 
    std::string proteinId;
    while (std::getline(ss, proteinId, ',')) { // split name by comma
      bool localErase = false;
      if (isDecoy) {
        std::string targetId = proteinId.substr(decoyPattern_.size());
        if (targetProts.find(targetId) != targetProts.end()) {
          localErase = true;
        } else {
          decoyProts.insert(proteinId);
        }
      } else {
        if (decoyProts.find(decoyPattern_ + proteinId) != decoyProts.end()) {
          localErase = true;
        } else {
          targetProts.insert(proteinId);
        }
      }
      if (localErase) globalErase = true;
    }
    if (globalErase) {
      if (isDecoy) --numberDecoyProteins_;
      else --numberTargetProteins_;
      proteins_.erase(it->first); // this is where the printed proteins are stored
      numErased += 1;
    } else {
      pickedProtIdProtPairs.push_back(*it);
    }
  }
  pickedProtIdProtPairs.swap(protIdProtPairs);
  
  if (numErased == 0) {
    std::cerr << "Warning: No target-decoy protein pairs found for the picked "
      << "protein strategy.\n  Check if the correct decoy pattern was specified "
      << "with the -P flag;\n  if the target protein is called \"protA\" and the "
      << "decoy protein \"decoy_protA\", use the option \"-P decoy_\"." 
      << std::endl << std::endl;
  }
}

void FisherInterface::estimatePEPs(
    std::vector<std::pair<std::string,Protein*> >& protIdProtPairs,
    std::vector<double>& peps) {
  std::vector<std::pair<double, bool> > combined;
  std::vector<double> pvals;
  switch (protInferenceMethod_) {
    case FISHER: { // if we have well calibrated p-values
      for (size_t i = 0; i < protIdProtPairs.size(); ++i) {
        double pValue = protIdProtPairs[i].second->getP();
        bool isDecoy = protIdProtPairs[i].second->getIsDecoy();
        combined.push_back(make_pair(pValue, !isDecoy));
        if (!isDecoy) {
          pvals.push_back(pValue);
        }
      }
      std::sort(combined.begin(), combined.end());
      std::sort(pvals.begin(), pvals.end());
      break;
    } case PEPPROD:
      case BESTPEPT: { // if we have some other type of score
      for (size_t i = 0; i < protIdProtPairs.size(); ++i) {
        double score = protIdProtPairs[i].second->getScore();
        bool isDecoy = protIdProtPairs[i].second->getIsDecoy();
        combined.push_back(make_pair(score, !isDecoy));
      }
      std::sort(combined.begin(), combined.end());
      PosteriorEstimator::getPValues(combined, pvals);
      break;
    }
  }
  
  if (usePi0_) {
    pi0_ = PosteriorEstimator::estimatePi0(pvals);
    if (VERB > 1) {
      std::cerr << "protein pi0 estimate = " << pi0_ << std::endl;
    }
  }
  
  bool includeNegativesInResult = true;
  PosteriorEstimator::setReversed(true);
  PosteriorEstimator::estimatePEP(combined, usePi0_, pi0_ * absenceRatio_, peps, includeNegativesInResult);
}

std::ostream& FisherInterface::printParametersXML(std::ostream &os) {
  return os;
}

string FisherInterface::printCopyright() {
  ostringstream oss;
  return oss.str();
}
