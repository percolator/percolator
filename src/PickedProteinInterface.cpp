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

#include "PickedProteinInterface.h"

using namespace PercolatorCrux;

PickedProteinInterface::PickedProteinInterface(const std::string& fastaDatabase,
    double pvalueCutoff, bool reportFragmentProteins, bool reportDuplicateProteins,
    bool trivialGrouping, double absenceRatio, bool outputEmpirQval, 
    std::string& decoyPattern, double specCountQvalThreshold) :
      ProteinProbEstimator(trivialGrouping, absenceRatio, outputEmpirQval, 
                           decoyPattern, specCountQvalThreshold),
      fastaProteinFN_(fastaDatabase), maxPeptidePval_(pvalueCutoff),
      reportFragmentProteins_(reportFragmentProteins),
      reportDuplicateProteins_(reportDuplicateProteins),
      protInferenceMethod_(BESTPEPT) {
  if (absenceRatio == 1.0) usePi0_ = false;
}

PickedProteinInterface::~PickedProteinInterface() {}

bool PickedProteinInterface::initialize(Scores& peptideScores, const Enzyme* enzyme) {
  int min_peptide_length = 1000, max_peptide_length = 0;
  int max_miscleavages = 0, max_non_enzymatic_flanks = 0;
  int max_miscleavages_trypsinp = 0, max_non_enzymatic_flanks_trypsinp = 0;
  int total_non_enzymatic_flanks = 0;
  
  std::string enzymeString = enzyme->getStringEnzyme();
  
  ENZYME_T cruxEnzyme;
  bool tryTrypsinP = false;
  Enzyme* trypsinP = NULL;  
  switch (enzyme->getEnzymeType()) {
    case Enzyme::TRYPSIN:
      cruxEnzyme = TRYPSIN;
      tryTrypsinP = true;
      trypsinP = Enzyme::createEnzyme(Enzyme::TRYPSINP);
      break;
    case Enzyme::TRYPSINP:
      cruxEnzyme = TRYPSINP;
      break;
    case Enzyme::CHYMOTRYPSIN:
      cruxEnzyme = CHYMOTRYPSIN;
      break;
    case Enzyme::ELASTASE:
      cruxEnzyme = ELASTASE;
      break;
    case Enzyme::LYSN:
      cruxEnzyme = LYSN;
      break;
    case Enzyme::LYSC:
      cruxEnzyme = LYSC;
      break;
    case Enzyme::ARGC:
      cruxEnzyme = ARGC;
      break;
    case Enzyme::ASPN:
      cruxEnzyme = ASPN;
      break;
    case Enzyme::GLUC:
      cruxEnzyme = GLUC;
      break;
    case Enzyme::NO_ENZYME:
    default:
      if (enzyme->getEnzymeType() != Enzyme::NO_ENZYME) {
        if (VERB > 1) {
          std::cerr << "Warning: specified protease " << enzyme->getStringEnzyme()
                << " currently not supported. Using unspecific digestion to identify"
                << " duplicate and fragment proteins" << std::endl;
        }
        enzymeString = NoEnzyme::getString();
      }
      cruxEnzyme = NO_ENZYME;
      max_non_enzymatic_flanks = 2;
      break;
  }
  
  if (max_non_enzymatic_flanks != 2) {
    for (std::vector<ScoreHolder>::iterator shIt = peptideScores.begin(); 
             shIt != peptideScores.end(); ++shIt) {    
      std::string peptideSequenceFlanked = shIt->pPSM->getFullPeptideSequence();
      
      peptideSequenceFlanked = PSMDescription::removePTMs(peptideSequenceFlanked);
      std::string peptideSequence = PSMDescription::removeFlanks(peptideSequenceFlanked);
      
      int peptide_length = peptideSequence.size();
      min_peptide_length = std::min(min_peptide_length, peptide_length);
      max_peptide_length = std::max(max_peptide_length, peptide_length);
      
      int miscleavages = enzyme->countEnzymatic(peptideSequence);
      if (miscleavages > max_miscleavages) {
        if (VERB > 1) {
          std::cerr << "Miscleavage detected: " << peptideSequenceFlanked << std::endl;
        }
        max_miscleavages = std::max(max_miscleavages, miscleavages);
      }
      
      if (tryTrypsinP) {
        int miscleavages_trypsinp = trypsinP->countEnzymatic(peptideSequence);
        if (miscleavages_trypsinp > max_miscleavages_trypsinp) {
          max_miscleavages_trypsinp = std::max(max_miscleavages_trypsinp, miscleavages_trypsinp);
        }
      }
      
      int non_enzymatic_flanks = 0;
      int non_enzymatic_flanks_trypsinp = 0;
      if (!enzyme->isEnzymatic(peptideSequenceFlanked[0], 
                               peptideSequenceFlanked[2]) 
             && peptideSequenceFlanked.substr(0,1) != "M") { // ignore protein N-term methionine cleavage
        ++non_enzymatic_flanks;
        if (tryTrypsinP && !trypsinP->isEnzymatic(peptideSequenceFlanked[0], 
                               peptideSequenceFlanked[2]) 
             && peptideSequenceFlanked.substr(0,1) != "M") {
          ++non_enzymatic_flanks_trypsinp;
        }
      }
      if (!enzyme->isEnzymatic(peptideSequenceFlanked[peptideSequenceFlanked.size() - 3],
                               peptideSequenceFlanked[peptideSequenceFlanked.size() - 1])) {
        ++non_enzymatic_flanks;
        if (tryTrypsinP && !trypsinP->isEnzymatic(peptideSequenceFlanked[peptideSequenceFlanked.size() - 3], 
                               peptideSequenceFlanked[peptideSequenceFlanked.size() - 1])) {
          ++non_enzymatic_flanks_trypsinp;
        }
      }
      if (non_enzymatic_flanks > max_non_enzymatic_flanks) {
        if (VERB > 1) {
          std::cerr << "Non enzymatic flank detected: " << peptideSequenceFlanked << std::endl;
        }
        max_non_enzymatic_flanks = std::max(max_non_enzymatic_flanks, non_enzymatic_flanks);
      }
      total_non_enzymatic_flanks += non_enzymatic_flanks;
      
      if (tryTrypsinP && non_enzymatic_flanks_trypsinp > max_non_enzymatic_flanks_trypsinp) {
        max_non_enzymatic_flanks_trypsinp = std::max(max_non_enzymatic_flanks_trypsinp, non_enzymatic_flanks_trypsinp);
      }
      
    }
    
    if (tryTrypsinP && max_non_enzymatic_flanks_trypsinp < max_non_enzymatic_flanks) {
      if (VERB > 1) {
        std::cerr << "Detected TrypsinP as protease instead of Trypsin, allowing (R|K).P cleavages." << std::endl;
      }
      max_non_enzymatic_flanks = max_non_enzymatic_flanks_trypsinp;
      max_miscleavages = max_miscleavages_trypsinp;
      cruxEnzyme = TRYPSINP;
      enzymeString = TrypsinP::getString();
    }
    
    // if more than half of the flanks or non enzymatic, probably the protease is wrong
    if (total_non_enzymatic_flanks > peptideScores.size()) {
      std::cerr << "Warning: more than half of the cleavage sites are non enzymatic, "
                << "please verify that the right protease was specified." << std::endl;
    }
    
    if (trypsinP) {
      delete trypsinP;
    }
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
    cruxEnzyme = TRYPSIN; // Use trypsin to create candidate protein groups
    digestString = "non-specific";
  }
  if (VERB > 1) {
    std::cerr << "Protein digestion parameters for duplicate/fragment detection (detected from PSM input):\n"
              << " enzyme=" << enzymeString 
              << ", digestion=" << digestString
              << ", min-pept-length=" << min_peptide_length
              << ", max-pept-length=" << max_peptide_length
              << ", max-miscleavages=" << max_miscleavages << std::endl;
  }
  fisherCaller_.initConstraints(cruxEnzyme, digest, min_peptide_length, 
                                max_peptide_length, max_miscleavages);
  
  groupProteins(peptideScores);
  
  return true;
}

void PickedProteinInterface::groupProteins(Scores& peptideScores) {
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
  for (vector<ScoreHolder>::iterator peptideIt = peptideScores.begin(); 
          peptideIt != peptideScores.end(); ++peptideIt) {
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
      ProteinScoreHolder::Peptide peptide(peptideIt->pPSM->getPeptideSequence(), 
          peptideIt->isDecoy(), peptideIt->p, peptideIt->pep, peptideIt->q, peptideIt->score);
      if (proteinToIdxMap_.find(lastProteinId) == proteinToIdxMap_.end()) {
        if (proteinsInGroup.size() > 1) {
          groupProteinIds[lastProteinId] = proteinsInGroup;
        }
        ProteinScoreHolder newProtein(lastProteinId, peptideIt->isDecoy(),
            peptide, ++numGroups);
        proteinToIdxMap_[lastProteinId] = proteins_.size();
        proteins_.push_back(newProtein);
        if (lastProteinId.find(decoyPattern_) == std::string::npos) {
          ++numberTargetProteins_;
        } else {
          ++numberDecoyProteins_;
        }
      } else {
        proteins_.at(proteinToIdxMap_[lastProteinId]).addPeptide(peptide);
        if (proteinsInGroup.size() > 1) {
          groupProteinIds[lastProteinId].insert(proteinsInGroup.begin(), 
                                                 proteinsInGroup.end());
        }
      }
    }
  }
  
  /* Update protein group identifier to include fragment and duplicate protein identifiers */
  if (reportFragmentProteins_ || reportDuplicateProteins_) {
    std::map<std::string, std::set<std::string> >::iterator groupIt;
    std::map<std::string, size_t>::iterator representIt;
    for (groupIt = groupProteinIds.begin(); groupIt != groupProteinIds.end(); ++groupIt) {
      representIt = proteinToIdxMap_.find(groupIt->first);
      if (representIt != proteinToIdxMap_.end()) { /* these are the protein group representatives */
        std::string newName = "";
        for (std::set<std::string>::iterator proteinIt = groupIt->second.begin(); proteinIt != groupIt->second.end(); ++proteinIt) {
          /* some protein identifiers contain commas, replace them by the much 
             less used semicolon as the comma is used to separate protein 
             identifiers */
          std::string proteinId = *proteinIt;
          std::replace(proteinId.begin(), proteinId.end(), ',', ';'); 
          newName += proteinId + ",";
        }
        newName = newName.substr(0, newName.size() - 1); // remove last comma
        proteins_.at(representIt->second).setName(newName);
      }
    }
  }
}

void PickedProteinInterface::computeProbabilities(const std::string& fname) {
  if (VERB > 2) {
    std::cerr << "Computing protein probabilities for " 
              << proteins_.size() << " protein groups." << std::endl;
  }
  for (std::vector<ProteinScoreHolder>::iterator it = proteins_.begin(); 
        it != proteins_.end(); it++) {
    const std::vector<ProteinScoreHolder::Peptide> peptides = it->getPeptidesByRef();
    switch (protInferenceMethod_) {
      case FISHER: {
        double fisher = 0.0;
        int significantPeptides = 0;
        for (std::vector<ProteinScoreHolder::Peptide>::const_iterator itP = peptides.begin();
              itP != peptides.end(); itP++) {
          fisher += log(itP->p / maxPeptidePval_);
        }
        //double proteinPvalue = boost::math::gamma_q(peptides.size(), -1.0*fisher);
        double proteinPvalue = 0.0;
        if (proteinPvalue == 0.0) proteinPvalue = DBL_MIN;
        it->setP(proteinPvalue);
        it->setScore(proteinPvalue);
        break;
      } case PEPPROD: { // MaxQuant's strategy
        double logPepProd = 0.0;
        for (std::vector<ProteinScoreHolder::Peptide>::const_iterator itP = peptides.begin();
              itP != peptides.end(); itP++) {
          logPepProd += log(itP->pep);
        }
        it->setScore(logPepProd);
        break;
      } case BESTPEPT: {
        double maxScore = -1000.0;
        for (std::vector<ProteinScoreHolder::Peptide>::const_iterator itP = peptides.begin();
              itP != peptides.end(); itP++) {
          maxScore = std::max(maxScore, itP->score);
        }
        it->setScore(-1.0*maxScore); // lower scores are better
        break;
      }
    }
  }
  
  std::sort(proteins_.begin(), proteins_.end(), IntCmpScore());
  
  if (!usePi0_) {
    pickedProteinStrategy();
  }
  
  estimatePEPs();
}

bool PickedProteinInterface::pickedProteinCheckId(std::string& proteinId, bool isDecoy,
    std::set<std::string>& targetProts, std::set<std::string>& decoyProts) {
  bool found = false;
  if (isDecoy) {
    std::string targetId = proteinId.substr(decoyPattern_.size());
    if (targetProts.find(targetId) != targetProts.end()) {
      found = true;
    } else {
      decoyProts.insert(proteinId);
    }
  } else {
    if (decoyProts.find(decoyPattern_ + proteinId) != decoyProts.end()) {
      found = true;
    } else {
      targetProts.insert(proteinId);
    }
  }
  return found;
}

bool PickedProteinInterface::pickedProteinCheck(std::string& proteinName, bool isDecoy, 
    std::set<std::string>& targetProts, std::set<std::string>& decoyProts) {
  bool erase = false;
  if (reportFragmentProteins_ || reportDuplicateProteins_) {
    std::istringstream ss(proteinName); 
    std::string proteinId;
    while (std::getline(ss, proteinId, ',')) { // split name by comma
      erase = erase || pickedProteinCheckId(proteinId, isDecoy, 
                                          targetProts, decoyProts);
    }
  } else {
    erase = pickedProteinCheckId(proteinName, isDecoy, targetProts, decoyProts);
  }
  return erase;
}

/* Executes the picked protein-FDR strategy from Savitski et al. 2015
   For protein groups, if one of the corresponding proteins has been observed
   the whole group is eliminated */
void PickedProteinInterface::pickedProteinStrategy() {
  if (VERB > 1) {
    std::cerr << "Performing picked protein strategy" << std::endl;
  }
  
  std::vector<ProteinScoreHolder> pickedProtIdProtPairs;
  std::set<std::string> targetProts, decoyProts;
  std::vector<ProteinScoreHolder>::iterator it = proteins_.begin();
  size_t numErased = 0;
  // TODO: what about peptides with both target and decoy proteins?
  for (; it != proteins_.end(); ++it) {
    bool isDecoy = it->isDecoy();
    std::string proteinName = it->getName(); 
    
    bool erase = pickedProteinCheck(proteinName, isDecoy, targetProts, decoyProts);
    if (erase) {
      if (isDecoy) --numberDecoyProteins_;
      else --numberTargetProteins_;
      numErased += 1;
    } else {
      pickedProtIdProtPairs.push_back(*it);
    }
  }
  pickedProtIdProtPairs.swap(proteins_);
  
  if (numErased == 0) {
    std::cerr << "Warning: No target-decoy protein pairs found for the picked "
      << "protein strategy.\n  Check if the correct decoy pattern was specified "
      << "with the -P flag;\n  if the target protein is called \"protA\" and the "
      << "decoy protein \"decoy_protA\", use the option \"-P decoy_\"." 
      << std::endl << std::endl;
  }
  
  if (VERB > 1) {
    std::cerr << "Eliminated lower-scoring target-decoy protein: "
              << targetProts.size() << " target proteins and "
              << decoyProts.size() << " decoy proteins remaining." << std::endl;
  }
}

void PickedProteinInterface::estimatePEPs() {
  std::vector<std::pair<double, bool> > combined;
  std::vector<double> pvals;
  switch (protInferenceMethod_) {
    case FISHER: { // if we have well calibrated p-values
      for (size_t i = 0; i < proteins_.size(); ++i) {
        double pValue = proteins_.at(i).getP();
        bool isTarget = proteins_.at(i).isTarget();
        combined.push_back(make_pair(pValue, isTarget));
        if (isTarget) {
          pvals.push_back(pValue);
        }
      }
      std::sort(combined.begin(), combined.end());
      std::sort(pvals.begin(), pvals.end());
      break;
    } case PEPPROD:
      case BESTPEPT: { // if we have some other type of score
      for (size_t i = 0; i < proteins_.size(); ++i) {
        double score = proteins_.at(i).getScore();
        bool isTarget = proteins_.at(i).isTarget();
        combined.push_back(make_pair(score, isTarget));
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
  
  std::vector<double> peps;
  bool includeNegativesInResult = true;
  PosteriorEstimator::setReversed(true);
  PosteriorEstimator::estimatePEP(combined, usePi0_, pi0_, peps, includeNegativesInResult);
  
  for (size_t i = 0; i < proteins_.size(); ++i) {
    proteins_.at(i).setPEP(peps.at(i));
  }
}

std::ostream& PickedProteinInterface::printParametersXML(std::ostream &os) {
  return os;
}

string PickedProteinInterface::printCopyright() {
  ostringstream oss;
  return oss.str();
}
