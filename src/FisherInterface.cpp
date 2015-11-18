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

double get_pvalue(std::pair<std::string,Protein*> d) {
  return d.second->getP();
}

FisherInterface::FisherInterface(const std::string& fastaDatabase, 
    bool reportFragmentProteins, bool reportDuplicateProteins,
    bool trivialGrouping, double pi0, bool outputEmpirQval, 
    std::string& decoyPattern) :
  ProteinProbEstimator(trivialGrouping, pi0, outputEmpirQval, decoyPattern),
  fastaProteinFN_(fastaDatabase), reportFragmentProteins_(reportFragmentProteins),
  reportDuplicateProteins_(reportDuplicateProteins) {}

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
  
  for (vector<ScoreHolder>::iterator shIt = peptideScores_->begin(); shIt != peptideScores_->end(); ++shIt) {
    std::string peptideSequenceFlanked = shIt->pPSM->getFullPeptideSequence();
    
    peptideSequenceFlanked = PSMDescription::removePTMs(peptideSequenceFlanked);
    std::string peptideSequence = PSMDescription::removeFlanks(peptideSequenceFlanked);
    
    int peptide_length = peptideSequence.size();
    min_peptide_length = std::min(min_peptide_length, peptide_length);
    max_peptide_length = std::max(max_peptide_length, peptide_length);
    
    int miscleavages = Enzyme::countEnzymatic(peptideSequence);
    /*if (miscleavages > max_miscleavages) {
      std::cerr << "MC" << peptideSequenceFlanked << std::endl;
    }*/
    max_miscleavages = std::max(max_miscleavages, miscleavages);
    
    int non_enzymatic_flanks = 0;
    if (!Enzyme::isEnzymatic(peptideSequenceFlanked[0], 
                             peptideSequenceFlanked[2])) {
      ++non_enzymatic_flanks;
    }
    if (!Enzyme::isEnzymatic(peptideSequenceFlanked[peptideSequenceFlanked.size() - 3],
                             peptideSequenceFlanked[peptideSequenceFlanked.size() - 1])) {
      ++non_enzymatic_flanks;
    }
    if (non_enzymatic_flanks > max_non_enzymatic_flanks) {
      std::cerr << "SP " << peptideSequenceFlanked << std::endl;
    }
    max_non_enzymatic_flanks = std::max(max_non_enzymatic_flanks, non_enzymatic_flanks);
    total_non_enzymatic_flanks += max_non_enzymatic_flanks;
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
  fisherCaller_.setFastaDatabase(fastaProteinFN_);
  bool generateDecoys = false;
  fisherCaller_.getProteinFragmentsAndDuplicates(fragment_map, duplicate_map, generateDecoys);
  generateDecoys = true;
  fisherCaller_.getProteinFragmentsAndDuplicates(fragment_map, duplicate_map, generateDecoys);
  
  std::map<std::string, std::set<std::string> > groupProteinIds;
  unsigned int numGroups = 0;
  for (vector<ScoreHolder>::iterator peptideIt = peptideScores_->begin(); 
          peptideIt != peptideScores_->end(); ++peptideIt) {
    std::string lastProteinId;
    std::set<std::string> proteinsInGroup;
    bool isFirst = true, isShared = false;
    for (set<string>::iterator protIt = peptideIt->pPSM->proteinIds.begin(); 
            protIt != peptideIt->pPSM->proteinIds.end(); protIt++) {
      std::string proteinId = *protIt;
      
      /*
      // just hacking a bit for our simulation scripts
      if (proteinId.find(decoyPattern_) == std::string::npos) {
        size_t found = proteinId.find_first_of("_");
        proteinId = proteinId.substr(found + 1);
      }
      */
      
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
          peptideIt->pPSM->getPeptideSequence(),peptideIt->isDecoy(),
			    peptideIt->p, peptideIt->pep,peptideIt->q,peptideIt->p);
      if (proteins.find(lastProteinId) == proteins.end()) {
        if (proteinsInGroup.size() > 1) {
          groupProteinIds[lastProteinId] = proteinsInGroup;
        }
        Protein *newprotein = new Protein(lastProteinId,0.0,0.0,0.0,0.0,
            peptideIt->isDecoy(),peptide,++numGroups);
        proteins.insert(std::make_pair(lastProteinId,newprotein));
        if (lastProteinId.find(decoyPattern_) == std::string::npos) {
          ++numberTargetProteins_;
        } else {
          ++numberDecoyProteins_;
        }
      } else {
        proteins[lastProteinId]->setPeptide(peptide);
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
      representIt = proteins.find(groupIt->first);
      if (representIt != proteins.end()) {
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
  for(std::map<const std::string,Protein*>::iterator it = proteins.begin(); 
        it != proteins.end(); it++) {
    std::string proteinId = it->first;
    std::vector<Protein::Peptide*> peptides = it->second->getPeptides();
    double fisher = 0.0;
    for(std::vector<Protein::Peptide*>::const_iterator itP = peptides.begin();
          itP != peptides.end(); itP++) {
      fisher += log((*itP)->p);
    }
    double proteinPvalue = boost::math::gamma_q(peptides.size(), -1.0*fisher);
    if (proteinPvalue == 0.0) proteinPvalue = DBL_MIN;
    it->second->setP(proteinPvalue);
  }
  
  std::vector<std::pair<std::string,Protein*> > myvec(proteins.begin(), proteins.end());
  std::sort(myvec.begin(), myvec.end(), IntCmpPvalue());
  
  std::vector<std::pair<double, bool> > combined;
  std::vector<double> pvals;
  for (size_t i = 0; i < myvec.size(); ++i) {
    double pValue = myvec[i].second->getP();
    bool isDecoy = myvec[i].second->getIsDecoy();
    combined.push_back(make_pair(pValue, !isDecoy));
    if (!isDecoy) {
      pvals.push_back(pValue);
    }
  }
  pi0_ = PosteriorEstimator::estimatePi0(pvals);
  
  if (VERB > 1) {
    std::cerr << "protein pi0 estimate = " << pi0_ << std::endl;
  }
  std::sort(combined.begin(), combined.end());
  
  bool includeNegativesInResult = true;
  std::vector<double> peps;
  PosteriorEstimator::setReversed(true);
  PosteriorEstimator::estimatePEP(combined, usePi0_, pi0_, peps, includeNegativesInResult);
  
  for (size_t i = 0; i < myvec.size(); ++i) {
    pvalues.push_back(myvec[i].second->getP());
    pepProteinMap_.insert(std::make_pair(peps[i], std::vector<std::string>(1, myvec[i].first) ));
  }
}

std::ostream& FisherInterface::printParametersXML(std::ostream &os) {
  return os;
}

string FisherInterface::printCopyright() {
  ostringstream oss;
  return oss.str();
}
