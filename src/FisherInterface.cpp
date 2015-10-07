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

void FisherInterface::run() {
  // previously set in setTargetandDecoyNames, but we can trash that
  numberTargetProteins_ = 0;
  numberDecoyProteins_ = 0;
  proteins.clear();
  
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
    // for each protein
    //if (!peptideIt->isDecoy()) {
    if (true) {
      std::string lastProteinId;
      std::set<std::string> proteinsInGroup;
      bool isFirst = true, isShared = false;
      for (set<string>::iterator protIt = peptideIt->pPSM->proteinIds.begin(); 
              protIt != peptideIt->pPSM->proteinIds.end(); protIt++) {
        std::string proteinId = *protIt;
        
        // just hacking a bit for our simulation scripts
        if (proteinId.find(decoyPattern_) == std::string::npos) {
          size_t found = proteinId.find_first_of("_");
          proteinId = proteinId.substr(found + 1);
        }
        
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
        newName = newName.substr(0, newName.size() - 2);
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
  bool decoysPresent = false;
  for (size_t i = 0; i < myvec.size(); ++i) {
    double pValue = myvec[i].second->getP();
    bool isDecoy = myvec[i].second->getIsDecoy();
    combined.push_back(make_pair(pValue, !isDecoy));
    if (!isDecoy) {
      pvals.push_back(pValue);
    } else {
      decoysPresent = true;
    }
  }
  pi0_ = PosteriorEstimator::estimatePi0(pvals);
  
  if (VERB > 1) {
    std::cerr << "protein pi0 estimate = " << pi0_ << std::endl;
  }
  if (!decoysPresent) {
    size_t nDec = myvec.size();
    double step = 1.0 / 2.0 / (double)nDec;
    for (size_t ix = 0; ix < nDec; ++ix) {
      combined.push_back(std::make_pair(step * (1 + 2 * ix), false));
    }
  }
  std::sort(combined.begin(), combined.end());
  
  bool includeNegativesInResult = decoysPresent;
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
