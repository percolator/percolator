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
#ifndef PSMDESCRIPTION_H_
#define PSMDESCRIPTION_H_
#include <set>
#include <map>
#include <vector>
#include <string>
#include <iostream>

#include "Enzyme.h"

/*
* PSMDescription
*
* Here are some useful abbreviations:
* PSM - Peptide Spectrum Match
*
*/
class PSMDescription {
  
 public: 
  PSMDescription();
  PSMDescription(const std::string peptide, const double retTime);
  PSMDescription(double ort, double prt);
  
  virtual ~PSMDescription();
  
  static void deletePtr(PSMDescription* psm);
  void clear() { proteinIds.clear(); }
  double* getFeatures() { return features; }
  double* getRetentionFeatures() { return retentionFeatures; }
  static std::vector<double*> getRetFeatures(std::vector<PSMDescription> & psms);
  
  // TODO: move these 2 functions somewhere else
  static std::string removePTMs(const std::string& peptideSeq);
  static std::string removeFlanks(const std::string& peptideSeq) { 
    return peptideSeq.substr(2, peptideSeq.size()-4); 
  }
  
  std::string& getFullPeptide() { return getAParent()->peptide; }
  std::string getPeptideSequence() { return peptide.substr(2, peptide.size()-4); }
  std::string& getFullPeptideSequence() { return peptide; }
  std::string getFlankN() { return peptide.substr(0, 1); }    
  std::string getFlankC() { return peptide.substr(peptide.size()-1, peptide.size()); }
  
  PSMDescription* getAParent() {
    if (parentFragment) {
      return parentFragment->getAParent();
    } else {
      return this;
    }
  }
  
  double getUnnormalizedRetentionTime() { return unnormalize(retentionTime); }
  static bool isSubPeptide(std::string& child, std::string& parent);
  bool isNotEnzymatic() {
    std::string peptideSeq = removePTMs(peptide);
    std::string peptideSeqNoFlanks = removeFlanks(peptide);
    return !(Enzyme::isEnzymatic(peptideSeq[0], peptideSeq[2])
        && Enzyme::isEnzymatic(peptideSeq[peptideSeq.size() - 3],
                               peptideSeq[peptideSeq.size() - 1])
        && Enzyme::countEnzymatic(peptideSeqNoFlanks) == 0);
  }
  void checkFragmentPeptides(std::vector<PSMDescription*>::reverse_iterator other,
                             std::vector<PSMDescription*>::reverse_iterator theEnd);
  
  static void setRetentionTime(std::vector<PSMDescription*>& psms, 
                               std::map<int, double>& scan2rt);
  static double unnormalize(double normalizedTime);
  static void unnormalizeRetentionTimes(std::vector<PSMDescription> & psms);
  // set the norm and div for a set of peptides
  static void setPSMSet(std::vector<PSMDescription> & psms);
  // normalize retention times for a  set of peptides
  static void normalizeRetentionTimes(std::vector<PSMDescription> & psms);
  friend std::ostream& operator<<(std::ostream& out, PSMDescription& psm);
  void printProteins(std::ostream& out);
  double getRetentionTime() const { return retentionTime; }
  double getPredictedRetentionTime() { return predictedTime; }

  static double normDiv, normSub;
  
  double* features; // owned by a FeatureMemoryPool instance, no need to delete
  double* retentionFeatures;
  double retentionTime, predictedTime, massDiff, pI, expMass, calcMass;
  unsigned int scan;
  std::string id;
  std::string peptide;
  std::set<std::string> proteinIds;
  PSMDescription* parentFragment;
};

inline bool const operator<(PSMDescription const& one,
                            PSMDescription const& other) {
  if (one.peptide == other.peptide) {
    return one.retentionTime < other.retentionTime;
  }
  return one.peptide < other.peptide;
}

inline bool operator==(PSMDescription const& one,
                       PSMDescription const& other) {
  if (one.peptide == other.peptide) {
    return true;
  } else {
    return false;
  }
}

inline std::ostream& operator<<(std::ostream& out, PSMDescription& psm) {
  out << "Peptide: " << psm.peptide << endl;
  out << "Spectrum scan number: " << psm.scan << endl;
  out << "Retention time, predicted retention time: " << psm.retentionTime
      << ", " << psm.predictedTime;
  out << "Retention features: ";
  for (int i = 0; i < 62; ++i) {
    out << psm.retentionFeatures[i] << "  ";
  }
  out << endl;
  return out;
}

#endif /*PSMDESCRIPTION_H_*/
