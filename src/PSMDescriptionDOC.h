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
#ifndef PSM_DESCRIPTION_DOC_H_
#define PSM_DESCRIPTION_DOC_H_
#include <map>
#include <vector>
#include <string>
#include <iostream>

#include "PSMDescription.h"

/*
* PSMDescriptionDOC
*
* Here are some useful abbreviations:
* PSM - Peptide Spectrum Match
* DOC - Description Of Correct
*
*/
class PSMDescriptionDOC : public PSMDescription {
  
 public: 
  PSMDescriptionDOC();
  PSMDescriptionDOC(const std::string peptide, const double retTime);
  
  ~PSMDescriptionDOC();
  
  inline void setRetentionFeatures(double* retentionFeatures) { 
    retentionFeatures_ = retentionFeatures; 
  }
  inline double* getRetentionFeatures() { return retentionFeatures_; }
  void deleteRetentionFeatures();
  
  inline void setIsoElectricPoint(const double pI) { pI_ = pI; }
  inline double getIsoElectricPoint() const { return pI_; }
  
  inline void setRetentionTime(const double retentionTime) {
    retentionTime_ = retentionTime;
  }
  inline double getRetentionTime() const { return retentionTime_; }
  
  inline void setUnnormalizedRetentionTime(const double retentionTime) {
    retentionTime_ = normalize(retentionTime);
  }
  inline double getUnnormalizedRetentionTime() const { return unnormalize(retentionTime_); }
  
  inline void setPredictedRetentionTime(const double predictedTime) {
    predictedTime_ = predictedTime;
  }
  inline double getPredictedRetentionTime() const { return predictedTime_; }
  
  inline void setMassDiff(const double dm) {
    massDiff_ = dm;
  }
  inline double getMassDiff() const { return massDiff_; }
  
  std::string& getFullPeptide() { return getAParent()->peptide; }
  PSMDescription* getAParent() {
    if (parentFragment_) return parentFragment_->getAParent();
    else return this;
  }
  
  friend std::ostream& operator<<(std::ostream& out, PSMDescriptionDOC& psm);
  
  // static methods and members for retention time normalization
  static double normDivRT_, normSubRT_;
  
  static void setPSMSet(std::vector<PSMDescription*>& psms);
  static void normalizeRetentionTimes(std::vector<PSMDescription*>& psms);
  static double normalize(double unnormalizedTime);
  static void unnormalizeRetentionTimes(std::vector<PSMDescription*>& psms);
  static double unnormalize(double normalizedTime);
  
  static std::vector<double*> getRetFeatures(std::vector<PSMDescription*>& psms);
 private:
  double* retentionFeatures_;
  PSMDescription* parentFragment_;
  double pI_, massDiff_, predictedTime_, retentionTime_;
};

inline std::ostream& operator<<(std::ostream& out, PSMDescriptionDOC& psm) {
  out << "Peptide: " << psm.peptide << endl;
  out << "Spectrum scan number: " << psm.scan << endl;
  out << "Retention time, predicted retention time: " << psm.retentionTime_
      << ", " << psm.predictedTime_;
  out << "Retention features: ";
  for (int i = 0; i < 62; ++i) {
    out << psm.retentionFeatures_[i] << "  ";
  }
  out << endl;
  return out;
}

#endif /*PSM_DESCRIPTION_DOC_H_*/
