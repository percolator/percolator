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
#include <cmath>
#include <assert.h>

#include "Globals.h"
#include "PSMDescriptionDOC.h"

double PSMDescriptionDOC::normDivRT_ = -1.0;
double PSMDescriptionDOC::normSubRT_ = 0.0;

PSMDescriptionDOC::PSMDescriptionDOC() : 
    PSMDescription(), retentionFeatures_(NULL), parentFragment_(NULL),
    pI_(0.0), massDiff_(0.0), retentionTime_(0.0), predictedTime_(0.0) {}

PSMDescriptionDOC::PSMDescriptionDOC(const string pep, const double retTime) :
    PSMDescription(pep), retentionFeatures_(NULL), parentFragment_(NULL),
    pI_(0.0), massDiff_(0.0), retentionTime_(retTime), predictedTime_(0.0) {}

PSMDescriptionDOC::~PSMDescriptionDOC() {}

void PSMDescriptionDOC::deleteRetentionFeatures() {
  if (retentionFeatures_) {
    delete[] retentionFeatures_;
    retentionFeatures_ = NULL;
  }
}

void PSMDescriptionDOC::setPSMSet(vector<PSMDescription*> & psms) {
  double minRT = 1e10, maxRT = -1;
  vector<PSMDescription*>::iterator psm;
  for (psm = psms.begin(); psm != psms.end(); ++psm) {
    minRT = min(minRT, (*psm)->getRetentionTime());
    maxRT = max(maxRT, (*psm)->getRetentionTime());
  }
  normDivRT_ = (maxRT - minRT) / 2.;
  normSubRT_ = minRT + normDivRT_;
  if (normDivRT_ == 0.0) {
    normDivRT_ = 1.0;
  }
}

void PSMDescriptionDOC::setRetentionTime(vector<PSMDescription*>& psms, map<int, double>& scan2rt) {
  vector<PSMDescription*>::iterator psm = psms.begin();
  if (scan2rt.size() == 0) {
    if ((*psm)->getRetentionTime() > 0) {
      if (normDivRT_ < 0) {
        double minRT = 1e10, maxRT = -1;
        for (; psm != psms.end(); ++psm) {
          minRT = min(minRT, (*psm)->getRetentionTime());
          maxRT = max(maxRT, (*psm)->getRetentionTime());
        }
        psm = psms.begin();
        normDivRT_ = (maxRT - minRT) / 2.;
        normSubRT_ = minRT + normDivRT_;
        if (normDivRT_ == 0.0) {
          normDivRT_ = 1.0;
        }
      }
      for (; psm != psms.end(); ++psm) {
        (*psm)->setUnnormalizedRetentionTime((*psm)->getRetentionTime());
        if (psm == psms.begin()) {
          continue;
        }
        vector<PSMDescription*>::reverse_iterator rpsm(psm);
        (*psm)->checkFragmentPeptides(rpsm, psms.rend());
      }
    } else {
      if (VERB > 1) {
        cerr << "Approximating retention time with scan number." << endl;
      }
      if (normDivRT_ < 0) {
        double minRT = (double)(*psm)->scan, diffRT = (*psms.rbegin())->scan - (*psm)->scan;
        normDivRT_ = diffRT / 2.;
        normSubRT_ = minRT + normDivRT_;
        if (normDivRT_ == 0.0) {
          normDivRT_ = 1.0;
        }
      }
      for (; psm != psms.end(); ++psm) {
        (*psm)->setUnnormalizedRetentionTime((*psm)->scan);
        if (psm == psms.begin()) {
          continue;
        }
        vector<PSMDescription*>::reverse_iterator rpsm(psm);
        (*psm)->checkFragmentPeptides(rpsm, psms.rend());
      }
    }
  } else {
    if (normDivRT_ < 0) {
      double minRT = scan2rt.begin()->second, diffRT = scan2rt.rbegin()->second - minRT;
      normDivRT_ = diffRT / 2.;
      normSubRT_ = minRT + normDivRT_;
      if (normDivRT_ == 0.0) {
        normDivRT_ = 1.0;
      }
    }
    for (; psm != psms.end(); ++psm) {
      assert(scan2rt.count((*psm)->scan) > 0);
      (*psm)->setUnnormalizedRetentionTime(scan2rt[(*psm)->scan]);
      if (psm == psms.begin()) {
        continue;
      }
      vector<PSMDescription*>::reverse_iterator rpsm(psm);
      (*psm)->checkFragmentPeptides(rpsm, psms.rend());
    }
  }
}

void PSMDescriptionDOC::checkFragmentPeptides(
    std::vector<PSMDescription*>::reverse_iterator other,
    std::vector<PSMDescription*>::reverse_iterator theEnd) {
  for (; other != theEnd; ++other) {
    if (abs(getRetentionTime() - (*other)->getRetentionTime()) > 0.02) {
      return;
    }
    if (isSubPeptide(peptide, (*other)->getFullPeptide())) {
      if (parentFragment_ == NULL
          || parentFragment_->getFullPeptide().length()
              < (*other)->getFullPeptide().length()) {
        parentFragment_ = (*other)->getAParent();
        //        cerr << parentFragment_->getFullPeptide() << " " << peptide << endl;
      }
    }
    if (isSubPeptide((*other)->peptide, getFullPeptide())) {
      if ((*other)->getParentFragment() == NULL
          || (*other)->getParentFragment()->peptide.length()
              < getFullPeptide().length()) {
        (*other)->setParentFragment(getAParent());
        //          cerr << getFullPeptide() << " " << getFullPeptide().length() << " " << other->peptide << " " << other->peptide.length() << endl;
      }
    }
  }
}

bool PSMDescriptionDOC::isSubPeptide(string& child, string& parent) {
  size_t len = parent.length();
  if (!(Enzyme::isEnzymatic(parent[0], parent[2])
      && Enzyme::isEnzymatic(parent[len - 3], parent[len - 1]))) {
    return false;
  }
  string strippedChild = child.substr(2, child.length() - 4);
  size_t found = parent.find(strippedChild);
  if (found == string::npos) {
    return false;
  }
  return parent.length() > child.length();
}

void PSMDescriptionDOC::normalizeRetentionTimes(vector<PSMDescription*> & psms) {
  vector<PSMDescription*>::iterator psm;
  for (psm = psms.begin(); psm != psms.end(); ++psm) {
    (*psm)->setUnnormalizedRetentionTime((*psm)->getRetentionTime());
  }
}

void PSMDescriptionDOC::unnormalizeRetentionTimes(vector<PSMDescription*> & psms) {
  vector<PSMDescription*>::iterator psm;
  for (psm = psms.begin(); psm != psms.end(); ++psm) {
    (*psm)->setRetentionTime(unnormalize((*psm)->getRetentionTime()));
  }
}

double PSMDescriptionDOC::normalize(double unnormalizedTime) {
  return (unnormalizedTime - normSubRT_) / normDivRT_;
}

double PSMDescriptionDOC::unnormalize(double normalizedTime) {
  return normalizedTime * normDivRT_ + normSubRT_;
}

std::vector<double*> PSMDescriptionDOC::getRetFeatures(
    std::vector<PSMDescription*>& psms) {
  vector<double*> features;
  vector<PSMDescription*>::iterator psm;
  for (psm = psms.begin(); psm != psms.end(); ++psm) {
    features.push_back((*psm)->getRetentionFeatures());
  }
  return features;
}

