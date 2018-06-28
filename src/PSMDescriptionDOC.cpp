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

