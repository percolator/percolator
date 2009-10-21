/*******************************************************************************
    Copyright 2006-2009 Lukas Käll <lukas.kall@cbr.su.se>

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
#include <assert.h>
#include "Globals.h"
#include "PSMDescription.h"
#include "DescriptionOfCorrect.h"

PSMDescription::PSMDescription():  q(0.),pep(0.),features(NULL),retentionFeatures(NULL),
                                   retentionTime(0.),predictedTime(0.),massDiff(0.),pI(0.),
                                   scan(0),id(""),peptide("")
{
}

PSMDescription::PSMDescription(const string pep, const double retTime):  q(0.),pep(0.),features(NULL),retentionFeatures(NULL),
                                   retentionTime(retTime),predictedTime(0.),massDiff(0.),pI(0.),
                                   scan(0),id(""),peptide(pep)
{
}


PSMDescription::~PSMDescription()
{
}

double PSMDescription::normDiv=1.0;
double PSMDescription::normSub=0.0;

double PSMDescription::unnormalize(double normalizedTime) {
	return normalizedTime*normDiv+normSub;
}

void PSMDescription::setRetentionTime(vector<PSMDescription>& psms, map<int,double>& scan2rt) {
  vector<PSMDescription>::iterator psm = psms.begin();
  if (scan2rt.size() == 0) {
	if (psm->retentionTime>0) {
	  double minRT=1e10,maxRT=-1;
      for(; psm != psms.end(); ++psm) {
        minRT=min(minRT,psm->retentionTime);
        maxRT=max(maxRT,psm->retentionTime);
      }
      psm = psms.begin();
      normDiv=(maxRT-minRT)/2.;
      normSub=minRT+normDiv;
      if (normDiv==0.0) normDiv = 1.0;
      for(; psm != psms.end(); ++psm) {
        psm->retentionTime = (psm->retentionTime - normSub)/normDiv;
      }
    } else {
      if (VERB>1) cerr << "Approximating retention time with scan number." << endl;
      double minRT = (double) psm->scan, diffRT = psms.rbegin()->scan - psm->scan;
      normDiv=diffRT/2.;
      normSub=minRT+normDiv;
      if (normDiv==0.0) normDiv = 1.0;
      for(; psm != psms.end(); ++psm) {
        psm->retentionTime = ((double) psm->scan - normSub)/normDiv;
      }
	}
  } else {
    double minRT = scan2rt.begin()->second, diffRT = scan2rt.rbegin()->second - minRT;
    normDiv=diffRT/2.;
    normSub=minRT+normDiv;
    if (normDiv==0.0) normDiv = 1.0;
    for(; psm != psms.end(); ++psm) {
      assert(scan2rt.count(psm->scan)>0);
      psm->retentionTime = (scan2rt[psm->scan] - normSub)/normDiv;
    }
  }
}

void PSMDescription::setPSMSet(vector<PSMDescription> & psms)
{
	double minRT = 1e10, maxRT = -1;
	vector<PSMDescription>::iterator psm;

	for(psm = psms.begin(); psm != psms.end(); ++psm)
	{
		minRT = min(minRT, psm->retentionTime);
	    maxRT = max(maxRT, psm->retentionTime);
	}

	normDiv = (maxRT - minRT) / 2.;
	normSub = minRT + normDiv;
	if (normDiv == 0.0)
		normDiv = 1.0;
}

vector<double *> PSMDescription::getRetFeatures(vector<PSMDescription> & psms)
{
	vector<double*> features;
	for(int i = 0; i < psms.size(); ++i)
		features.push_back(psms[i].getRetentionFeatures());

	return features;
}

void PSMDescription::normalizeRetentionTimes(vector<PSMDescription> & psms)
{
	vector<PSMDescription>::iterator psm;

	for(psm = psms.begin(); psm != psms.end(); ++psm)
		psm->retentionTime = (psm->retentionTime - normSub)/normDiv;
}

