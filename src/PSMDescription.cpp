#include <assert.h>
#include "Globals.h"
#include "PSMDescription.h"
#include "DescriptionOfCorrect.h"

PSMDescription::PSMDescription():  q(0.),pep(0.),features(NULL),retentionFeatures(NULL),
                                   retentionTime(0.),predictedTime(0.),massDiff(0.),pI(0.),
                                   scan(0),id(""),peptide("")
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
        maxRT=min(maxRT,psm->retentionTime);
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
