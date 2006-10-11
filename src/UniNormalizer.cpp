#include <vector>
#include <iostream>
#include <math.h>
using namespace std;
#include "DataSet.h"
#include "Normalizer.h"
#include "UniNormalizer.h"
#include "IsoChargeSet.h"

UniNormalizer::UniNormalizer()
{
  sub.resize(NUM_FEATURE,0.0);
  div.resize(NUM_FEATURE,0.0);
}

UniNormalizer::~UniNormalizer()
{
}

void UniNormalizer::normalize(const double *in,double* out){
  for (int ix=0;ix<NUM_FEATURE;ix++) {
  	out[ix]=(in[ix]-sub[ix])/div[ix];
  }
}

void UniNormalizer::unnormalizeweight(const double *in,double* out){
  double sum = 0;
  int i=0;
  for (;i<NUM_FEATURE;i++) {
  	out[i]=in[i]/div[i];
  	sum=sub[i]/div[i];
  }
  out[i]=in[i]-sum;
}

void UniNormalizer::normalizeweight(const double *in,double* out){
  double sum = 0;
  int i=0;
  for (;i<NUM_FEATURE;i++) {
  	out[i]=in[i]*div[i];
  	sum=sub[i]/div[i];
  }
  out[i]=in[i]+sum;
}

void UniNormalizer::setSet(IsoChargeSet * set){
  int setPos=0;
  int ixPos=-1;
  const double * features;
  int ix;
  vector<double> mins(NUM_FEATURE,1e+100);
  vector<double> maxs(NUM_FEATURE,-1e+100);
  while((features=set->getNext(setPos,ixPos))!=NULL) {
	for (ix=0;ix<NUM_FEATURE;ix++) {
	  mins[ix]=min(features[ix],mins[ix]);
	  maxs[ix]=max(features[ix],maxs[ix]);
	}
  }
  for (ix=0;ix<NUM_FEATURE;ix++) {
  	sub[ix]=mins[ix];
  	div[ix]=maxs[ix]-mins[ix];
  	if (div[ix]<=0)
  	  div[ix]=1.0;
  }
}
