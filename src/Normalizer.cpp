#include <math.h>
#include "DataSet.h"
#include "Normalizer.h"
#include "IsoChargeSet.h"

Normalizer::Normalizer()
{
  avg.resize(NUM_FEATURE,0.0);
  stdv.resize(NUM_FEATURE,0.0);
}

Normalizer::~Normalizer()
{
}

void Normalizer::normalize(const double *in,double* out){
  for (int ix=0;ix<NUM_FEATURE;ix++) {
  	out[ix]=(in[ix]-avg[ix])/stdv[ix];
  }
}

void Normalizer::unnormalizeweight(const double *in,double* out){
  double sum = 0;
  int i=0;
  for (;i<NUM_FEATURE;i++) {
  	out[i]=in[i]/stdv[i];
  	sum=avg[i]/stdv[i];
  }
  out[i]=in[i]-sum;
}

void Normalizer::setSet(IsoChargeSet * set){
  int n=0;
  int setPos=0;
  int ixPos=-1;
  const double * features;
  while((features=set->getNext(setPos,ixPos))!=NULL) {
	n++;
	for (int ix=0;ix<NUM_FEATURE;ix++) {
	  avg[ix]+=features[ix];
	}
  }
  for (int ix=0;ix<NUM_FEATURE;ix++) {
  	avg[ix]/=n;
  }
  setPos=0;
  ixPos=-1;
  while((features=set->getNext(setPos,ixPos))!=NULL) {
    for (int ix=0;ix<NUM_FEATURE;ix++) {
      double d = features[ix]-avg[ix];
      stdv[ix]+=d*d;
    }
  }
  for (int ix=0;ix<NUM_FEATURE;ix++) {
  	stdv[ix]=sqrt(stdv[ix]/n);
  }
}
