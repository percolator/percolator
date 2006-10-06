#include <math.h>
#include "Normalizer.h"

Normalizer::Normalizer()
{
  avg.resize(NUM_FEATURE,0.0);
  stdv.resize(NUM_FEATURE,0.0);
}

Normalizer::~Normalizer()
{
}

void Normalizer::normalize(double *in,double* out){
  for (int ix=0;ix<NUM_FEATURE;ix++) {
  	out[ix]=(in[ix]-avg[ix])/stdv[ix];
  }
}

void Normalizer::setSets(IsoChargeSet *set){
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
