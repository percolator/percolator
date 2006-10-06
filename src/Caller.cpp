#include <iostream>
#include "DataSet.h"
#include "Normalizer.h"
#include "IsoChargeSet.h"
#include "Caller.h"
#include "ssl.h"

Caller::Caller()
{
}

Caller::~Caller()
{
}


void Caller::step(double *w) {
  	scores.calcScores(w,*pSet);
  	vector<int> ixs;
  	vector<int>::iterator it;
  	scores.getPositiveTrainingIxs(0.01,ixs);
  	int rows = ixs.size();
  	int negatives = pSet->getSubSet(1)->getIsoChargeSize(pSet->getCharge());
  	rows += negatives;
  	int nz = (NUM_FEATURE+1)*rows;
  	double* VAL = new double[nz];
    int* C = new int[nz];
    int* R = new int[rows];
  	int r=0;
  	int pos=0;
  	Normalizer *norm =pSet->getNormalizer();
  	double * featUnl = pSet->getSubSet(0)->getFeature();
  	for(it=ixs.begin();it!=ixs.end();it++) {
  	  R[r++]=pos;
  	  norm->normalize(&featUnl[ROW(*it)],&VAL[pos]);
  	  for(int i=0;i<=NUM_FEATURE;i++,pos++) {
  	    C[pos]=i;
  	  } 
  	  VAL[pos-1]=1.0;        	  
  	}
  	int positives = r;
//  	double * featNeg = pSet->getSubSet(1)->getFeature();
//  	getNext...
//  	int negatives = pSet->getSubSet(1)->getSize();
    int rowInSet = -1;
    int set =1;
    while(const double * featNeg = pSet->getNext(set,rowInSet)) {
  	  R[r++]=pos;
  	  norm->normalize(featNeg,&VAL[pos]);
  	  for(int i=0;i<=NUM_FEATURE;i++,pos++) {
  	    C[pos]=i;
  	  } 
  	  VAL[pos-1]=1.0;        	  
  	}
  	cout << "Calling with " << positives << " positives and " << negatives << " negatives\n";
  	struct data *Data = new data[1];
  	Data->n=NUM_FEATURE+1;
    Data->m=rows;
    Data->l=rows;
    Data->u=0;
    Data->val=VAL;
    Data->rowptr=R;
    Data->colind=C;
    Data->nz = nz;
    Data->Y = new double[rows];
    Data->C = new double[rows];
  	for(int a=0;a<positives;a++) {
  		Data->Y[a]=1;
  		Data->C[a]=10;  		
  	}
  	for(int b=positives;b<negatives+positives;b++) {
  		Data->Y[b]=-1;
  		Data->C[b]=10;  		
  	}
  	// Setup options
    struct options *Options = new options[1];
  	Options->algo = 1;
    Options->lambda=1.0;
    Options->lambda_u=1.0;
    Options->S=10000;
    Options->R=0.5;
    Options->epsilon=EPSILON;
    Options->cgitermax=CGITERMAX;
    Options->mfnitermax=MFNITERMAX;
    Options->Cp = 10;
    Options->Cn = 10;
    
    struct vector_double *Weights = new vector_double[1];
    Weights->d = NUM_FEATURE+1;
    Weights->vec = new double[Weights->d];
    struct vector_double *Outputs = new vector_double[1];
  	ssl_train(Data, Options,Weights,Outputs);
  	norm->unnormalizeweight(Weights->vec,w);
}

int main(int argc, char **argv){
  char *forwardFN = "/var/noble/data/tandem-ms/maccoss/rphase/2006-05-06/forward/050606-pps-2-01-forward-no-norm.sqt";
  char *randomFN = "/var/noble/data/tandem-ms/maccoss/rphase/2006-05-06/random/050606-pps-2-01-random-no-norm.sqt";
//  char *random2FN = "/var/noble/data/tandem-ms/maccoss/rphase/2006-05-06/random2/050606-pps-2-01-random2-nonorm.sqt";

  cout << endl << " ReadData\n";
  vector<DataSet> set(2); 
  set[0].read_sqt(forwardFN);
  set[0].setLabel(0);
  set[1].read_sqt(randomFN);
  set[1].setLabel(-1);
  
  IsoChargeSet all(2);
  all.setSet(&set);
  Caller caller;
  caller.setSet(&all);
  double w[9];
  w[3]=1;
  w[8]=1;
  for(int i=0;i<10;i++) {
  	caller.step(w);
  }
}
