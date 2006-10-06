#include "Caller.h"
#include "ssl.h"

Caller::Caller()
{
}

Caller::~Caller()
{
}

void Caller::step(double *w) {
  	scores.calcScores(w,all);
  	vector<int> ixs;
  	vector<int>::iterator it;
  	scores.getPositiveTrainingIxs(0.01,ixs);
  	int rows = ixs.getSize();
  	rows += pSet->getSubSetSize(1);
  	int nz = (NUM_FEATURE+1)*rows;
  	double* VAL = new double[nz];
    int* C = new int[nz];
    int* R = new int[NUM_FEATURE+1];
  	int s=0;
  	int pos=0;
  	double ** featUnl = pSet->getSubSet(0)->getFeature();
  	for(it=ixs.begin();it!=ixs.end;it++) {
  	  double * f = featUnl[*it];
  	  for(int i=0;i<NUM_FEATURE;i++,pos++) {
  	    VAL[pos]=f[ix];
  	  } 
  	  VAL[pos++]=1;
  	}
  	struct data *Data = new data[1];
  	
}

int main() {
  char *forwardFN = "/var/noble/data/tandem-ms/maccoss/rphase/2006-05-06/forward/050606-pps-2-01-forward-no-norm.sqt";
  char *randomFN = "/var/noble/data/tandem-ms/maccoss/rphase/2006-05-06/random/050606-pps-2-01-random-no-norm.sqt";
//  char *random2FN = "/var/noble/data/tandem-ms/maccoss/rphase/2006-05-06/random2/050606-pps-2-01-random2-nonorm.sqt";

  cout << endl << " ReadData\n";
  vector<DataSet> sets(2); 
  set[0].read_sqt(forwardFN);
  set[0].setLabel(0);
  set[1].read_sqt(randomFN);
  set[1].setLabel(-1);
  
  IsoChargeSet all(2);
  all.setSet(&sets);
  Caller caller();
  caller.setSet(&all)
  double w[9];
  w[3]=1;
  w[8]=1;
  for(int i=0;i<10;i++) {
  	caller.step(w);
  }
}
