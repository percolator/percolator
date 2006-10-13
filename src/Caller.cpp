#include <iostream>
#include <fstream>
#include <getopt.h>
#include <vector>
#include <string>
using namespace std;
#include "DataSet.h"
#include "Normalizer.h"
#include "Scores.h"
#include "Normalizer.h"
#include "IsoChargeSet.h"
#include "Caller.h"
#include "ssl.h"

Caller::Caller()
{
  forwardFN = "/var/noble/data/tandem-ms/maccoss/rphase/2006-05-06/forward/050606-pps-2-01-forward-no-norm.sqt";
  shuffledFN = "/var/noble/data/tandem-ms/maccoss/rphase/2006-05-06/random/050606-pps-2-01-random-no-norm.sqt";
  shuffled2FN = "/var/noble/data/tandem-ms/maccoss/rphase/2006-05-06/random2/050606-pps-2-01-random2-nonorm.sqt";
  doRoc = false;
  rocFN = "";
  gistFN = "";
  weightFN = "";
  DataSet::setQuadraticFeatures(false);
  fdr=0.01;
}

Caller::~Caller()
{
}

bool Caller::parseOptions(int argc, char **argv){
   int c;
   opterr = 0;
 
   while ((c = getopt (argc, argv, "f:s:S:r:ug:w:qF:")) != -1)
   switch (c)
   {
   case 'f':
     forwardFN = optarg;
     break;
   case 's':
     shuffledFN = optarg;
     break;
   case 'S':
     shuffled2FN = optarg;
     break;
   case 'g':
     gistFN = optarg;
     break;
   case 'w':
     weightFN = optarg;
     break;
   case 'r':
     doRoc = true;
     rocFN = optarg;
     break;
   case 'u':
     Normalizer::setType(Normalizer::UNI);
     break;
   case 'q':
     DataSet::setQuadraticFeatures(true);
     break;
   case 'F':
     fdr=atof(optarg);
     break;
   case '?':
     if (isprint (optopt))
	   fprintf (stderr, "Unknown option `-%c'.\n", optopt);
     else
       fprintf (stderr, "Unknown option character `\\x%x'.\n",optopt);
     return false;
   default:
     abort ();
   }
     
/*       for (index = optind; index < argc; index++)
         printf ("Non-option argument %s\n", argv[index]);
       return 0;
*/
  return true;
}
void Caller::step(double *w) {
  	scores.calcScores(w,*pSet);
  	vector<int> ixs;
  	vector<int>::iterator it;
  	w[DataSet::getNumFeatures()]=scores.getPositiveTrainingIxs(fdr,ixs);
  	int rows = ixs.size();
  	int negatives = pSet->getSubSet(1)->getIsoChargeSize(100);
  	rows += negatives;
  	int nz = (DataSet::getNumFeatures()+1)*rows;
  	double* VAL = new double[nz];
    int* C = new int[nz];
    int* R = new int[rows+1];
  	int r=0;
  	int pos=0;
  	Normalizer *norm =pSet->getNormalizer();
  	double * featUnl = pSet->getSubSet(0)->getFeature();
  	for(it=ixs.begin();it!=ixs.end();it++) {
  	  R[r++]=pos;
  	  norm->normalize(&featUnl[DataSet::rowIx(*it)],&VAL[pos]);
  	  for(int i=0;i<=DataSet::getNumFeatures();i++,pos++) {
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
  	  for(int i=0;i<=DataSet::getNumFeatures();i++,pos++) {
  	    C[pos]=i;
  	  } 
  	  VAL[pos-1]=1.0;        	  
  	}
  	R[r]=nz;
/*  	int ixy;
  	for (ixy=0;ixy<20;ixy++){
  	  cout << ixy << " " << VAL[ixy] << " " << C[ixy] << " " << R[ixy] << "\n";
  	}
  	for (ixy=nz-3;ixy<nz;ixy++){
  	  cout << ixy << " " << VAL[ixy] << " " << C[ixy] << " " << R[ixy-nz+rows] << "\n";
  	} */
  	cout << "Calling with " << positives << " positives and " << negatives << " negatives\n";
  	struct data *Data = new data;
  	Data->n=DataSet::getNumFeatures()+1;
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
    struct options *Options = new options;
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
    
    struct vector_double *Weights = new vector_double;
    Weights->d = DataSet::getNumFeatures()+1;
    Weights->vec = new double[Weights->d];
    for(int ix=0;ix<Weights->d;ix++) Weights->vec[ix]=0;
    
    struct vector_double *Outputs = new vector_double;
    Outputs->vec = new double[rows];
    Outputs->d = rows;
    for(int ix=0;ix<Outputs->d;ix++) Outputs->vec[ix]=0;

//    norm->normalizeweight(w,Weights->vec);
//    Weights->vec[DataSet::getNumFeatures()] = 0;
    L2_SVM_MFN(Data,Options,Weights,Outputs,0);
  	norm->unnormalizeweight(Weights->vec,w);

  	Clear(Data);
  	delete [] Weights->vec;
  	delete Weights;
  	delete [] Outputs->vec;
  	delete Outputs;
    delete Options;
}


int Caller::run() {
  bool doShuffled2 = shuffled2FN.size()>0;
  DataSet forward,shuffled,shuffled2;
  forward.read_sqt(forwardFN);
  forward.setLabel(0);
  shuffled.read_sqt(shuffledFN);
  shuffled.setLabel(-1);
  if (doShuffled2) {
    shuffled2.read_sqt(shuffled2FN);
    shuffled2.setLabel(-1);
  }
  
  vector<DataSet *> set1(2); 
  set1[0]=&forward;
  set1[1]=&shuffled;

  vector<DataSet *> set2(2); 
  set2[0]=&forward;
  set2[1]=(doShuffled2?&shuffled2:&shuffled);
  
  IsoChargeSet train,test;
  train.setSet(&set1);
  test.setSet(&set2);
  setSet(&train);
  if (gistFN.length()>0) {
    pSet->gistWrite(gistFN);
  }
  double w[DataSet::getNumFeatures()+1];
  for(int ix=0;ix<DataSet::getNumFeatures();ix++) w[ix]=0;
  w[3]=1;
  w[DataSet::getNumFeatures()]=1;
  for(int i=0;i<20;i++) {
  	step(w);
  }
  if (weightFN.length()>0) {
     ofstream weightStream(weightFN.data(),ios::out);
     for(int ix=0;ix<DataSet::getNumFeatures()+1;ix++) {
        weightStream << w[ix] << " ";
     }
     weightStream << endl;
     weightStream.close(); 
  }
  if (doRoc) {
    if(doShuffled2) {
      Scores testScores;
      testScores.calcScores(w,test);
      testScores.printRoc(rocFN);
    } else {
      scores.printRoc(rocFN);
    }
  }
  
  return 0;
}

int main(int argc, char **argv){
  Caller caller;
  if(caller.parseOptions(argc,argv))
    return caller.run();
  return 1;
}	






