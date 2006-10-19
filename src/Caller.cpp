#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>
#include <cstdlib>
#include <vector>
#include <string>
using namespace std;
#include "argvparser.h"
using namespace CommandLineProcessing;
#include "DataSet.h"
#include "Normalizer.h"
#include "Scores.h"
#include "Normalizer.h"
#include "SetHandler.h"
#include "Caller.h"
#include "ssl.h"

Caller::Caller()
{
  forwardFN = "/var/noble/data/tandem-ms/maccoss/rphase/2006-05-06/forward/050606-pps-2-01-forward-no-norm.sqt";
  shuffledFN = "/var/noble/data/tandem-ms/maccoss/rphase/2006-05-06/random/050606-pps-2-01-random-no-norm.sqt";
  shuffled2FN = "/var/noble/data/tandem-ms/maccoss/rphase/2006-05-06/random2/050606-pps-2-01-random2-nonorm.sqt";
  modifiedFN = "";
  rocFN = "";
  gistFN = "";
  weightFN = "";
  fdr=0.01;
  niter = 10;
  Cpos=10;
  Cneg=10;
}

Caller::~Caller()
{
}

string Caller::extendedGreeter() {
  ostringstream oss;
  oss << greeter();
  oss << "Issued command:" << endl << call;
  oss << "Started " << ctime(&startTime) << " on " << getenv("HOST") << endl;
  oss << 
  oss << "fdr=" << fdr << " (defining positive set)" << endl;
  oss << "Cpos=" << Cpos << " ,Cneg=" << Cneg << endl;
  oss << "maxNiter=" << niter << endl;
  return oss.str();
}

string Caller::greeter() {
  ostringstream oss;
  oss << "percolator (c) 2006 Lukas Käll, University of Washington" << endl;
  oss << "Version  " << __DATE__ << " " << __TIME__ << endl;
  return oss.str();
}

bool Caller::parseOptions(int argc, char **argv){
  ostringstream callStream;
  callStream << argv[0];
  for (int i=1;i<argc;i++) callStream << " " << argv[i];
  callStream << endl;
  call = callStream.str();
  ostringstream intro;
  intro << greeter() << endl << "Usage:" << endl;
  intro << "   percolator [-huq] [-g trunc_fn] [-F val] [-i val] [-w fn] \\" << endl;
  intro << "           [-r fn] [-o sqt_fn] forward shuffled [shuffled2]" << endl << endl;
  intro << "   where forward is the normal sqt-file," << endl;
  intro << "         shuffle the shuffled sqt-file," << endl;
  intro << "         and shuffle2 is a possible second shuffled sqt-file for validation" << endl;
  // init
  ArgvParser cmd;
  cmd.setIntroductoryDescription(intro.str());
  cmd.defineOption("o",
    "Output predictions to a modified sqt-file",
    ArgvParser::OptionRequiresValue);
  cmd.defineOption("p",
    "Cpos, penalizing factor for misstakes made on positive examples. Default is 10",
    ArgvParser::OptionRequiresValue);
  cmd.defineOption("n",
    "Cneg, penalizing factor for misstakes made on negative examples. Default is 10",
    ArgvParser::OptionRequiresValue);
  cmd.defineOption("F",
    "The FDR filter value. Default is 0.01",
    ArgvParser::OptionRequiresValue);
  cmd.defineOption("i",
    "Maximal number of iteratins",
    ArgvParser::OptionRequiresValue);
  cmd.defineOption("g",
    "Output test features to a Gist format files with the given trunc filename",
    ArgvParser::OptionRequiresValue);
  cmd.defineOption("w",
    "Output final weights to the given filename",
    ArgvParser::OptionRequiresValue);
  cmd.defineOption("r",
    "Output result file (score ranked labels) to given filename",
    ArgvParser::OptionRequiresValue);
  cmd.defineOption("u",
    "Use unit normalization [0-1] instead of standard deviation normalization");
  cmd.defineOption("q",
    "Calculate quadratic feature terms");

  //define error codes
  cmd.addErrorCode(0, "Success");
  cmd.addErrorCode(-1, "Error");

  cmd.setHelpOption("h", "help", "Print this help page");

  // finally parse and handle return codes (display help etc...)
  int result = cmd.parse(argc, argv);

  if (result != ArgvParser::NoParserError)
  {
    cout << cmd.parseErrorDescription(result);
    exit(-1);
  }
  
  // now query the parsing results
  if (cmd.foundOption("o"))
    modifiedFN = cmd.optionValue("o");
  if (cmd.foundOption("p")) {
    Cpos = atof(cmd.optionValue("p").c_str());
    if (Cpos<=0.0 || Cpos > 1e127) {
      cerr << "-p option requres a positive float value" << endl;
      cerr << cmd.usageDescription();
      exit(-1); 
    }
  }
  if (cmd.foundOption("n")) {
    Cneg = atof(cmd.optionValue("n").c_str());
    if (Cneg<=0.0 || Cneg > 1e127) {
      cerr << "-n option requres a positive float value" << endl;
      cerr << cmd.usageDescription();
      exit(-1); 
    }
  }
  if (cmd.foundOption("g"))
    gistFN = cmd.optionValue("g");
  if (cmd.foundOption("w"))
    weightFN = cmd.optionValue("w");
  if (cmd.foundOption("r"))
    rocFN = cmd.optionValue("r");
  if (cmd.foundOption("u"))
    Normalizer::setType(Normalizer::UNI);
  if (cmd.foundOption("q"))
    DataSet::setQuadraticFeatures(true);
  if (cmd.foundOption("i")) {
    niter = atoi(cmd.optionValue("n").c_str());
  }
  if (cmd.foundOption("F")) {
    fdr = atof(cmd.optionValue("F").c_str());
    if (fdr<=0.0 || fdr > 1.0) {
      cerr << "-F option requres a positive number < 1" << endl;
      cerr << cmd.usageDescription();
      exit(-1); 
    }
  }
  if (cmd.arguments()>3) {
      cerr << "Too many arguments given" << endl;
      cerr << cmd.usageDescription();
      exit(1);   
  }
  if (cmd.arguments()>0)
    forwardFN = cmd.argument(0);
  if (cmd.arguments()>1)
     shuffledFN = cmd.argument(1);
  if (cmd.arguments()>2)
     shuffled2FN = cmd.argument(2);
  return true;
}

void Caller::step(double *w) {
  	scores.calcScores(w,*pSet);
  	vector<int> ixs;
  	vector<int>::iterator it;
  	w[DataSet::getNumFeatures()]=scores.getPositiveTrainingIxs(fdr,ixs);
  	int rows = ixs.size();
  	int negatives = pSet->getSubSet(1)->getSize();
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
  		Data->C[a]=Cpos;  		
  	}
  	for(int b=positives;b<negatives+positives;b++) {
  		Data->Y[b]=-1;
  		Data->C[b]=Cneg;  		
  	}
  	// Setup options
    struct options *Options = new options;
    Options->lambda=1.0;
    Options->lambda_u=1.0;
    Options->epsilon=EPSILON;
    Options->cgitermax=CGITERMAX;
    Options->mfnitermax=MFNITERMAX;
    
    struct vector_double *Weights = new vector_double;
    Weights->d = DataSet::getNumFeatures()+1;
    Weights->vec = new double[Weights->d];
//    for(int ix=0;ix<Weights->d;ix++) Weights->vec[ix]=w[ix];
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
  time(&startTime);
  startClock=clock();
  cerr << greeter();
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
  
  
  SetHandler train,test;
  train.setSet(forward,shuffled);
  test.setSet(forward,(doShuffled2?shuffled2:shuffled));
  setSet(&train);
  if (gistFN.length()>0) {
    pSet->gistWrite(gistFN);
  }
  double w[DataSet::getNumFeatures()+1];
  for(int ix=0;ix<DataSet::getNumFeatures();ix++) w[ix]=0;
  w[3]=1;
  w[DataSet::getNumFeatures()]=1;
  for(int i=0;i<niter;i++) {
    cout << "Iteration " << i+1 << " : ";
  	step(w);
  }
  time_t end;
  time (&end);
  double dif = difftime (end,startTime);
  
  ostringstream timerValues;
  timerValues << "Processing took " << (clock()-startClock)/CLOCKS_PER_SEC;
  timerValues << " cpu seconds or " << dif << " wall time seconds" << endl; 
  cout << timerValues.str();
  Scores testScores;
  testScores.calcScores(w,test);
  if (modifiedFN.size()>0) {
    vector<double> sc(forward.getSize(),0.0),fdr(forward.getSize(),0.0);
    testScores.getScoreAndFdr(0,sc,fdr);
    forward.modify_sqt(modifiedFN,sc,fdr,extendedGreeter()+timerValues.str());
  }
  if (weightFN.size()>0) {
     ofstream weightStream(weightFN.data(),ios::out);
     for(int ix=0;ix<DataSet::getNumFeatures()+1;ix++) {
        weightStream << w[ix] << " ";
     }
     weightStream << endl;
     weightStream.close(); 
  }
  if (rocFN.size()>0) {
      testScores.printRoc(rocFN);
  }
  
  return 0;
}

int main(int argc, char **argv){
  Caller caller;
  if(caller.parseOptions(argc,argv))
    return caller.run();
  return -1;
}	






