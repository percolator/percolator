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
  forwardFN = "";
  shuffledFN = "";
  shuffled2FN = "";
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
  char * host = getenv("HOST");
  oss << greeter();
  oss << "Issued command:" << endl << call;
  oss << "Started " << ctime(&startTime);
  oss.seekp(-1, ios_base::cur);
  oss << " on " << host << endl;
  oss << "Hyperparameters fdr=" << fdr;
  oss << ", Cpos=" << Cpos << ", Cneg=" << Cneg << ", maxNiter=" << niter << endl;
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
  intro << "   percolator [options] forward shuffled [shuffled2]" << endl << endl;
  intro << "   where forward is the normal sqt-file," << endl;
  intro << "         shuffle the shuffled sqt-file," << endl;
  intro << "         and shuffle2 is a possible second shuffled sqt-file for validation" << endl;
  // init
  ArgvParser cmd;
  cmd.setIntroductoryDescription(intro.str());
  cmd.defineOption("o",
    "Create an SQT file with the given name, in which the XCorr value has been replaced with the learned score.",
    ArgvParser::OptionRequiresValue);
  cmd.defineOptionAlternative("o","sqt-out");
  cmd.defineOption("p",
    "Cpos, penalizing factor for misstakes made on positive examples. Default is 10",
    ArgvParser::OptionRequiresValue);
  cmd.defineOption("n",
    "Cneg, penalizing factor for misstakes made on negative examples. Default is 10",
    ArgvParser::OptionRequiresValue);
  cmd.defineOption("F",
    "Use the specified false discovery rate threshold to define positive examples.  Default=0.01.",
    ArgvParser::OptionRequiresValue);
  cmd.defineOptionAlternative("F","fdr");
  cmd.defineOption("i",
    "Maximal number of iteratins",
    ArgvParser::OptionRequiresValue);
  cmd.defineOptionAlternative("i","maxiter");
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
  cmd.defineOption("y",
    "Turn off calculation of tryptic features");
  cmd.defineOption("aaron",
    "Debuging option: use Aarons files");

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
  if (cmd.foundOption("y"))
    DataSet::setTrypticFeatures(false);
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

  if (cmd.foundOption("aaron")) {
    forwardFN = "/var/noble/data/tandem-ms/maccoss/rphase/2006-05-06/forward/050606-pps-2-01-forward-no-norm.sqt";
    shuffledFN = "/var/noble/data/tandem-ms/maccoss/rphase/2006-05-06/random/050606-pps-2-01-random-no-norm.sqt";
    shuffled2FN = "/var/noble/data/tandem-ms/maccoss/rphase/2006-05-06/random2/050606-pps-2-01-random2-nonorm.sqt";
  }
  return true;
}

void Caller::readFile(const string fn, const int label, vector<DataSet *> & sets) {
  ifstream fileIn(fn.data(),ios::in);
  if (!fileIn) {
    cerr << "Could not open file " << fn << endl;
    exit(-1);
  }
  string line;
  if (!getline(fileIn,line)) {
    cerr << "Could not read file " << fn << endl;
    exit(-1);
  }
  fileIn.close();
  if (line.size()>1 && line[0] == 'H' && line[1]=='\t') {
    if (line.find("SQTGenerator")==string::npos) {
      cerr << "SQT file not generated by SEQUEST: " << fn << endl;
      exit(-1);  
    }
    DataSet * pSet = new DataSet();
    pSet->setLabel(label);
    pSet->read_sqt(fn);
    sets.push_back(pSet);
  } else {
    // we hopefully found a meta file
    ifstream meta(fn.data(),ios::in);
    while(getline(fileIn,line)) {
      if (line.size()>0 && line[0] != '#')
        readFile(line,label,sets);
    }
  }
}

void Caller::modifyFile(const string fn, vector<DataSet *> & sets, Scores &sc , const string greet) {
  string line;
  ifstream fileIn(fn.data(),ios::in);
  if (sets.size()>1 && (!fileIn)) {
    cerr << "More than one input file, and file " << fn << 
            " do not contain listing to a set of files" << endl;
    exit(-1);
  }
  if (sets.size()>1 && getline(fileIn,line) && line.size()>1 && 
        line[0] == 'H' && line[1]=='\t') {
    cerr << "More than one input file, and file " << fn << 
            " do not contain listing to a set of files" << endl;
    exit(-1);
  }
  if (!!fileIn)
    fileIn.close();
  if (sets.size()==1 ) {
    vector<double> s,q;
    s.resize(sets[0]->getSize(),-100);
    q.resize(sets[0]->getSize(),-100);
    sc.getScoreAndQ(0,s,q);
    sets[0]->modify_sqt(fn,s,q,greet);
    return;
  }
  unsigned int ix=0;
  fileIn.open(fn.data(),ios::in);
  while(getline(fileIn,line)) {
    if(line.size()>0 && line[0]!='#') {
      vector<double> s,q;
      s.resize(sets[ix]->getSize(),-100);
      q.resize(sets[ix]->getSize(),-100);
      sc.getScoreAndQ(ix,s,q);
      sets[ix++]->modify_sqt(fn,s,q,greet);
    }    
  }
  fileIn.close();
}


void clean(vector<DataSet *> & set) {
  for(unsigned int ix=0;ix<set.size();ix++) {
    if (set[ix]!=NULL) 
      delete set[ix];
    set[ix]=NULL;
  }
}

void Caller::step(double *w) {
  	scores.calcScores(w,*pSet,fdr);
  	vector<int> pos_ix;
    vector<int> pos_set;
    scores.getPositiveTrainingIxs(fdr,pos_set,pos_ix);
  	int rows = pos_ix.size();
  	int negatives = pSet->getNegativeSize();
  	rows += negatives;
  	int nz = (DataSet::getNumFeatures()+1)*rows;
  	double* VAL = new double[nz];
    int* C = new int[nz];
    int* R = new int[rows+1];
  	int r=0;
  	int pos=0;
  	Normalizer *norm =pSet->getNormalizer();
  	for(unsigned int ix=0;ix<pos_ix.size();ix++) {
  	  R[r++]=pos;
  	  norm->normalize(pSet->getFeatures(pos_set[ix],pos_ix[ix]),&VAL[pos]);
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
    int set =-1;
    while(pSet->getLabel(++set)!=-1);
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
  cout << extendedGreeter();
  bool doShuffled2 = shuffled2FN.size()>0;
  vector<DataSet *> forward,shuffled,shuffled2;
  readFile(forwardFN,1,forward);
  readFile(shuffledFN,-1,shuffled);
  if (doShuffled2)
    readFile(shuffled2FN,-1,shuffled2);
  
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
  double diff = difftime (end,startTime);
  
  ostringstream timerValues;
  timerValues.precision(4);
  timerValues << "Processing took " << ((double)(clock()-startClock))/(double)CLOCKS_PER_SEC;
  timerValues << " cpu seconds or " << diff << " seconds wall time" << endl; 
  cout << timerValues.str();
  Scores testScores;
  testScores.calcScores(w,test);
  if (modifiedFN.size()>0) {
    modifyFile(modifiedFN,forward,testScores,extendedGreeter()+timerValues.str());
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
  clean(forward);clean(shuffled);clean(shuffled2);
  return 0;
}

int main(int argc, char **argv){
  Caller caller;
  if(caller.parseOptions(argc,argv))
    return caller.run();
  return -1;
}	






