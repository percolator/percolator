#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>
#include <cstdlib>
#include <vector>
#include <set>
#include <map>
#include <string>
using namespace std;
#include "Option.h"
#include "DataSet.h"
#include "IntraSetRelation.h"
#include "Normalizer.h"
#include "Scores.h"
#include "Normalizer.h"
#include "SetHandler.h"
#include "ssl.h"
#include "Caller.h"
#include "Globals.h"

const unsigned int Caller::xval_fold = 3;

Caller::Caller()
{
  forwardFN = "";
  shuffledFN = "";
  shuffled2FN = "";
  modifiedFN = "";
  modifiedShuffledFN = "";
  shuffledWC = "";
  rocFN = "";
  gistFN = "";
  weightFN = "";
  selectionfdr=0.01;
  test_fdr=0.01;
  niter = 10;
  selectedCpos=0;
  selectedCneg=0;
  xv_type=EACH_STEP;
  gistInput=false;
  pNorm=NULL;
  svmInput=NULL;
}

Caller::~Caller()
{
    if (pNorm)
      delete pNorm;
    pNorm=NULL;
    if (svmInput)
      delete svmInput;
}

string Caller::extendedGreeter() {
  ostringstream oss;
  char * host = getenv("HOST");
  oss << greeter();
  oss << "Issued command:" << endl << call;
  oss << "Started " << ctime(&startTime);
  oss.seekp(-1, ios_base::cur);
  oss << " on " << host << endl;
  oss << "Hyperparameters fdr=" << selectionfdr;
  oss << ", Cpos=" << selectedCpos << ", Cneg=" << selectedCneg << ", maxNiter=" << niter << endl;
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
  intro << "   percolator [options] forward shuffled [shuffled2]" << endl;
  intro << "or percolator [options] -P pattern mormal_and_shuffled.sqt" << endl;
  intro << "or percolator [options] -g gist.data gist.label" << endl << endl;
  intro << "   where forward is the normal sqt-file," << endl;
  intro << "         shuffle the shuffled sqt-file," << endl;
  intro << "         and shuffle2 is an otional second shuffled sqt-file for q-value calculation" << endl;
  // init
  CommandLineParser cmd(intro.str());
  cmd.defineOption("o","sqt-out",
    "Remake an sqt file out of the normal sqt-file with the given name, \
in which the XCorr value has been replaced with the learned score \
and Sp has been replaced with the negated q-value.","filename");
  cmd.defineOption("s","shuffled",
    "Remake an SQT file out of the test (shuffled or if present shuffled2) sqt-file \
with the given name, \
in which the XCorr value has been replaced with the learned score \
and Sp has been replaced with the negated q-value.",
    "filename");
  cmd.defineOption("P","pattern",
    "Option for single sqt file mode defining the name pattern used for shuffled data base. Typically set to random_seq",
    "pattern");
  cmd.defineOption("p","Cpos",
    "Cpos, penalty for mistakes made on positive examples. Set by cross validation if not specified.",
    "value");
  cmd.defineOption("n","Cneg",
    "Cneg, penalty for mistakes made on negative examples. Set by cross validation if not specified or -p not specified.",
    "value");
  cmd.defineOption("F","trainFDR",
    "False discovery rate threshold to define positive examples in training. Set by cross validation if 0. Default is 0.01.",
    "value");
  cmd.defineOption("t","testFDR",
    "False discovery rate threshold for evaluating best cross validation result and the reported end result. Default is 0.01.",
    "value");
  cmd.defineOption("i","maxiter",
    "Maximal number of iteratins",
    "number");
  cmd.defineOption("m","matches",
    "Maximal number of matches to take in consideration per spectrum when using sqt-files",
    "number");
  cmd.defineOption("r","train-ratio",
    "Fraction of the negative data set to be used as train set when only providing one negative set, remaining examples will be used as test set. Set to 0.7 by default.",
    "value");
  cmd.defineOption("G","gist-out",
    "Output test features to a Gist format files with the given trunc filename",
    "trunc name");
  cmd.defineOption("g","gist-in",
    "Input files are given as gist files, first argument should be a file name of the data file, the second the label file. \
Label 1 is interpreted as positive, -1 negative in train set, -2 negative in test set.","",TRUE_IF_SET);
  cmd.defineOption("w","weights",
    "Output final weights to the given filename",
    "filename");
  cmd.defineOption("v","verbose",
    "Set verbosity of output: 0=no processing info, 5=all, default is 2",
    "level");
  cmd.defineOption("r","result",
    "Output result file (score ranked labels) to given filename",
    "filename");
  cmd.defineOption("u","unitnorm",
    "Use unit normalization [0-1] instead of standard deviation normalization","",TRUE_IF_SET);
  cmd.defineOption("q","quadratic",
    "Calculate quadratic feature terms","",TRUE_IF_SET);
  cmd.defineOption("y","notryptic",
    "Turn off calculation of tryptic/chymo-tryptic features.","",TRUE_IF_SET);
  cmd.defineOption("c","chymo",
    "Replace tryptic features with chymo-tryptic features.","",TRUE_IF_SET);
  cmd.defineOption("x","whole-xval",
    "Select hyper parameter cross validation to be performed on whole itterating procedure, rather than on each iteration step."
    ,"",TRUE_IF_SET);

  // finally parse and handle return codes (display help etc...)
  cmd.parseArgs(argc, argv);

  
  // now query the parsing results
  if (cmd.optionSet("o"))
    modifiedFN = cmd.options["o"];
  if (cmd.optionSet("s"))
    modifiedShuffledFN = cmd.options["s"];
  if (cmd.optionSet("P"))
    shuffledWC = cmd.options["P"];
  if (cmd.optionSet("p")) {
    selectedCpos = cmd.getDouble("p",0.0,1e127);
  }
  if (cmd.optionSet("n")) {
    selectedCneg = cmd.getDouble("n",0.0,1e127);
  }
  if (cmd.optionSet("G"))
    gistFN = cmd.options["G"];
  if (cmd.optionSet("g")) {
    gistInput=true;
    if (cmd.arguments.size()!=2) {
      cerr << "Provide exactly two arguments when using gist input" << endl;
      exit(-1);
    }   
  }
  if (cmd.optionSet("w"))
    weightFN = cmd.options["w"];
  if (cmd.optionSet("f")) {
    double frac = cmd.getDouble("f", 0.0 ,1.0);
    Scores::setTrainRatio(frac);
  }
  if (cmd.optionSet("r"))
    rocFN = cmd.options["r"];
  if (cmd.optionSet("u"))
    Normalizer::setType(Normalizer::UNI);
  if (cmd.optionSet("q"))
    DataSet::setQuadraticFeatures(true);
  if (cmd.optionSet("y"))
    DataSet::setTrypticFeatures(false);
  if (cmd.optionSet("c"))
    DataSet::setChymoTrypticFeatures(true);
  if (cmd.optionSet("x"))
    xv_type=WHOLE;
  if (cmd.optionSet("i")) {
    niter = cmd.getInt("i",0,100000000);
  }
  if (cmd.optionSet("m")) {
    DataSet::setHitsPerSpectrum(cmd.getInt("m",1,30000));
  }
  if (cmd.optionSet("v")) {
    Globals::getInstance()->setVerbose(cmd.getInt("v",0,10));
  }
  if (cmd.optionSet("F")) {
    selectionfdr = cmd.getDouble("F",0.0,1.0);
  }
  if (cmd.optionSet("t")) {
    test_fdr = cmd.getDouble("t",0.0,1.0);
  }
  if (cmd.arguments.size()>3) {
      cerr << "Too many arguments given" << endl;
      cmd.help();
  }
  if (cmd.arguments.size()>0)
    forwardFN = cmd.arguments[0];
  if (cmd.arguments.size()>1)
     shuffledFN = cmd.arguments[1];
  if (cmd.arguments.size()>2)
     shuffled2FN = cmd.arguments[2];
  return true;
}

void Caller::printWeights(ostream & weightStream, double * w) {
  weightStream << "# first line contains normalized weights, second line the raw weights" << endl;  
  weightStream << DataSet::getFeatureNames() << "\tm0" << endl;
  weightStream.precision(3);
  weightStream << w[0];
  for(int ix=1;ix<DataSet::getNumFeatures()+1;ix++) {
    weightStream << "\t" << w[ix];
  }
  weightStream << endl;
  double ww[DataSet::getNumFeatures()+1];
  pNorm->unnormalizeweight(w,ww);
  weightStream << ww[0];
  for(int ix=1;ix<DataSet::getNumFeatures()+1;ix++) {
    weightStream << "\t" << ww[ix];
  }
  weightStream << endl;
}


void Caller::step(Scores& train,double * w, double Cpos, double Cneg, double fdr) {
  	train.calcScores(w,test_fdr);
    train.generateTrainingSet(*svmInput,fdr,Cpos,Cneg);
//    int nex=train.getTrainingSetSize();
    int negatives = train.negSize();
    int positives = train.posNowSize();
  	if (VERB>1) cerr << "Calling with " << positives << " positives and " << negatives << " negatives\n";
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
    Outputs->vec = new double[positives+negatives];
    Outputs->d = positives+negatives;
    for(int ix=0;ix<Outputs->d;ix++) Outputs->vec[ix]=0;

//    norm->normalizeweight(w,Weights->vec);
//    Weights->vec[DataSet::getNumFeatures()] = 0;
    L2_SVM_MFN(*svmInput,Options,Weights,Outputs);
    for(int i= DataSet::getNumFeatures()+1;i--;)
      w[i]=Weights->vec[i];
  	delete [] Weights->vec;
  	delete Weights;
  	delete [] Outputs->vec;
  	delete Outputs;
    delete Options;
}

void Caller::trainEm(double * w) {
  // iterate
  for(int i=0;i<niter;i++) {
    if(VERB>1) cerr << "Iteration " << i+1 << " : ";
    if (xv_type!=EACH_STEP)
      step(trainset,w,selectedCpos,selectedCneg,selectionfdr);
    else
      xvalidate_step(w);
    if(VERB>2) {cerr<<"Obtained weights" << endl;printWeights(cerr,w);}
  }
  if(VERB==2 ) {cerr<<"Obtained weights" << endl;printWeights(cerr,w);}
}

void Caller::xvalidate_step(double *w) {
  Globals::getInstance()->decVerbose();
  int bestTP = 0;
  double best_fdr,best_cpos,best_cneg;
  vector<double>::iterator fdr,cpos,cfrac;
  for(fdr=xv_fdrs.begin();fdr!=xv_fdrs.end();fdr++) {
    for(cpos=xv_cposs.begin();cpos!=xv_cposs.end();cpos++) {
      for(cfrac=xv_cfracs.begin();cfrac!=xv_cfracs.end();cfrac++) {
        if(VERB>1) cerr << "-cross validation with cpos=" << *cpos <<
          ", cfrac=" << *cfrac << ", fdr=" << *fdr << endl;
        int tp=0;
        double ww[DataSet::getNumFeatures()+1];
        for (unsigned int i=0;i<xval_fold;i++) {
          if(VERB>2) cerr << "cross calidation - fold " << i+1 << " out of " << xval_fold << endl;
          for(int ix=0;ix<DataSet::getNumFeatures()+1;ix++) ww[ix]=w[ix];
          step(xv_train[i],ww,*cpos,(*cpos)*(*cfrac),*fdr);
          tp += xv_test[i].calcScores(ww,test_fdr);
          if(VERB>2) cerr << "Cumulative # of positives " << tp << endl;
        }
        if(VERB>1) cerr << "- cross validation found " << tp << " positives over " << test_fdr*100 << "% FDR level" << endl;
        if (tp>bestTP) {
          if(VERB>1) cerr << "Better than previous result, store this" << endl;
          bestTP = tp;
          best_fdr=*fdr;
          best_cpos = *cpos;
          best_cneg = (*cpos)*(*cfrac);          
        }
      }     
    }
  }
  Globals::getInstance()->incVerbose();
  if(VERB>0) cerr << "cross validation found " << bestTP << " positives over " << test_fdr*100 << "% FDR level for hyperparameters Cpos=" << best_cpos
                  << ", Cneg=" << best_cneg << ", fdr=" << best_fdr << endl;
  
  step(trainset,w,best_cpos,best_cneg,best_fdr);
}

void Caller::xvalidate(double *w) {
  Globals::getInstance()->decVerbose();
  int bestTP = 0;
  double ww[DataSet::getNumFeatures()+1],www[DataSet::getNumFeatures()+1];
  vector<double>::iterator fdr,cpos,cfrac;
  for(fdr=xv_fdrs.begin();fdr!=xv_fdrs.end();fdr++) {
    for(cpos=xv_cposs.begin();cpos!=xv_cposs.end();cpos++) {
      for(cfrac=xv_cfracs.begin();cfrac!=xv_cfracs.end();cfrac++) {
        if(VERB>1) cerr << "-cross validation with cpos=" << *cpos <<
          ", cfrac=" << *cfrac << ", fdr=" << *fdr << endl;
        int tp=0;
        for (unsigned int i=0;i<xval_fold;i++) {
          if(VERB>2) cerr << "cross calidation - fold " << i+1 << " out of " << xval_fold << endl;
          for(int ix=0;ix<DataSet::getNumFeatures()+1;ix++) ww[ix]=w[ix];
          for(int k=0;k<niter;k++) {
            step(xv_train[i],ww,*cpos,(*cpos)*(*cfrac),*fdr);
          }
          tp += xv_test[i].calcScores(ww,test_fdr);
          if(VERB>2) cerr << "Cumulative # of positives " << tp << endl;
        }
        if(VERB>1) cerr << "- cross validation found " << tp << " positives over " << test_fdr*100 << "% FDR level" << endl;
        if (tp>bestTP) {
          if(VERB>1) cerr << "Better than previous result, store this" << endl;
          bestTP = tp;
          selectionfdr=*fdr;
          selectedCpos = *cpos;
          selectedCneg = (*cpos)*(*cfrac);          
          for(int ix=0;ix<DataSet::getNumFeatures()+1;ix++) www[ix]=ww[ix];
        }
      }     
    }
  }
  Globals::getInstance()->incVerbose();
  if(VERB>0) cerr << "cross validation found " << bestTP << " positives over " << test_fdr*100 << "% FDR level for hyperparameters Cpos=" << selectedCpos
                  << ", Cneg=" << selectedCneg << ", fdr=" << selectionfdr << endl << "Now train on all data" << endl;
  trainEm(www);
}

int Caller::run() {
  time(&startTime);
  startClock=clock();
  if(VERB>0)  cerr << extendedGreeter();
  bool doShuffled2 = !shuffled2FN.empty();
  if (gistInput) {
    normal.readGist(forwardFN,shuffledFN,1);
    shuffled.readGist(forwardFN,shuffledFN,-1);
    shuffled2.readGist(forwardFN,shuffledFN,-2);
    if (shuffled2.getSize()>0) {
      doShuffled2=true;
    }
  } else if (shuffledWC.empty()) {
    normal.readFile(forwardFN,1);
    shuffled.readFile(shuffledFN,-1);    
    if (doShuffled2)
      shuffled2.readFile(shuffled2FN,-1);
  } else {
    normal.readFile(forwardFN,shuffledWC,false);  
    shuffled.readFile(forwardFN,shuffledWC,true);  
  }
  if (doShuffled2) {
    trainset.fillFeatures(normal,shuffled);
    testset.fillFeatures(normal,shuffled2);
  } else {
  	Scores::fillFeatures(trainset,testset,normal,shuffled);
  }
  if (VERB>1) {
    cerr << "Train set contains " << trainset.posSize() << " positives and " << trainset.negSize() << " negatives." << endl;
    cerr << "Test set contains " << testset.posSize() << " positives and " << testset.negSize() << " negatives." << endl;
  }
  if (gistFN.length()>0) {
    SetHandler::gistWrite(gistFN,normal,shuffled,shuffled2);
  }
  //Normalize features
  set<DataSet *> all;
  all.insert(normal.getSubsets().begin(),normal.getSubsets().end());
  all.insert(shuffled.getSubsets().begin(),shuffled.getSubsets().end());
  if (doShuffled2) 
    all.insert(shuffled2.getSubsets().begin(),shuffled2.getSubsets().end());
  pNorm=Normalizer::getNew();
  pNorm->setSet(all);
  pNorm->normalizeSet(all);
  
  svmInput = new AlgIn(trainset.size(),DataSet::getNumFeatures()+1); // One input set, to be reused multiple times
  
  // Set up a first guess of w
  double w[DataSet::getNumFeatures()+1];
  for(int ix=0;ix<DataSet::getNumFeatures()+1;ix++) w[ix]=0;
  w[3]=1;
//  w[DataSet::getNumFeatures()]=1;
  
  if (selectionfdr<=0 || selectedCpos<=0 || selectedCneg <= 0) {
  	xv_train.resize(xval_fold); xv_test.resize(xval_fold);
  	trainset.createXvalSets(xv_train,xv_test,xval_fold);
    if (selectionfdr > 0) {
      xv_fdrs.push_back(selectionfdr);
    } else {
      xv_fdrs.push_back(0.01);xv_fdrs.push_back(0.03); xv_fdrs.push_back(0.07);
      selectionfdr=test_fdr;
      if(VERB>0) cerr << "selecting fdr by cross validation" << endl;
    }
    if (selectedCpos > 0) {
      xv_cposs.push_back(selectedCpos);
    } else {
      xv_cposs.push_back(10);xv_cposs.push_back(1);xv_cposs.push_back(0.1);
      if(VERB>0) cerr << "selecting cpos by cross validation" << endl;
    }
    if (selectedCpos > 0 && selectedCneg > 0) {
      xv_cfracs.push_back(selectedCneg/selectedCpos);
    } else  {
      xv_cfracs.push_back(1);xv_cfracs.push_back(3);xv_cfracs.push_back(10);
      if(VERB>0) cerr << "selecting cneg by cross validation" << endl;  
    }
  } else {
    xv_type = NO_XV;
  }
  if(VERB>0) cerr << "---Training with Cpos=" << selectedCpos <<
          ", Cneg=" << selectedCneg << ", fdr=" << selectionfdr << endl;
  if (xv_type==WHOLE)
    xvalidate(w);
  else  
    trainEm(w);
  time_t end;
  time (&end);
  double diff = difftime (end,startTime);
  
  ostringstream timerValues;
  timerValues.precision(4);
  timerValues << "Processing took " << ((double)(clock()-startClock))/(double)CLOCKS_PER_SEC;
  timerValues << " cpu seconds or " << diff << " seconds wall time" << endl; 
  if (VERB>1) cerr << timerValues.str();
  int overFDR = testset.calcScores(w,selectionfdr);
  if (VERB>0) cerr << "Found " << overFDR << " peptides scoring over " << selectionfdr*100 << "% FDR level on testset" << endl;
  double ww[DataSet::getNumFeatures()+1];
  pNorm->unnormalizeweight(w,ww);    
  normal.modifyFile(modifiedFN,ww,testset,extendedGreeter()+timerValues.str());
  if (doShuffled2)
    shuffled2.modifyFile(modifiedShuffledFN,ww,testset,extendedGreeter()+timerValues.str());
  else
    shuffled.modifyFile(modifiedShuffledFN,ww,testset,extendedGreeter()+timerValues.str());
  if (weightFN.size()>0) {
     ofstream weightStream(weightFN.data(),ios::out);
     printWeights(weightStream,w);
     weightStream.close(); 
  }
  if (rocFN.size()>0) {
      testset.printRoc(rocFN);
  }
  normal.print(testset);
  return 0;
}

int main(int argc, char **argv){
  Caller *pCaller=new Caller();
  int retVal=-1;
  if(pCaller->parseOptions(argc,argv))
  {
    retVal=pCaller->run();
  }
  delete pCaller;
  return retVal;
}	





