#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>
#include <cstdlib>
#include <vector>
#include <set>
#include <string>
using namespace std;
#include "argvparser.h"
using namespace CommandLineProcessing;
#include "DataSet.h"
#include "VirtualSet.h"
#include "IntraSetRelation.h"
#include "Normalizer.h"
#include "Scores.h"
#include "Normalizer.h"
#include "SetHandler.h"
#include "Caller.h"
#include "ssl.h"
#include "Globals.h"

Caller::Caller()
{
  forwardFN = "";
  shuffledFN = "";
  shuffled2FN = "";
  modifiedFN = "";
  rocFN = "";
  gistFN = "";
  weightFN = "";
  selectedfdr=0;
  niter = 10;
  selectedCpos=0;
  selectedCneg=0;
  pNorm=NULL;
}

Caller::~Caller()
{
    if (pNorm)
      delete pNorm;
    pNorm=NULL;
}

string Caller::extendedGreeter() {
  ostringstream oss;
  char * host = getenv("HOST");
  oss << greeter();
  oss << "Issued command:" << endl << call;
  oss << "Started " << ctime(&startTime);
  oss.seekp(-1, ios_base::cur);
  oss << " on " << host << endl;
  oss << "Hyperparameters fdr=" << selectedfdr;
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
  cmd.defineOption("v",
    "Set verbosity of output: 0=no processing info, 5=all, default is 2",
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
    selectedCpos = atof(cmd.optionValue("p").c_str());
    if (selectedCpos<=0.0 || selectedCpos > 1e127) {
      cerr << "-p option requres a positive float value" << endl;
      cerr << cmd.usageDescription();
      exit(-1); 
    }
  }
  if (cmd.foundOption("n")) {
    selectedCneg = atof(cmd.optionValue("n").c_str());
    if (selectedCneg<=0.0 || selectedCneg > 1e127) {
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
    niter = atoi(cmd.optionValue("i").c_str());
  }
  if (cmd.foundOption("v")) {
    Globals::getInstance()->setVerbose(atoi(cmd.optionValue("v").c_str()));
  }
  if (cmd.foundOption("F")) {
    selectedfdr = atof(cmd.optionValue("F").c_str());
    if (selectedfdr<=0.0 || selectedfdr > 1.0) {
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


void Caller::readFile(const string fn, const int label, vector<DataSet *> & sets, IntraSetRelation *intra, bool calc) {
  ifstream fileIn(fn.c_str(),ios::in);
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
    pSet->read_sqt(fn,intra);
    sets.push_back(pSet);
  } else {
    // we hopefully found a meta file
    string line2;
    ifstream meta(fn.data(),ios::in);
    while(getline(meta,line2)) {
      if (line2.size()>0 && line2[0] != '#')
//        cout << "0:" << line2[0] << " 1:" << line2[1] << " e:" << line2[line2.size()-1] << endl;
        readFile(line2,label,sets,intra,false);
    }
    meta.close();
  }
  if (calc) {
    vector<DataSet *>::iterator it;
    for (it=sets.begin();it!=sets.end();it++) {
      (*it)->computeIntraSetFeatures();
    }
  }  
}

void Caller::modifyFile(const string fn, vector<DataSet *> & sets, double *w, Scores &sc , const string greet) {
  string line;
  ifstream fileIn(fn.c_str(),ios::in);
  if (sets.size()>1 && (!fileIn)) {
    cerr << "More than one input file, and file " << fn << 
            " that should contain listing to a set of files is not readable" << endl;
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
    sets[0]->modify_sqt(fn,w,&sc,greet);
    return;
  }
  unsigned int ix=0;
  fileIn.open(fn.c_str(),ios::in);
  while(getline(fileIn,line)) {
    if(line.size()>0 && line[0]!='#') {
      sets[ix++]->modify_sqt(line,w,&sc,greet);
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

void Caller::step(SetHandler & train,double * w, double Cpos, double Cneg, double fdr) {
  	scores.calcScores(w,train);
    train.generateTrainingSet(fdr,Cpos,Cneg,scores);
    int nex=train.getTrainingSetSize();
    int negatives = train.getNegativeSize();
    int positives = nex - negatives;
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
    Outputs->vec = new double[nex];
    Outputs->d = nex;
    for(int ix=0;ix<Outputs->d;ix++) Outputs->vec[ix]=0;

//    norm->normalizeweight(w,Weights->vec);
//    Weights->vec[DataSet::getNumFeatures()] = 0;
    L2_SVM_MFN(train,Options,Weights,Outputs);
    for(int i= DataSet::getNumFeatures()+1;i--;)
      w[i]=Weights->vec[i];
  	delete [] Weights->vec;
  	delete Weights;
  	delete [] Outputs->vec;
  	delete Outputs;
    delete Options;
}

void Caller::trainEm(SetHandler & set,double * w, double Cpos, double Cneg, double fdr) {
  // iterate
  for(int i=0;i<niter;i++) {
    if(VERB>1) cerr << "Iteration " << i+1 << " : ";
    step(set,w,Cpos,Cneg,fdr);
    if(VERB>2) {cerr<<"Obtained weights" << endl;printWeights(cerr,w);}
  }
  if(VERB<=2) {cerr<<"Obtained weights" << endl;printWeights(cerr,w);}
}

void Caller::xvalidate(vector<DataSet *> &forward,vector<DataSet *> &shuffled, double *w) {
  Globals::getInstance()->decVerbose();
  vector<vector<DataSet *> > forwards(xval_fold),shuffleds(xval_fold);
  for(unsigned int j=0;j<xval_fold;j++) {
    forwards[j].resize(forward.size());
    for(unsigned int i=0;i<forward.size();i++) {
      forwards[j][i]=new VirtualSet(*(forward[i]),xval_fold,j);
    }  
  }
  for(unsigned int j=0;j<xval_fold;j++) {
    shuffleds[j].resize(shuffled.size());
    for(unsigned int i=0;i<shuffled.size();i++) {
      shuffleds[j][i]=new VirtualSet(*(shuffled[i]),xval_fold,j);
    }  
  }
  vector<SetHandler> train(xval_fold),test(xval_fold);
  for(unsigned int j=0;j<xval_fold;j++) {
    vector<DataSet *> ff(0),ss(0);
    for(unsigned int i=0;i<xval_fold;i++) {
      if (i==j)
        continue;
      ff.insert(ff.end(),forwards[i].begin(),forwards[i].end());
      ss.insert(ss.end(),shuffleds[i].begin(),shuffleds[i].end());
    }
    train[j].setSet(ff,ss);
    test[j].setSet(forwards[j],shuffleds[j]);
  }  
  vector<double> fdrs,cposs,cfracs;
  if (selectedfdr > 0) {
    fdrs.push_back(selectedfdr);
  } else {
    fdrs.push_back(0.01);fdrs.push_back(0.03); fdrs.push_back(0.07);
    if(VERB>0) cerr << "selecting fdr by cross validation" << endl;
  }
  if (selectedCpos > 0) {
    cposs.push_back(selectedCpos);
  } else {
    cposs.push_back(10);cposs.push_back(1);cposs.push_back(0.1);
    if(VERB>0) cerr << "selecting cpos by cross validation" << endl;
  }
  if (selectedCpos > 0 && selectedCneg > 0) {
    cfracs.push_back(selectedCneg/selectedCpos);
  } else  {
    cfracs.push_back(1);cfracs.push_back(3);cfracs.push_back(10);
    if(VERB>0) cerr << "selecting cneg by cross validation" << endl;  
  }
  int bestTP = 0;
  vector<double>::iterator fdr,cpos,cfrac;
  for(fdr=fdrs.begin();fdr!=fdrs.end();fdr++) {
    for(cpos=cposs.begin();cpos!=cposs.end();cpos++) {
      for(cfrac=cfracs.begin();cfrac!=cfracs.end();cfrac++) {
        if(VERB>0) cerr << "-cross validation with cpos=" << *cpos <<
          ", cfrac=" << *cfrac << ", fdr=" << *fdr << endl;
        int tp=0;
        double ww[DataSet::getNumFeatures()+1];
        for (unsigned int i=0;i<xval_fold;i++) {
          if(VERB>0) cerr << "cross calidation - fold " << i+1 << " out of " << xval_fold << endl;
          for(int ix=0;ix<DataSet::getNumFeatures()+1;ix++) ww[ix]=0;
          ww[3]=1;
          trainEm(train[i],ww,*cpos,(*cpos)*(*cfrac),*fdr);
          Scores sc;
          tp += sc.calcScores(ww,test[i],test_fdr);
          if(VERB>2) cerr << "Cumulative # of positives " << tp << endl;
        }
        if(VERB>0) cerr << "- cross validation found " << tp << " positives" << endl;
        if (tp>bestTP) {
          if(VERB>1) cerr << "Better than previous result, store this" << endl;
          bestTP = tp;
          selectedfdr=*fdr;
          selectedCpos = *cpos;
          selectedCneg = (*cpos)*(*cfrac);          
          for(int ix=0;ix<DataSet::getNumFeatures();ix++) w[ix]=ww[ix];
        }
      }     
    }
  }
  for(unsigned int j=0;j<xval_fold;j++) {
    for(unsigned int i=0;i<forward.size();i++) {
      delete forwards[j][i];
    }  
  }
  for(unsigned int j=0;j<xval_fold;j++) {
    for(unsigned int i=0;i<shuffled.size();i++) {
      delete shuffleds[j][i];
    }  
  }
  Globals::getInstance()->incVerbose();
}


int Caller::run() {
  time(&startTime);
  startClock=clock();
  if(VERB>0)  cerr << extendedGreeter();
  bool doShuffled2 = shuffled2FN.size()>0;
  vector<DataSet *> forward,shuffled,shuffled2;
  IntraSetRelation forRel,shuRel,shu2Rel;
  readFile(forwardFN,1,forward,&forRel);
  readFile(shuffledFN,-1,shuffled,&shuRel);
  if (doShuffled2)
    readFile(shuffled2FN,-1,shuffled2,&shu2Rel);
  
  SetHandler trainset,testset;
  trainset.setSet(forward,shuffled);
  testset.setSet(forward,(doShuffled2?shuffled2:shuffled));
  if (gistFN.length()>0) {
    trainset.gistWrite(gistFN);
  }
  //Normalize features
  vector<DataSet *> all;
  all.assign(forward.begin(),forward.end());
  all.insert(all.end(),shuffled.begin(),shuffled.end());  
  all.insert(all.end(),shuffled2.begin(),shuffled2.end());  
  pNorm=Normalizer::getNew();
  pNorm->setSet(all);
  pNorm->normalizeSet(all);
  
  // Set up a first guess of w
  double w[DataSet::getNumFeatures()+1];
  for(int ix=0;ix<DataSet::getNumFeatures()+1;ix++) w[ix]=0;
  w[3]=1;
//  w[DataSet::getNumFeatures()]=1;
  
  if (selectedfdr<=0 || selectedCpos<=0 || selectedCneg <= 0) {
    xvalidate(forward,shuffled,w);  
  }
  if(VERB>0) cerr << "---Training with Cpos=" << selectedCpos <<
          ", Cneg=" << selectedCneg << ", fdr=" << selectedfdr << endl;
    
  trainEm(trainset,w,selectedCpos,selectedCneg,selectedfdr);
  time_t end;
  time (&end);
  double diff = difftime (end,startTime);
  
  ostringstream timerValues;
  timerValues.precision(4);
  timerValues << "Processing took " << ((double)(clock()-startClock))/(double)CLOCKS_PER_SEC;
  timerValues << " cpu seconds or " << diff << " seconds wall time" << endl; 
  if (VERB>1) cerr << timerValues.str();
  Scores testScores;
  testScores.calcScores(w,testset,selectedfdr);
  if (modifiedFN.size()>0) {
    double ww[DataSet::getNumFeatures()+1];
    pNorm->unnormalizeweight(w,ww);    
    modifyFile(modifiedFN,forward,ww,testScores,extendedGreeter()+timerValues.str());
  }
  if (weightFN.size()>0) {
     ofstream weightStream(weightFN.data(),ios::out);
     printWeights(weightStream,w);
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






