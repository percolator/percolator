#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <set>
#include <map>
#include <utility>
#include <vector>
#include <string>
using namespace std;
#include "DataSet.h"
#include "IntraSetRelation.h"
#include "Scores.h"
#include "Globals.h"

int DataSet::numFeatures = maxNumRealFeatures;
int DataSet::numRealFeatures = maxNumRealFeatures;
int DataSet::hitsPerSpectrum = 1;
bool DataSet::calcQuadraticFeatures = false;
bool DataSet::calcAAFrequencies = true;
bool DataSet::calcTrypticFeatures = true;
bool DataSet::chymoInsteadOfTryptic = false;
bool DataSet::calcIntraSetFeatures = true;
string DataSet::featureNames = "";
string DataSet::aaAlphabet = "ACDEFGHIKLMNPQRSTVWY";

DataSet::DataSet()
{
   feature = NULL;
   n_examples=0;
   sqtFN = "";
   pattern="";
   doPattern=false;
   matchPattern=true;
}

DataSet::~DataSet()
{
	if (feature) {
//		for(int i=0;i<n_feature;i++)
//	   		delete feature[i];
		delete [] feature;
		feature=NULL;
	}
}

double * DataSet::getNext(int& pos) {
  pos++;
  if (pos<0)
    pos=0;
  if(pos>=getSize())
    return NULL;
  return &feature[DataSet::rowIx(pos)];
}

const double * DataSet::getFeatures(const int pos) const {
  return &feature[DataSet::rowIx(pos)];
}

bool DataSet::getGistDataRow(int & pos,string &out){
  ostringstream s1;
  double * feature = NULL;
//  while (!feature || charge[pos] !=2)  { //tmp fix
  if ((feature = getNext(pos)) == NULL) return false; 
//  }
  s1 << ids[pos];
  for (int ix = 0;ix<DataSet::getNumFeatures();ix++) {
    s1 << '\t' << feature[ix];
  }
  s1 << endl;
  out = s1.str();
  return true;
}

void DataSet::print_features() {
   for(int i=0;i<getSize();i++) {
       for(int j=0;j<DataSet::getNumFeatures();j++) {
          cerr << j+1 << ":" << feature[DataSet::rowIx(i)+j] << " ";
       }
       cerr << endl;
   }
}

void DataSet::print_10features() {
   cerr << DataSet::getFeatureNames() << endl;
   for(int i=0;i<10;i++) {
       for(int j=0;j<DataSet::getNumFeatures();j++) {
          cerr << feature[DataSet::rowIx(i)+j] << "\t";
       }
       cerr << endl;
   }
}

void DataSet::print(Scores& test, vector<pair<double,string> > &outList) {
  ostringstream out;
  int ix =-1;
  while (double * features=getNext(ix)) {
    double score = test.calcScore(features);
    double q = test.getQ(score);
    out << ids[ix] << "\t" << score << "\t" << q;
    if ((int)pepSeq.size()>ix)
        out << "\t" << pepSeq[ix];
    if ((int)proteinIds.size()>ix) {
      set<string> prots = proteinIds[ix];
      set<string>::const_iterator it = prots.begin();
      for(;it!=prots.end();it++) {
        out << "\t" << *it;
      }
    }
    pair<double,string> op(score,out.str());    
    outList.push_back(op);
    out.str("");
  }
}

void DataSet::setNumFeatures() {
  numRealFeatures= maxNumRealFeatures
                 - (calcTrypticFeatures?0:2) 
                 - (calcIntraSetFeatures?0:3)
                 - (calcAAFrequencies?0:aaAlphabet.size());
  numFeatures=(calcQuadraticFeatures?numRealFeatures*(numRealFeatures+1)/2:numRealFeatures);
}


double DataSet::isTryptic(const char n,const char c) {
  return (
  (((n=='K' || n=='R') && c != 'P') ||
  n=='-' || c=='-')
  ?1.0:0.0);
}

// [FHWY].[^P]
double DataSet::isChymoTryptic(const char n,const char c) {
  return (
  (((n=='F' || n=='H' || n=='W' || n=='Y') && c!= 'P') ||
  n=='-' || c=='-')
  ?1.0:0.0);
}

unsigned int DataSet::cntEnz(const string& peptide) {
    unsigned int pos=2, cnt=0;
    char n = peptide.at(pos++);
    while (pos<peptide.size()-2) {
      char c = peptide.at(pos++);
      if (isEnz(n,c))
        cnt++;
      n=c;
    }
    return cnt;
}

void DataSet::readGistData(ifstream & is, vector<unsigned int> ixs) {
  string tmp,line;
  is.clear();
  is.seekg(0,ios::beg);
  getline(is,line);
  getline(is,line);
  unsigned int m=0,n=ixs.size();
  istringstream buff(line);
  while(true) {
    buff >> tmp;
    if (!buff) break;
    m++;
  }
  if (m<3) {
    cerr << "To few features in Gist data file";
    exit(-1);
  }
  m--; // remove id line
  DataSet::numFeatures = m;

  proteinIds.resize(0);
  pepSeq.resize(0);
  string seq;

  feature = new double[n*DataSet::getNumFeatures()];
  ids.resize(n,"");
  n_examples=n;

  is.clear();
  is.seekg(0,ios::beg);
  getline(is,line); // id line

  unsigned int ix = 0;
  getline(is,line);
  for(unsigned int i=0;i<n;i++) {
    while (ix<ixs[i]){
    	getline(is,line);
    	ix++;
    }
    buff.str(line);
    buff.clear();
    buff >> ids[i];
    double *featureRow=&feature[rowIx(i)];
    for(register unsigned int j=0;j<m;j++) {
      buff >> featureRow[j];
    } 
  } 
}

string DataSet::modifyRec(const string record,const set<int>& theMs, const double *w, Scores * pSc, bool dtaSelect) {
  double feat[DataSet::numFeatures];
//  if (mLines>3) mLines=3;
  vector<pair<double,string> > outputs;
  istringstream in(record);
  ostringstream out,outtmp;
  string line,tmp,lineRem;
  getline(in,line);
  out << line << endl;
  set<string> proteinsTmp;
  string tmpPepSeq;
  double rSp,mass;
  set<int>::const_iterator it;
  for(it=theMs.begin();it!=theMs.end();it++) {
    proteinsTmp.clear();
    readFeatures(record,feat,*it,proteinsTmp,tmpPepSeq,true);
    double score = 0;
    for (int i=DataSet::numFeatures;i--;)
      score += feat[i]*w[i];
    score += w[DataSet::numFeatures];
    double q = pSc->getQ(score);
    in >> tmp >> tmp >> rSp >> mass;
//    outtmp << "M\t%i\t" << rSp << "\t" << mass << "\t";
    outtmp << "M\t" << tmp << "\t" << rSp << "\t" << mass << "\t";
    // deltCn and XCorr and Sp
    in >> tmp >> tmp >> tmp;
    outtmp << "%6.4g\t" << score << "\t" << -q;
    getline(in,lineRem);
    outtmp << lineRem << endl;
    // L lines
    while(in.peek()=='L' && getline(in,line)) {      
      assert(line[0]=='L');
      outtmp << line << endl;
    }
    outputs.push_back(pair<double,string>(score,outtmp.str()));
    outtmp.clear();
    outtmp.str("");
  }
  if(dtaSelect) {
    outputs.push_back(pair<double,string>(-10,
      "M\t600\t600\t1\t%6.4g\t-10\t-1\t0\t0\tI.AMINVALI.D\tU\nL\tPlaceholder satisfying DTASelect\n"));
  }
  sort(outputs.begin(),outputs.end());
  reverse(outputs.begin(),outputs.end());
  double x0=0,delt;
  char buf[1024];
  for(unsigned int i=0;i<outputs.size();i++) {
    string aStr = outputs[i].second;
    double score = outputs[i].first;
    if (i==0) x0=score;
    if (x0 != 0) {
      delt = (x0-score)/x0;
    } else delt = 0;
//    sprintf(buf,outputs[i].second.c_str(),i+1,delt);
    sprintf(buf,aStr.c_str(),delt);
    out << buf;
  }
  return out.str();
}

void DataSet::modify_sqt(const string & outFN, const double *w, Scores * pSc ,const string greet, bool dtaSelect) {
  string line;
  ifstream sqtIn(sqtFN.data(),ios::in);
  ofstream sqtOut(outFN.data(),ios::out);
  istringstream greetStream(greet);
  sqtOut.precision(5);
  sqtOut << "H\tfile massaged by" << endl;
  while(getline(greetStream,line)) {
    sqtOut << "H\t" << line << endl;
  }
  sqtOut << "H\t" << "InputFile: " << sqtFN << endl;
  sqtOut << "H\t" << "OutputFile: " << outFN << endl;
  if (VERB>1) cerr << "Writing Output to sqt file " << outFN << endl;
  
  ostringstream buff;
  istringstream lineParse;
  int lines=0,ms=0,charge=0;
  string tmp,lineRem;
  set<int> theMs;
  while (getline(sqtIn,line)) {
    if (line[0]=='H')
      sqtOut << line <<endl;
    if (line[0]=='S') {
      if(lines>1 && charge<=3) {
        if (ms>hitsPerSpectrum)
          ms=hitsPerSpectrum;
        string record=buff.str();
        sqtOut << modifyRec(record,theMs, w, pSc, dtaSelect);
      }
      buff.str("");
      buff.clear();
      lines=1;
      buff << line << endl;
      lineParse.str(line);
      lineParse >> tmp >> tmp >> tmp >> charge;
      ms=0;
      theMs.clear();
    }
    if (line[0]=='M') {
      ++ms;
      ++lines;
      buff << line << endl;    
    }
    if (line[0]=='L') {
      ++lines;
      buff << line << endl;
      if((int)theMs.size()<hitsPerSpectrum && (!doPattern || (line.find(pattern,0)!= string::npos)==matchPattern)) {
        theMs.insert(ms-1);
      }
    }
  }
  if(lines>1 && charge<=3) {
    string record=buff.str();
    sqtOut << modifyRec(record,theMs, w, pSc, dtaSelect);
  }
  sqtIn.close();
  sqtOut.close();
}
   
string DataSet::getFeatureNames() {
  if (featureNames.empty()) {
    ostringstream oss;
    oss << "RankSp\tdeltaMa\tdeltCn\tXcorr\tSp\tIonFrac\tMass\tPepLen\tCharge1\tCharge2\tCharge3";
    if (calcTrypticFeatures)
      oss << "\tenzN\tenzC";
    if (calcAAFrequencies) {
      for (string::const_iterator it=aaAlphabet.begin();it!=aaAlphabet.end();it++)
        oss << "\t" << *it << "Freq";
    }
    if (calcIntraSetFeatures)
      oss << "\tnumPep\tnumProt\tpepSite";
    featureNames = oss.str();
  }
  return featureNames;
}

void DataSet::readFeatures(const string &in,double *feat,int match,set<string> & proteins, string & pep, bool getIntra) {
  istringstream instr(in),linestr;
  string line,tmp;
  int charge;
  double mass,deltCn,otherXcorr=0,xcorr=0;
  bool gotL=true,gotDeltCn=(match==0);
  int ms=0;
  
  while (getline(instr,line)) {
    if (line[0]=='S') {
      linestr.clear();
      linestr.str(line);
      linestr >> tmp >> tmp >> tmp >> charge >> tmp >> tmp >> mass;
    }
    if (line[0]=='M' && !gotDeltCn) {
      linestr.clear();
      linestr.str(line);
      linestr >> tmp >> tmp >> tmp >> tmp >> deltCn >> otherXcorr;
      gotDeltCn = true;
    }
    if ((line[0]=='M') && (match==ms++)) {
      double rSp,cMass,sp,matched,expected;
      linestr.clear();
      linestr.str(line);
      linestr >> tmp >> tmp >> rSp >> cMass >> tmp >> xcorr >> sp >> matched >> expected >> pep;
      
      feat[0]=rSp;      // rank by Sp
      feat[1]=mass-cMass; // obs - calc mass
      feat[2]=0.0;                         // deltCn (leave until next M line)
      feat[3]=xcorr;      // Xcorr
      feat[4]=sp;      // Sp
      feat[5]=matched/expected; //Fraction matched/expected ions
      feat[6]=mass;                        // Observed mass
      feat[7]=pep.size()-4;                // Peptide length
      feat[8]=(charge==1?1.0:0.0);     // Charge
      feat[9]=(charge==2?1.0:0.0);
      feat[10]=(charge==3?1.0:0.0);
      int nxtFeat=11;
      if (calcTrypticFeatures) {
        feat[nxtFeat++]=isEnz(pep.at(0),pep.at(2));        
        feat[nxtFeat++]=isEnz(pep.at(pep.size()-3),pep.at(pep.size()-1));
        feat[nxtFeat++]=(double)cntEnz(pep);
      }
      if (calcAAFrequencies) {
        computeAAFrequencies(pep,&feat[nxtFeat]);
        nxtFeat += aaAlphabet.size();
      }
      gotDeltCn = (match!=0);
      gotL = false;
    }
    if (line[0]=='L' && !gotL) {
      if (instr.peek() != 'L') gotL=true;
      string p;
      linestr.clear();
      linestr.str(line);
      linestr >> tmp >> p;      
      proteins.insert(p.c_str());
    }
  }
//  feat[2]=deltCn;
  if (xcorr>0)
    feat[2]=(xcorr-otherXcorr)/xcorr;
  if (!isfinite(feat[2])) cerr << in;
  if (getIntra)
    computeIntraSetFeatures(feat,pep,proteins);
}

void DataSet::computeAAFrequencies(const string& pep, double *feat) {
  string::size_type pos = aaAlphabet.size();
  for (;pos--;) {feat[pos]=0.0;}
  int len=0;
  for (string::const_iterator it=pep.begin();it!=pep.end();it++) {
    pos=aaAlphabet.find(*it);
    if (pos!=string::npos) feat[pos]++;
    len++;
  }
  for (pos = aaAlphabet.size();pos--;) {feat[pos]/=len;}

}

void DataSet::computeIntraSetFeatures(double * feat,string &pep,set<string> &prots) {
  if (DataSet::calcIntraSetFeatures) {
    feat[numRealFeatures-3]=
      log((double)intra->getNumPep(pep));
    feat[numRealFeatures-2]=
      log((double)intra->getNumProt(prots));
    feat[numRealFeatures-1]=
      log((double)intra->getPepSites(prots));
  }
  if (DataSet::calcQuadraticFeatures) {
    int ix = numRealFeatures;
    for (int ixf1=1;ixf1<numRealFeatures;ixf1++){
      double f1 = feat[ixf1];
      for (int ixf2=0;ixf2<ixf1;ixf2++,ix++){
        double f2 = feat[ixf2];
        double fp=f1*f2;
        double newFeature;
        if (fp>=0.0) {
          newFeature=sqrt(fp);
        } else {
          newFeature=-sqrt(-fp);
        }
        feat[ix]=newFeature;
      }        
    }
    assert(ix==numFeatures);    
  }
}

void DataSet::computeIntraSetFeatures() {
  for(int row=0;row<n_examples;row++) {
    computeIntraSetFeatures(&feature[rowIx(row)],pepSeq[row],proteinIds[row]);
  }
  return;
}

void DataSet::read_sqt(const string fname, IntraSetRelation * intraRel,const string & wild, bool match) {
  intra=intraRel;
  matchPattern=match;
  pattern=wild;
  doPattern = !wild.empty();
  setNumFeatures();
  sqtFN.assign(fname);
  int n = 0,charge=0,ms=0;
  string line,tmp,prot;
  istringstream lineParse;  
  ifstream sqtIn;
  sqtIn.open(sqtFN.data(),ios::in);
  if (!sqtIn) {
  	cerr << "Could not open file " << sqtFN << endl;
  	exit(-1);
  }
  bool look = false;
  while (getline(sqtIn,line)) {
    if (line[0]=='S' && sqtIn.peek() != 'S') {
         lineParse.str(line);  
         lineParse >> tmp >> tmp >> tmp >> charge;       
         if (charge <= 3) look=true;
         ms=0;
    }
    if (look & line[0]=='L' && ms < hitsPerSpectrum) {
         lineParse.str(line);  
         lineParse >> tmp >> prot;
         if(!doPattern || ((line.find(wild,0)!= string::npos)==match)) {
             ++ms;
             ++n;
         }
    }
  }
  if (VERB>1) cerr << n << " records in file " << sqtFN << endl;
  sqtIn.clear();
  sqtIn.seekg(0,ios::beg);
  
  ids.resize(n,"");
  proteinIds.resize(n);
  pepSeq.resize(n);
  string seq;
  
  string fileId(fname);
  unsigned int spos = fileId.rfind('/');
  if (spos!=string::npos)
    fileId.erase(0,spos+1);
  spos = fileId.find('.');
  if (spos!=string::npos)
    fileId.erase(spos);
  
  feature = new double[n*DataSet::getNumFeatures()];
  ostringstream buff,id;
  n_examples=n;
  
  int ix=0,lines=0;
  string scan;
  set<int> theMs;
  while (getline(sqtIn,line)) {
    if (line[0]=='S') {
      if(lines>1 && charge<=3) {
        string record=buff.str();
        string idstr = id.str();
        set<int>::const_iterator it;
        for(it=theMs.begin();it!=theMs.end();it++) {
          id.str("");id << idstr << '_' << (*it +1);
          ids[ix]=id.str();
          readFeatures(record,&feature[rowIx(ix)],*it,proteinIds[ix],pepSeq[ix],false);
          intra->registerRel(pepSeq[ix],proteinIds[ix]);
          ix++;
        }
      }
      buff.str("");
      buff.clear();
      id.str("");
      lines=1;
      buff << line << endl;
      lineParse.clear();
      lineParse.str(line);
      lineParse >> tmp >> tmp >> scan >> charge;
      id << fileId << '_' << scan << '_' << charge;
      ms=0;
      theMs.clear();
    }
    if (line[0]=='M') {
      ++ms;
      ++lines;
      buff << line << endl;    
    }
    if (line[0]=='L') {
      ++lines;
      buff << line << endl;
      if((int)theMs.size()<hitsPerSpectrum && (!doPattern || (line.find(wild,0)!= string::npos)==match)) {
        theMs.insert(ms-1);
      }
    }
  }
  if(lines>1 && charge<=3) {
    string record=buff.str();
    string idstr = id.str();
    set<int>::const_iterator it;
    for(it=theMs.begin();it!=theMs.end();it++) {
      id.str("");id << idstr << '_' << *it;
      ids[ix]=id.str();
      readFeatures(record,&feature[rowIx(ix)],*it,proteinIds[ix],pepSeq[ix],false);
      intra->registerRel(pepSeq[ix],proteinIds[ix]);
      ix++;
    }
  }
  sqtIn.close();
//  cout << "Read File" << endl;
}


