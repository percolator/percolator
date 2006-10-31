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
bool DataSet::calcQuadraticFeatures = false;
bool DataSet::calcTrypticFeatures = true;
bool DataSet::chymoInsteadOfTryptic = false;
bool DataSet::calcIntraSetFeatures = true;

DataSet::DataSet()
{
   feature = NULL;
   n_examples=0;
   sqtFN = "";
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

void DataSet::setNumFeatures() {
  numRealFeatures= maxNumRealFeatures
                 - (calcTrypticFeatures?0:2) 
                 - (calcIntraSetFeatures?0:3);
  numFeatures=(calcQuadraticFeatures?numRealFeatures*(numRealFeatures+1)/2:numRealFeatures);
}

int DataSet::line2fields(string & line, vector<string> * words) {
  istringstream iss;
  iss.str(line);
  string word;
//  iss >> word;
  unsigned int n_w=0;
  while(iss.good()) {
    iss >> (*words)[n_w++];
  }
  unsigned int i = n_w;
  while (i<words->size())
    (*words)[i++]="";
  return n_w;
}
void DataSet::print_features() {
   for(int i=0;i<getSize();i++) {
	   for(int j=0;j<DataSet::getNumFeatures();j++) {
	      cout << j+1 << ":" << feature[DataSet::rowIx(i)+j] << " ";
	   }
	   cout << endl;
   }
}

void DataSet::print_10features() {
   cerr << getFeatureNames() << endl;
   for(int i=0;i<10;i++) {
       for(int j=0;j<DataSet::getNumFeatures();j++) {
          cerr << feature[DataSet::rowIx(i)+j] << "\t";
       }
       cerr << endl;
   }
}

double DataSet::isTryptic(const string & str) {
  assert(str[1]=='.');
  return (
  (((str[0]=='K' || str[0]=='R') && str[2]!= 'P') ||
  str[0]=='-' || str[2]=='-')
  ?1.0:0.0);
}

// [FHWY].[^P]
double DataSet::isChymoTryptic(const string & str) {
  assert(str[1]=='.');
  return (
  (((str[0]=='F' || str[0]=='H' || str[0]=='W' || str[0]=='Y') && str[2]!= 'P') ||
  str[0]=='-' || str[2]=='-')
  ?1.0:0.0);
}


double * DataSet::getNext(int& pos) {
  pos++;
  if (pos<0)
    pos=0;
  if(pos>=getSize())
    return NULL;
  return &feature[DataSet::rowIx(pos)];
}

const double * DataSet::getFeatures(const int pos) {
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

void DataSet::modifyRec(const string record,string &outStr,int ix, int mLines) {
  double feat[DataSet::numFeatures+1];
  if (mLines>3) mLines=3;
  istringstream in(record);
  ostringstream out;
  string line,tmp,lineRem;
  getLine(in,line);
  out << line;
  for(int m=0;m<mLines;m++) {
    readFeatures(record,feat,m,proteinIds[ix],pepSeq[ix],true);
    for(int a=0;a<4;a++) {
      in >> tmp;
      out << tmp << "\t";
    }
    // deltCn and XCorr and Sp
    in >> tmp >> tmp >> tmp;
    out << 0.0 << "\t" << sc[ix] << "\t" << -fdr[ix];
    getline(in,lineRem);
    out << lineRem << endl;
    // L lines
    while(in.peek()=='L' && getline(in,line)) {      
      assert(line[0]=='L');
      out << line << endl;
    }
  }
}

void DataSet::modify_sqt(const string outFN, vector<double> & sc, vector<double> & fdr,const string greet) {
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

  ostringstream buff;
  istringstream lineParse;
  int ix=0,lines=0,ms=0,charge=0;
  string tmp,lineRem;
  bool print = true;
  while (getline(sqtIn,line)) {
    if (line[0]=='S') {
      if(lines>1 && charge<=3) {
        string record=buff.str();
        readFeatures(record,&feature[rowIx(ix)],0,proteinIds[ix],pepSeq[ix],false);
        ix++;
        buff.str("");
        buff.clear();
      }
      lines=1;
      buff << line << endl;
      lineParse.str(line);
      lineParse >> tmp >> tmp >> tmp >> charge;
    }
    if (line[0]=='L'||(line[0]=='M' && ++ms)) {
      lines++;
      buff << line << endl;
    }
  }
  if(lines>1 && charge<=3) {
    string record=buff.str();
    readFeatures(record,&feature[rowIx(ix)],0,proteinIds[ix],pepSeq[ix],true);
  }

  while(getline(sqtIn,line)) {
    if(!print && line[0]!= 'M' && line[0] != 'L')
      print = true;
    while (line[0] == 'S') {
        char c;
        c=sqtIn.peek();
        if (c!='S')
          break;
        getline(sqtIn,line);
    }
    if (print)
      sqtOut << line << endl;
    if (line[0] == 'S') {
      string tmp,charge,scan;
      ix++;
      istringstream iss;
      iss.str(line);
      iss >> tmp >> tmp >> scan >> charge;
      string id = charge + '_' + scan;
      while (id!=ids[ix]){cout << "dropping " << ids[ix] <<endl;ix++;}
//      cout << id << " " << ids[ix] << endl;
//      assert(id ==ids[ix]);
      getline(sqtIn,line); // Get M line
      assert(line[0]=='M');
      iss.str(line);
      for(int a=0;a<4;a++) {
        iss >> tmp;
        sqtOut << tmp << "\t";
      }
      // deltCn and XCorr and Sp
      iss >> tmp >> tmp >> tmp;
      sqtOut << 0.0 << "\t" << sc[ix] << "\t" << -fdr[ix];
      getline(iss,lineRem);
      sqtOut << lineRem << endl;
      // L line
      getline(sqtIn,line);      
      assert(line[0]=='L');
      sqtOut << line << endl;
      sqtOut << "M\t2\t15\t600.0\t" << (sc[ix]+10)/(sc[ix]==0.0?1:sc[ix]) << "\t-10.0\t-1.0\t3\t10\tK.IAMAFAK.E\tU" << endl;
      sqtOut << "L\tBogusin" << endl;
      print=false;
    } 
  }
  sqtIn.close();
  sqtOut.close();
}
   
string DataSet::getFeatureNames() {
  ostringstream oss;
  oss << "RankSp\tdeltaMa\tdeltCn\tXcorr\tSp\tIonFrac\tMass\tPepLen\tCharge1\tCharge2\tCharge3";
  if (calcTrypticFeatures)
    oss << "\tenzN\tenzC";
  if (calcIntraSetFeatures)
    oss << "\tnumPep\tnumProt\tpepSite";
  return oss.str();
}

void DataSet::readFeatures(string &in,double *feat,int match,set<string> & proteins, string & pep, bool getIntra) {
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
      string sub1=pep.substr(0,3);
      string sub2=pep.substr(pep.size()-3);
      if (calcTrypticFeatures) {
        feat[11]=isEnz(sub1);        
        feat[12]=isEnz(sub2);
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
        double f2 = feature[ixf2];
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

void DataSet::read_sqt(const string fname, IntraSetRelation * intraRel) {
  intra=intraRel;
  setNumFeatures();
  sqtFN.assign(fname);
  int n = 0,charge=0;
  string line,tmp;
  istringstream lineParse;  
  ifstream sqtIn;
  sqtIn.open(sqtFN.data(),ios::in);
  if (!sqtIn) {
  	cerr << "Could not open file " << sqtFN << endl;
  	exit(-1);
  }
  while (getline(sqtIn,line)) {
    if (line[0]=='S' && sqtIn.peek() != 'S') {
         lineParse.str(line);  
         lineParse >> tmp >> tmp >> tmp >> charge;       
         if (charge <= 3) n++;
    }
  }
  if (VERB>1) cerr << n << " records in file " << sqtFN << endl;
  sqtIn.clear();
  sqtIn.seekg(0,ios::beg);
  
  proteinIds.resize(n);
  pepSeq.resize(n);
  string seq;

  feature = new double[n*DataSet::getNumFeatures()];
  ostringstream buff;
  ids.resize(n,"");
  n_examples=n;
  int ix=0,lines=0;
  string id,scan,pep;
  while (getline(sqtIn,line)) {
    if (line[0]=='S') {
      if(lines>1 && charge<=3) {
        string record=buff.str();
        readFeatures(record,&feature[rowIx(ix)],0,proteinIds[ix],pepSeq[ix],false);
        intra->registerRel(pepSeq[ix],proteinIds[ix]);
        ix++;
        buff.str("");
        buff.clear();
      }
      lines=1;
      buff << line << endl;
      lineParse.str(line);
      lineParse >> tmp >> tmp >> tmp >> charge;
    }
    if (line[0]=='M' || line[0]=='L') {
      lines++;
      buff << line << endl;
    }
  }
  if(lines>1 && charge<=3) {
    string record=buff.str();
    readFeatures(record,&feature[rowIx(ix)],0,proteinIds[ix],pepSeq[ix],false);
    intra->registerRel(pepSeq[ix],proteinIds[ix]);
  }
  sqtIn.close();
//  cout << "Read File" << endl;
}


