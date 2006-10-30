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

void DataSet::modify_sqt(const string outFN, vector<double> & sc, vector<double> & fdr,const string greet) {
  int ix=-1;
  string line,lineRem;
  bool print = true;
  ifstream sqtIn(sqtFN.data(),ios::in);
  ofstream sqtOut(outFN.data(),ios::out);
  istringstream greetStream(greet);
  sqtOut.precision(5);
  sqtOut << "H\tfile massaged by" << endl;
  while(getline(greetStream,line)) {
    sqtOut << "H\t" << line << endl;
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

void DataSet::computeIntraSetFeatures() {
  if (DataSet::calcIntraSetFeatures) {
    for(int ix=0;ix<n_examples;ix++) {
      feature[rowIx(ix)+numRealFeatures-3]=
        log((double)intra->getNumPep(pepSeq[ix]));
      feature[rowIx(ix)+numRealFeatures-2]=
        log((double)intra->getNumProt(proteinIds[ix]));
      feature[rowIx(ix)+numRealFeatures-1]=
        log((double)intra->getPepSites(proteinIds[ix]));
//      cerr << ix << " " << feature[rowIx(ix)+numRealFeatures-3] << " " << feature[rowIx(ix)+numRealFeatures-2] << " "  << feature[rowIx(ix)+numRealFeatures-1] << endl;
    }
  }
  if (DataSet::calcQuadraticFeatures) {
    for (int r=0;r<getSize();r++){
      int ix = numRealFeatures;
      for (int ixf1=1;ixf1<numRealFeatures;ixf1++){
        double f1 = feature[rowIx(r)+ixf1];
        for (int ixf2=0;ixf2<ixf1;ixf2++,ix++){
          double f2 = feature[rowIx(r)+ixf2];
          double fp=f1*f2;
          double newFeature;
          if (fp>=0.0) {
            newFeature=sqrt(fp);
          } else {
            newFeature=-sqrt(-fp);
          }
          feature[rowIx(r)+ix]=newFeature;
        }        
      }
      assert(ix==numFeatures);    
    }
  }
  
  return;
}

void DataSet::read_sqt(const string fname, IntraSetRelation * intraRel) {
  intra=intraRel;
  setNumFeatures();
  sqtFN.assign(fname);
  int n = 0;
  vector<string> fields;
  fields.resize(25,"");

  string line;
  ifstream sqtIn;
  sqtIn.open(sqtFN.data(),ios::in);
  if (!sqtIn) {
  	cerr << "Could not open file " << sqtFN << endl;
  	exit(-1);
  }
  while (getline(sqtIn,line)) {
    if (line[0]=='S') {
         getline(sqtIn,line);  // Protect our selves against double S-line errors
         n++;
    }
  }
  if (VERB>1) cerr << n << " records in file " << sqtFN << endl;
  sqtIn.clear();
  sqtIn.seekg(0,ios::beg);
  
  set<string> proteins;
  proteinIds.resize(n);
  pepSeq.resize(n);
  string seq;

  feature = new double[n*DataSet::getNumFeatures()];
  ids.resize(n,"");
  charge.resize(n,0);
  n_examples=n;
  int ix=-1,chrg;
  bool gotL = true,gotDeltCn=true;
  string id;
  double mass;
  while (getline(sqtIn,line)) {
    if (line[0]=='S') {
      line2fields(line,&fields);
      id = fields[3];
      id += '_';
      id += fields[2];
      chrg=atoi(fields[3].data());
      mass=atof(fields[6].data());
      gotL = 0;
    }
    if (line[0]=='M' && !gotDeltCn) {
      line2fields(line,&fields);
      feature[DataSet::rowIx(ix)+2]=atof(fields[4].data());
      gotDeltCn = true;
    }
    if (line[0]=='M' && !gotL) {
      ids[++ix]=id;
      charge[ix]=chrg;
      proteins.clear();
      line2fields(line,&fields);
      feature[rowIx(ix)+0]=atof(fields[2].data());      // rank by Sp
      feature[rowIx(ix)+1]=mass-atof(fields[3].data()); // obs - calc mass
      feature[rowIx(ix)+2]=0.0;                         // deltCn (leave until next M line)
      feature[rowIx(ix)+3]=atof(fields[5].data());      // Xcorr
      feature[rowIx(ix)+4]=atof(fields[6].data());      // Sp
      feature[rowIx(ix)+5]=atof(fields[7].data())/atof(fields[8].data()); //Fraction matched/expected ions
      feature[rowIx(ix)+6]=mass;                        // Observed mass
      seq=fields[9];
      pepSeq[ix]=fields[9];
      feature[rowIx(ix)+7]=seq.size()-4;                // Peptide length
      feature[rowIx(ix)+8]=(charge[ix]==1?1.0:0.0);     // Charge
      feature[rowIx(ix)+9]=(charge[ix]==2?1.0:0.0);
      feature[rowIx(ix)+10]=(charge[ix]==3?1.0:0.0);
      string sub1=seq.substr(0,3);
      string sub2=seq.substr(seq.size()-3);
      if (calcTrypticFeatures) {
        feature[rowIx(ix)+11]=isEnz(sub1);        
        feature[rowIx(ix)+12]=isEnz(sub2);
      }
      gotDeltCn = false;
    }
    if (line[0]=='L' && !gotL) {
      if (sqtIn.peek() != 'L') gotL=true;
      line2fields(line,&fields);
      proteins.insert(fields[1]);
      if (gotL) {
        proteinIds[ix].insert(proteins.begin(),proteins.end());
        intra->registerRel(seq,proteins);
      }
    }
  }
  sqtIn.close();
//  cout << "Read File" << endl;
}


