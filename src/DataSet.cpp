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
#include "Scores.h"

int DataSet::numFeatures = numRealFeatures;
bool DataSet::calcQuadraticFeatures = false;

DataSet::DataSet()
{
   delim.assign(" \t\n");
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

int DataSet::line2fields(string & line, vector<string> * words) {
  words->clear();
  istringstream iss;
  iss.str(line);
  string word;
//  iss >> word;
  while(iss.good()) {
    iss >> word;
    words->push_back(word);
//    iss >> word;
  }
  return words->size();
}
void DataSet::print_features() {
   for(int i=0;i<getSize();i++) {
	   for(int j=0;j<DataSet::getNumFeatures();j++) {
	      cout << j+1 << ":" << feature[DataSet::rowIx(i)+j] << " ";
	   }
	   cout << endl;
   }
}

double DataSet::isTryptic(const string & str) {
  assert(str[1]=='.');
  return (((str[0]=='K' || str[0]=='R') && str[2]!= 'P')?1.0:0.0);
}


double * DataSet::getNext(int& pos) {
  pos++;
  if (pos<0)
    pos=0;
  if(pos>=getSize())
    return NULL;
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

void DataSet::modify_sqt(string & outFN, vector<double> & sc, vector<double> & fdr) {
  int ix=-1;
  string line,lineRem;
  bool print = true;
  ifstream sqtIn(sqtFN.data(),ios::in);
  ofstream sqtOut(outFN.data(),ios::out);
  sqtOut.precision(5);
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
      print=false;
    } 
  }
  sqtIn.close();
  sqtOut.close();
}

void DataSet::read_sqt(string & fname) {
  sqtFN.assign(fname);
  int n = 0;
  vector<string> fields;

  string line;
  ifstream sqtIn;
  sqtIn.open(sqtFN.data(),ios::in);
  if (sqtIn.fail()) {
  	cerr << "Could not open file " << sqtFN << endl;
  	exit(-1);
  }
  while (getline(sqtIn,line)) {
    if (line[0]=='S') {
         getline(sqtIn,line);  // Protect ourselfs against double S-line errors
         n++;
    }
  }
  cout << n << " records in file " << sqtFN << endl;
  sqtIn.clear();
  sqtIn.seekg(0,ios::beg);

  map<string, vector<int> > protids2ix;
  feature = new double[n*DataSet::getNumFeatures()];
  ids.resize(n,"");
  charge.resize(n,0);
  vector<string> ix2seq(n,"");
  n_examples=n;
  int ix=-1,gotL = 1,gotDeltCn=1,chrg;
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
      gotDeltCn = 1;
    }
    if (line[0]=='M' && !gotL) {
      ids[++ix]=id;
      charge[ix]=chrg;
      line2fields(line,&fields);
      feature[DataSet::rowIx(ix)+0]=atof(fields[2].data());
      feature[DataSet::rowIx(ix)+1]=mass - atof(fields[3].data());
//      feature[DataSet::rowIx(ix)+2]=atof(fields[4].data());
      feature[DataSet::rowIx(ix)+2]=0.0;
      feature[DataSet::rowIx(ix)+3]=atof(fields[5].data());
      feature[DataSet::rowIx(ix)+4]=atof(fields[6].data());
      feature[DataSet::rowIx(ix)+5]=atof(fields[7].data())/atof(fields[8].data());
      string seq(fields[9]);
      string sub1=seq.substr(0,3);
      string sub2=seq.substr(seq.size()-3);
      ix2seq[ix].assign(seq);
      feature[DataSet::rowIx(ix)+6]=isTryptic(sub1);
      feature[DataSet::rowIx(ix)+7]=isTryptic(sub2);
      feature[DataSet::rowIx(ix)+8]=(charge[ix]==1?1.0:0.0);
      feature[DataSet::rowIx(ix)+9]=(charge[ix]==2?1.0:0.0);
      feature[DataSet::rowIx(ix)+10]=(charge[ix]==3?1.0:0.0);
      gotDeltCn = 0;
    }
    if (line[0]=='L' && !gotL) {
      gotL=1;
      line2fields(line,&fields);
//      string * prot_id= new string(fields[1]);
//      protids2ix[*prot_id].push_back(ix);
      protids2ix[fields[1]].push_back(ix);
    }
  }
  sqtIn.close();
  map<string,int> seqfreq;
  map<string, vector<int> >::iterator ixvec;
  for( ixvec = protids2ix.begin(); ixvec != protids2ix.end(); ixvec++ ) {
  	 seqfreq.clear();
     double f1=(double)ixvec->second.size();
     for (unsigned int i=0;i<ixvec->second.size();i++) {
       seqfreq[ix2seq[ixvec->second[i]]]=(seqfreq.count(ix2seq[ixvec->second[i]])>0?seqfreq[ix2seq[ixvec->second[i]]]+1:1);
     }
     for (unsigned int i=0;i<ixvec->second.size();i++) {
       feature[DataSet::rowIx(ixvec->second[i])+11]=log((float)seqfreq.size());
       feature[DataSet::rowIx(ixvec->second[i])+12]=log((float)seqfreq[ix2seq[ixvec->second[i]]]);
       feature[DataSet::rowIx(ixvec->second[i])+13]=f1;
     }
     
  }            
//  cout << "Read File" << endl;
  if (DataSet::calcQuadraticFeatures) {
    for (int r=0;r<getSize();r++){
      int ix = DataSet::numRealFeatures;
      for (int ixf1=1;ixf1<DataSet::numRealFeatures;ixf1++){
        double f1 = feature[DataSet::rowIx(r)+ixf1];
        for (int ixf2=0;ixf2<ixf1;ixf2++,ix++){
          double f2 = feature[DataSet::rowIx(r)+ixf2];
          double fp=f1*f2;
          double newFeature;
          if (fp>=0.0) {
            newFeature=sqrt(fp);
          } else {
            newFeature=-sqrt(-fp);
          }
          feature[DataSet::rowIx(r)+ix]=newFeature;
        }        
      }
      assert(ix==numFeatures);    
    }
  }
}


