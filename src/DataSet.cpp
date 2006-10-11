#include <iostream>
#include <stdio.h>
#include <math.h>
#include <set>
#include <map>
#include <utility>
#include <vector>
#include <string>
using namespace std;
#include "DataSet.h"
//using namespace std;


DataSet::DataSet()
{
   delim.assign(" \t\n");
   feature = NULL;
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

int DataSet::line2fields(char * str, vector<string> *fields) {
  char *brkt,*wrd;
  int ix = 0;
  wrd=strtok_r(str, delim.data(), &brkt);
  while(wrd) {
//    fields[ix] = (char*)malloc(sizeof(char)*100);//tmp.copy();
    (*fields)[ix++].assign(wrd);
    wrd=strtok_r(NULL, delim.data(), &brkt);
  }
  return ix;
}
void DataSet::print_features() {
   for(int i=0;i<n_feature;i++) {
	   for(int j=0;j<NUM_FEATURE;j++) {
	      cout << j+1 << ":" << feature[ROW(i)+j] << " ";
	   }
	   cout << endl;
   }
}

bool DataSet::isTryptic(const string & str) {
  assert(str[1]=='.');
  return ((str[0]=='K' || str[0]=='R') && str[2]!= 'P');
}

int DataSet::getIsoChargeSize(int c){
  int n=0;
  for (int ix=0;ix<(signed int)charge.size();ix++) {
    if (charge[ix]==c) n++;
  }
  return n;
}

double * DataSet::getNext(const int c, int& pos) {
  pos++;
  if (pos<0)
    pos=0;
  while((pos<n_feature) && (charge[pos]!=c))
    pos++;
  if(pos>=n_feature)
    return NULL;
  return &feature[ROW(pos)];
}

void DataSet::read_sqt(string & fname) {
  FILE *fp1;
  int n = 0;
  size_t len1 = 0;
  vector<string> fields(30,"");

  char * str = (char *) calloc(1023,sizeof(char));

  if ((fp1=fopen(fname.data(),"r"))==NULL) {
  	cerr << "Could not open file " << fname << endl;
  	exit(-1);
  }
  while (getline(&str, &len1, fp1) != -1) {
    if (str[0]=='S') {
/*      line2fields(str,&fields);
      int charge = atoi(fields[3].data());
      if (charge == DataSet::ONLY_CHARGE) */
         n++;
    }
  }
  cout << n << " records in file" << endl;

  rewind(fp1);

  map<string, vector<int> > protids2ix;
  feature = new double[n*NUM_FEATURE];
  ids.resize(n,"");
  charge.resize(n,0);
  vector<string> ix2seq(n,"");
//  for(int i=0;i<n;i++) {
//    feature[i]=new double[NUM_FEATURE];
//  }
  n_feature=n;
  int ix=-1,gotL = 1,gotDeltCn=1;
  double mass;
  while (getline(&str, &len1, fp1) != -1) {
    if (str[0]=='S') {
      line2fields(str,&fields);
      ids[++ix]+=fields[3];
      ids[ix]+='_';
      ids[ix]+=fields[2];
      mass=atof(fields[6].data());
      gotL = 0;
      charge[ix]=atoi(fields[3].data());
    }
    if (str[0]=='M' && !gotDeltCn) {
      line2fields(str,&fields);
      feature[ROW(ix)+2]=atof(fields[4].data());
      gotDeltCn = 1;
    }
    if (str[0]=='M' && !gotL) {
      line2fields(str,&fields);
      feature[ROW(ix)+0]=atof(fields[2].data());
      feature[ROW(ix)+1]=mass - atof(fields[3].data());
//      feature[ROW(ix)+2]=atof(fields[4].data());
      feature[ROW(ix)+2]=0.0;
      feature[ROW(ix)+3]=atof(fields[5].data());
      feature[ROW(ix)+4]=atof(fields[6].data());
      feature[ROW(ix)+5]=atof(fields[7].data())/atof(fields[8].data());
      string seq(fields[9]);
      string sub1=seq.substr(0,3);
      string sub2=seq.substr(seq.size()-3);
      ix2seq[ix].assign(seq);
      feature[ROW(ix)+6]=isTryptic(sub1);
      feature[ROW(ix)+7]=isTryptic(sub2);
      gotDeltCn = 0;
    }
    if (str[0]=='L' && !gotL) {
      gotL=1;
      line2fields(str,&fields);
      string * prot_id= new string(fields[1]);
      protids2ix[*prot_id].push_back(ix);
    }
  }
  fclose(fp1);
  map<string,int> seqfreq;
  map<string, vector<int> >::iterator ixvec;
  for( ixvec = protids2ix.begin(); ixvec != protids2ix.end(); ixvec++ ) {
  	 seqfreq.clear();
     double f1=(double)ixvec->second.size();
     for (unsigned int i=0;i<ixvec->second.size();i++) {
       seqfreq[ix2seq[ixvec->second[i]]]=(seqfreq.count(ix2seq[ixvec->second[i]])>0?seqfreq[ix2seq[ixvec->second[i]]]+1:1);
     }
     for (unsigned int i=0;i<ixvec->second.size();i++) {
       feature[ROW(ixvec->second[i])+8]=log((float)seqfreq.size());
       feature[ROW(ixvec->second[i])+9]=log((float)seqfreq[ix2seq[ixvec->second[i]]]);
       feature[ROW(ixvec->second[i])+10]=f1;
     }
     
  }            
  cout << "Read File" << endl;
  delete str;
}


