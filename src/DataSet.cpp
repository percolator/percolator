#include <stdio.h>
#include <iostream>
#include <set>
#include <map>
#include <utility>
#include "DataSet.h"
//using namespace std;



DataSet::DataSet()
{
   delim.assign(" \t\n");
}

DataSet::~DataSet()
{
	if (feature) {
		for(int i=0;i<n_feature;i++)
	   		delete feature[i];
		delete feature;
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
	      cout << j+1 << ":" << feature[i][j] << " ";
	   }
	   cout << endl;
   }
}

double * DataSet::getNext(const int c, int& pos) {
  pos++;
  if (pos<0)
    pos=0;
  while((pos<n_feature) && (charge[pos]!=c))
    (pos)++;
  if(pos>=n_feature)
    return NULL;
  return feature[pos];
}

void DataSet::read_sqt(char* fname) {
  FILE *fp1;
  int n = 0;
  size_t len1 = 0;
  vector<string> fields(30,"");

  char * str = (char *) calloc(1023,sizeof(char));

  if ((fp1=fopen(fname,"r"))==NULL) {printf("Could not open file %s\n",fname);}
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
  feature = new double*[n];
  ids.resize(n,"");
  charge.resize(n,0);
  vector<string> ix2seq(n,"");
  for(int i=0;i<n;i++) {
    feature[i]=new double[NUM_FEATURE];
  }
  n_feature=n;
  int ix=-1,gotL = 1;
  double mass;
  string seq;
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
    if (str[0]=='M' && !gotL) {
      line2fields(str,&fields);
      feature[ix][0]=atof(fields[2].data());
      feature[ix][1]=mass - atof(fields[3].data());
      feature[ix][2]=atof(fields[4].data());
      feature[ix][3]=atof(fields[5].data());
      feature[ix][4]=atof(fields[6].data());
      feature[ix][5]=atof(fields[7].data())/atof(fields[8].data());
      ix2seq[ix].assign(fields[9]);
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
       feature[ixvec->second[i]][6]=f1;
       feature[ixvec->second[i]][7]=seqfreq[ix2seq[ixvec->second[i]]];
     }
     
  }            
  cout << "Read File" << endl;
  delete str;
}

int main() {
  char *forwardFN = "/var/noble/data/tandem-ms/maccoss/rphase/2006-05-06/forward/050606-pps-2-01-forward-no-norm.sqt";
  char *randomFN = "/var/noble/data/tandem-ms/maccoss/rphase/2006-05-06/random/050606-pps-2-01-random-no-norm.sqt";
  char *random2FN = "/var/noble/data/tandem-ms/maccoss/rphase/2006-05-06/random2/050606-pps-2-01-random2-nonorm.sqt";

  cout << endl << " ReadData\n";

  DataSet * d = new DataSet(); 
  d->read_sqt(forwardFN);
  d->print_features();
  d->read_sqt(randomFN);
  d->read_sqt(random2FN);
  delete d;
  
}

