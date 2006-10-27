#include<iostream>
#include<fstream>
#include <vector>
#include <string>
using namespace std;
#include "DataSet.h"
#include "Normalizer.h"
#include "SetHandler.h"
#include "Scores.h"
#include "Globals.h"

SetHandler::SetHandler() {
//	charge=c;
//	norm=Normalizer::getNew();
    n_examples=0;
    n_pos=0;
    n_neg=0;
    labels = NULL;
    c_vec = NULL;
}

SetHandler::~SetHandler()
{
//	if (norm)
//      delete norm;
//	norm=NULL;
    if (labels)
      delete [] labels;
    labels=NULL;
    if (c_vec)
      delete [] c_vec;
    c_vec=NULL;
}

void SetHandler::generateTrainingSet(const double fdr,const double cpos, const double cneg,const Scores & sc) {
  double tp=0,fp=0;
  unsigned int ix=0;
  examples.clear();
  bool underCutOff = true;
  vector<ScoreHolder>::const_iterator it;
  for(it=sc.begin();it!=sc.end();it++) {
    if (it->label==-1) fp++; else tp++;
    if (underCutOff && fdr<(fp/(tp+fp)))
      underCutOff=false;
    if (it->label==-1 || underCutOff) {
      examples.push_back(getFeatures(it->set,it->index));
      labels[ix]=it->label;
      c_vec[ix++]=(it->label!=-1?cpos:cneg);
    }
  }
}


const double * SetHandler::getNext(int& setPos,int& ixPos) {
  double * features = subsets[setPos]->getNext(ixPos);
  if (features) return features;
  if (++setPos>=((signed int)subsets.size()))
    return NULL;
  ixPos=-1;
  return subsets[setPos]->getNext(ixPos);
}

const double * SetHandler::getFeatures(const int setPos,const int ixPos) {
  return subsets[setPos]->getFeatures(ixPos);
}

int const SetHandler::getLabel(int setPos) {
  assert(setPos>=0 && setPos<(signed int)subsets.size());
  return subsets[setPos]->getLabel();
}



void SetHandler::setSet(vector<DataSet *> & pos,vector<DataSet *> &neg){
    subsets.clear();
    subsets.assign(pos.begin(),pos.end());
    subsets.insert(subsets.end(),neg.begin(),neg.end());
	n_examples=0;
    n_pos=0;
    n_neg=0;
	int i=0,j=-1;
	while(getNext(i,j)) {
	  n_examples++;
      if (getLabel(i)==-1) n_neg++;
      else n_pos++;
    }
    if(!labels) labels= new double[n_examples];
    if(!c_vec) c_vec = new double[n_examples];
    if (VERB>3) {
      int pos=0,neg=0;
      for (unsigned int i=0;i<subsets.size();i++) {
        if (subsets[i]->getLabel()==1) pos++; else neg++;
      }
      cerr << "Set up a SetHandler with " << pos << " positive DataSet:s and " << n_pos << " examples" << endl;
      cerr << "and " << neg << " negative DataSet:s and " << n_neg << " examples" << endl;
      if (VERB>4) {
        for (unsigned int i=0;i<subsets.size();i++) {
          cerr << "First 10 lines of " << i+1 << " set with " << subsets[i]->getLabel() << " label" << endl;
          subsets[i]->print_10features();
        }
      }
    }
}

void SetHandler::gistWrite(const string & fileNameTrunk) {
  string dataFN = fileNameTrunk + ".data";
  string labelFN = fileNameTrunk + ".label";
  ofstream dataStream(dataFN.data(),ios::out);
  ofstream labelStream(labelFN.data(),ios::out);
  labelStream << "SpecId\tLabel\n"; 
  dataStream << "SpecId\t" << DataSet::getFeatureNames();
  string str;
  for (int setPos=0;setPos< (signed int)subsets.size();setPos++) {
    int ixPos=-1;
    while (subsets[setPos]->getGistDataRow(ixPos,str)) {
      dataStream << str;
      labelStream << str.substr(0,str.find('\t')+1) << (subsets[setPos]->getLabel()==-1?-1:+1) << endl; 
    }
  }    
  dataStream.close();
  labelStream.close();
  
}
