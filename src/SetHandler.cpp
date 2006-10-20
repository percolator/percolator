#include<fstream>
#include <vector>
#include <string>
using namespace std;
#include "DataSet.h"
#include "Normalizer.h"
#include "SetHandler.h"

SetHandler::SetHandler() {
//	charge=c;
	norm=Normalizer::getNew();
    n_examples=0;
    n_pos=0;
    n_neg=0;
}

SetHandler::~SetHandler()
{
	delete norm;
	norm=NULL;
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
	norm->setSet(this);
}

void SetHandler::gistWrite(const string & fileNameTrunk) {
  string dataFN = fileNameTrunk + ".data";
  string labelFN = fileNameTrunk + ".label";
  ofstream dataStream(dataFN.data(),ios::out);
  ofstream labelStream(labelFN.data(),ios::out);
  labelStream << "SpecId\tLabel\n"; 
  dataStream << "SpecId\tF1\tF2\tF3\tF4\tF5\tF6\tF7\tF8\tF9\tF10\tF11\tF12\tF13\tF14\n";
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
