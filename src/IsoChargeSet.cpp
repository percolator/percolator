#include<fstream>
#include <vector>
#include <string>
using namespace std;
#include "DataSet.h"
#include "Normalizer.h"
#include "IsoChargeSet.h"

IsoChargeSet::IsoChargeSet() {
//	charge=c;
	norm=Normalizer::getNew();
    n_points=0;
    pSet=NULL;
}

IsoChargeSet::~IsoChargeSet()
{
	delete norm;
	norm=NULL;
}

const double * IsoChargeSet::getNext(int& setPos,int& ixPos) {
  DataSet *set = (*pSet)[setPos];
  double * features = set->getNext(0,ixPos);
  if (features) return features;
  if (++setPos>=((signed int)pSet->size()))
    return NULL;
  ixPos=-1;
  set = (*pSet)[setPos];
  return set->getNext(0,ixPos);
}

int const IsoChargeSet::getLabel(int *setPos) {
  assert(pSet!=NULL);
  assert(*setPos>=0 && *setPos<(signed int)pSet->size());
  return (*pSet)[*setPos]->getLabel();
}

void IsoChargeSet::setSet(vector<DataSet *> *set){
	n_points=0;
	pSet=set;
	int i=0,j=-1;
	while(getNext(i,j))
	  n_points++;
	norm->setSet(this);
}

void IsoChargeSet::gistWrite(const string & fileNameTrunk) {
  string dataFN = fileNameTrunk + ".data";
  string labelFN = fileNameTrunk + ".label";
  ofstream dataStream(dataFN.data(),ios::out);
  ofstream labelStream(labelFN.data(),ios::out);
  labelStream << "SpecId\tLabel\n"; 
  dataStream << "SpecId\tF1\tF2\tF3\tF4\tF5\tF6\tF7\tF8\tF9\tF10\tF11\tF12\tF13\tF14\n";
  string str;
  for (int setPos=0;setPos< (signed int)pSet->size();setPos++) {
    int ixPos=-1;
    while ((*pSet)[setPos]->getGistDataRow(ixPos,str)) {
      dataStream << str;
      labelStream << str.substr(0,str.find('\t')+1) << ((*pSet)[setPos]->getLabel()==-1?-1:+1) << endl; 
    }
  }    
  dataStream.close();
  labelStream.close();
  
}
