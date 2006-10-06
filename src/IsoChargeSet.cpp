#include "DataSet.h"
#include "Normalizer.h"
#include "IsoChargeSet.h"

IsoChargeSet::IsoChargeSet(int c) {charge=c;norm=new Normalizer();}

IsoChargeSet::~IsoChargeSet()
{
	delete norm;
	norm=NULL;
}

const double * IsoChargeSet::getNext(int& setPos,int& ixPos) {
  DataSet *set = &(*pSet)[setPos];
  double * features = set->getNext(charge,ixPos);
  if (features) return features;
  if (++setPos>=((signed int)pSet->size()))
    return NULL;
  ixPos=-1;
  set = &(*pSet)[setPos];
  return set->getNext(charge,ixPos);
}

int const IsoChargeSet::getLabel(int *setPos) {
  return (*pSet)[*setPos].getLabel();
}

void IsoChargeSet::setSet(vector<DataSet> *set){
	pSet=set;
	int i=0,j=-1;
	while(getNext(i,j))
	  n_points++;
	norm->setSet(this);
}
