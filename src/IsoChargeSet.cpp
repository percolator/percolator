#include "IsoChargeSet.h"

IsoChargeSet::IsoChargeSet(int c) {charge=c;}

IsoChargeSet::~IsoChargeSet()
{
}

const double * const IsoChargeSet::getNext(int& setPos,int& ixPos) {
  DataSet set = (*pSet)[setPos];
  double * features = set.getNext(charge,ixPos);
  if (features) return features;
  if (++setPos>=((signed int)pSet->size()))
    return NULL;
  ixPos=-1;
  set = (*pSet)[setPos];
  return set.getNext(charge,ixPos);
}

int const IsoChargeSet::getLabel(int *setPos) {
  return (*pSet)[*setPos].getLabel();
}

void IsoChargeSet::setSet(vector<DataSet> *set){
	pSet=set;
	int i=-1,j=0;
	while(getNext(i,j))
	  n_points++;
	norm.setSet(this);
}
