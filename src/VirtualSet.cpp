#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <string>
#include <math.h>
using namespace std;
#include "DataSet.h"
#include "VirtualSet.h"

VirtualSet::VirtualSet()
{
   feature = NULL;
   n_examples=0;
   sqtFN = "";
}


VirtualSet::VirtualSet(const VirtualSet& mother,int fold,int ix)
{   bool *normalizedFlag;
    
    assert(ix<fold);
	label=mother.label;
	double nex=mother.n_examples/((double)fold);
	int start=(int)round(nex*ix);
	int end=(int)round(nex*(ix+1));
	n_examples=end-start;
    double * f = mother.feature;
	feature=&f[DataSet::rowIx(start)];
	intra=mother.intra;
}

VirtualSet::~VirtualSet()
{
  feature=NULL;
  intra=NULL;
}

double * VirtualSet::getNext(int& pos) {
  pos++;
  if (pos<0)
    pos=0;
  if(pos>=getSize())
    return NULL;
  return &feature[DataSet::rowIx(pos)];
}

const double * VirtualSet::getFeatures(const int pos) {
  return &feature[DataSet::rowIx(pos)];
}

bool VirtualSet::getGistDataRow(int & pos,string &out){
  ostringstream s1;
  double * feature = NULL;
//  while (!feature || charge[pos] !=2)  { //tmp fix
  if ((feature = getNext(pos)) == NULL) return false; 
//  }
  s1 << pos;
  for (int ix = 0;ix<DataSet::getNumFeatures();ix++) {
    s1 << '\t' << feature[ix];
  }
  s1 << endl;
  out = s1.str();
  return true;
}

void VirtualSet::print_features() {
   for(int i=0;i<getSize();i++) {
	   for(int j=0;j<DataSet::getNumFeatures();j++) {
	      cout << j+1 << ":" << feature[DataSet::rowIx(i)+j] << " ";
	   }
	   cout << endl;
   }
}

void VirtualSet::print_10features() {
   cerr << DataSet::getFeatureNames() << endl;
   for(int i=0;i<10;i++) {
       for(int j=0;j<DataSet::getNumFeatures();j++) {
          cerr << feature[DataSet::rowIx(i)+j] << "\t";
       }
       cerr << endl;
   }
}
