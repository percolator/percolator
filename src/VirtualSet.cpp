#include <vector>
#include <set>
#include <string>
#include <math.h>
using namespace std;
#include "DataSet.h"
#include "VirtualSet.h"

VirtualSet::VirtualSet(DataSet & mother,int fold,int ix) : DataSet()
{
    assert(ix<fold);
	label=mother.getLabel();
	double nex=mother.getSize()/((double)fold);
	int start=(int)round(nex*ix);
	int end=(int)round(nex*(ix+1));
	n_examples=end-start;
    double * f = mother.getFeature();
	feature=&f[DataSet::rowIx(start)];
}

VirtualSet::~VirtualSet()
{
  feature=NULL;
}
