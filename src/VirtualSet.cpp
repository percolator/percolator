#include <vector>
#include <string>
#include <math.h>
using namespace std;
#include "DataSet.h"
#include "VirtualSet.h"

VirtualSet::VirtualSet(DataSet &mother,int fold,int ix) : DataSet()
{
	label=mother.label;
	double nex=mother.n_examples/((double)fold);
	int start=round(nex*ix);
	int end=round(nex*(ix+1));
	n_examples=end-start;
	feature=&mother.feature[start];
}

VirtualSet::~VirtualSet()
{
  feature=NULL;
}
