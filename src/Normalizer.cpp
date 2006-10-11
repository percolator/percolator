#include <vector>
using namespace std;
#include "Normalizer.h"
#include "StdvNormalizer.h"
#include "UniNormalizer.h"

Normalizer::Normalizer()
{
}

Normalizer::~Normalizer()
{
}

int Normalizer::subclass_type = STDV;

Normalizer * Normalizer::getNew(){
	if (subclass_type == UNI)
	  return new UniNormalizer();
	return new StdvNormalizer();
}

void Normalizer::setType(int type){
	assert(type==UNI || type==STDV);
	subclass_type = type;
}
