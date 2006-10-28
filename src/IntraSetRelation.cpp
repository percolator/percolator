#include <set>
#include <vector>
#include <map>
#include <string>
using namespace std;
#include "IntraSetRelation.h"

IntraSetRelation::IntraSetRelation()
{
}

IntraSetRelation::~IntraSetRelation()
{
}
void IntraSetRelation::registerRel(const string pep, const vector<string> &prot){
  if(numPeptides.count(pep)==0) numPeptides[pep]=1;
  else numPeptides[pep]++;
  vector<string>::iterator iProt;
  for(iProt=prot.begin();iProt!=prot.end();iProt++) {
    if(numProteins.count(*iProt)==0) numProteins[*iProt]=1;
    else numProteins[*iProt]++;
    prot2pep[*iProt].insert(pep);
  }
}

int IntraSetRelation::getPepSites(const vector<string> &prot) {
  int maxn=0;
  vector<string>::iterator iProt;
  for(iProt=prot.begin();iProt!=prot.end();iProt++) {
  	int n=prot2pep[*iProt].size();
  	if (n>maxn) maxn=n;
  }
  return maxn;
}

