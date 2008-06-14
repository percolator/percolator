#ifndef PSMDESCRIPTION_H_
#define PSMDESCRIPTION_H_
#include <set>
#include <string>
using namespace std;

class PSMDescription
{
public:
  PSMDescription();
  virtual ~PSMDescription();
  void clear() {proteinIds.clear();}

  double q,pep;
  double * features;
  double retentionTime,massDiff,pI;
  unsigned int scan; 
  set<string> proteinIds;
  string peptide;  
};

#endif /*PSMDESCRIPTION_H_*/
