#ifndef PSMDESCRIPTION_H_
#define PSMDESCRIPTION_H_
#include <set>
#include <map>
#include <vector>
#include <string>
using namespace std;

class PSMDescription
{
public:
  PSMDescription();
  virtual ~PSMDescription();
  void clear() {proteinIds.clear();}
  double * getFeatures() {return features;}
  double * getRetentionFeatures() {return retentionFeatures;}
  string& getPeptide() {return peptide;}
  double getUnnormalizedRetentionTime() { return unnormalize(retentionTime);}
  static void setRetentionTime(vector<PSMDescription>& psms, map<int,double>& scan2rt);
  static double unnormalize(double normalizedTime);

  static double normDiv,normSub;

  double q,pep;
  double * features;
  double * retentionFeatures;
  double retentionTime,predictedTime,massDiff,pI;
  unsigned int scan;
  string id;
  string peptide;
  set<string> proteinIds;
};

#endif /*PSMDESCRIPTION_H_*/
