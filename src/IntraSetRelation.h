#ifndef INTRASETRELATION_H_
#define INTRASETRELATION_H_

class IntraSetRelation
{
public:
  IntraSetRelation();
  virtual ~IntraSetRelation();
  void registerRel(const string pep, const vector<string> &prot);
  int getNumProt(const string & prot) {return (numProteins.count(prot)?0:numProteins[prot]);}
  int getNumPep(const string & pep) {return (numPeptides.count(pep)?0:numPeptides[pep]);}
  int getPepSites(const vector<string> &prot);
protected:
  map<string,int> numProteins;
  map<string,int> numPeptides;
  map<string,set<string> > prot2pep;
};

#endif /*INTRASETRELATION_H_*/
