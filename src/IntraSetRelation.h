#ifndef INTRASETRELATION_H_
#define INTRASETRELATION_H_

class IntraSetRelation
{
public:
  IntraSetRelation();
  virtual ~IntraSetRelation();
  void registerRel(string pep, set<string> &prot);
  int getNumProt(string & prot) {return (numProteins.count(prot)?0:numProteins[prot]);}
  int getNumPep(string & pep) {return (numPeptides.count(pep)?0:numPeptides[pep]);}
  int getPepSites(set<string> &prot);
protected:
  map<string,int> numProteins;
  map<string,int> numPeptides;
  map<string,set<string> > prot2pep;
};

#endif /*INTRASETRELATION_H_*/
