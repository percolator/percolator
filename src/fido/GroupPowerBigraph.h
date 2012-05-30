// Written by Oliver Serang 2009
// see license for more information

#ifndef _GroupPowerBigraph_H
#define _GroupPowerBigraph_H
#include "ProteinIdentifier.h"
#include "BasicGroupBigraph.h"
#include "PivdoSplitter.h"
#include "Scores.h" // from Percolator

using namespace std;


class GroupPowerBigraph : public ProteinIdentifier
{
  
public:
  
 static double LOG_MAX_ALLOWED_CONFIGURATIONS;

 GroupPowerBigraph(Scores* fullset,double __alpha, double __beta, double __gamma, 
		   bool __groupProteins = false, bool __noseparate = false , bool __noprune = false) :
  ProteinIdentifier(),
  gm(__alpha,__beta,__gamma),
  groupProteins(__groupProteins), 
  noseparate(__noseparate), 
  noprune(__noprune)
      {
	setAlphaBetaGamma(__alpha, __beta, __gamma);
	read(fullset);
      }
  
  ~GroupPowerBigraph();
  Array<double> proteinProbs();
  void printProteinWeights() const;
  std::multimap<double, std::vector<std::string> > getProteinProbsPercolator() const;
  void getProteinProbsAndNames(std::vector<std::vector<std::string> > &names, std::vector<double> &probs) const;
  void getProteinProbs();
  Array<string> peptideNames() const;
  Array<BasicBigraph> iterativePartitionSubgraphs(BasicBigraph & bb, double newPeptideThreshold );
  double getLogNumberStates() const;
  double logLikelihoodAlphaBetaGivenD( const GridModel & myGM ) const;
  void outputPivdo(ostream & os) const;
  pair<Array<Array<string> >, Array<double> > getDescendingProteinsAndWeights() const;
  void setAlphaBetaGamma(double alpha, double beta, double gamma);
  Array<std::string> getSeveredProteins();
  //NOTE to clone object
  //GroupPowerBigraph *clone();
protected:

  void initialize();
  Array<double> probabilityXAndEGivenD();
  Array<double> probabilityXGivenD();
  
  void getGroupProtNames();
  void read(Scores* fullset);

  Model gm;
  Numerical zeroChecker;
  int numberClones;
  bool groupProteins;
  bool noseparate;
  bool noprune;
  Array<std::string> severedProteins;
  Array<double> probabilityR;
  Array<Array<std::string> > groupProtNames;
  Array<BasicGroupBigraph> subgraphs;
};


ostream & operator <<(ostream & os, pair<double,double> rhs);


#endif

