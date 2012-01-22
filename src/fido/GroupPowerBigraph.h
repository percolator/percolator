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
		   bool __groupProteins = true, bool __noseparate = false , bool __noprune = false) :
  ProteinIdentifier(),
  gm(RealRange(__alpha,1,__alpha),RealRange(__beta,1,__beta),__gamma),
  groupProteins(__groupProteins), 
  noseparate(__noseparate), 
  noprune(__noprune),
  sumLogLikelihoodOverAllAlphaBetaCachedFunctor( & GroupPowerBigraph::sumLogLikelihoodOverAllAlphaBeta, "sumLogLikelihoodOverAllAlphaBeta")
      {
	setAlphaBetaGamma(__alpha, __beta, __gamma);
	read(fullset);
      }
      
  Array<double> proteinProbs(const GridModel & myGM);
  Array<double> proteinProbsOverAllAlphaBeta();
  void printProteinWeights() const;
  std::multimap<double, std::vector<std::string> > getProteinProbsPercolator() const;
  void getProteinProbs();
  void getProteinProbsOverAllAlphaBeta();
  Array<string> peptideNames() const;
  double probabilityAlphaBetaGivenD();
  void scanAlphaBeta();
  Array<BasicBigraph> iterativePartitionSubgraphs(BasicBigraph & bb, double newPeptideThreshold );
  double getLogNumberStates() const;
  double likelihoodAlphaBetaGivenD( const GridModel & myGM ) const;
  double logLikelihoodAlphaBetaGivenD( const GridModel & myGM ) const;
  // only pass this so that it can be cached easily 
  // even though you are passing a member variable
  double sumLogLikelihoodOverAllAlphaBeta( const GridModel & myGM ) const;
  // only pass this so that it can be cached easily 
  // even though you are passing a member variable
  double probabilityAlphaBetaGivenD( const GridModel & myGM ) const;
  void gridScan() const;
  void outputPivdo(ostream & os) const;
  pair<Array<Array<string> >, Array<double> > getDescendingProteinsAndWeights() const;
  void setAlphaBetaGamma(double alpha, double beta, double gamma);
  Array<std::string> getSeveredProteins();

protected:

  void initialize();
  Array<double> probabilityXAndEGivenD();
  Array<double> probabilityXGivenD();
  void getProteinProbsGivenAlphaBeta();
  // cached functors
  LastCachedMemberFunction<GroupPowerBigraph, double, GridModel> sumLogLikelihoodOverAllAlphaBetaCachedFunctor;
  void getGroupProtNames();
  void read(Scores* fullset);

  GridModel gm;
  double alpha,beta,gamma;
  Numerical zeroChecker;
  int numberClones;
  bool groupProteins;
  bool noseparate;
  bool noprune;
  //NOTE make protected
  Array<std::string> severedProteins;
  Array<double> probabilityR;
  Array<Array<std::string> > groupProtNames;
  Array<BasicGroupBigraph> subgraphs;
};


ostream & operator <<(ostream & os, pair<double,double> rhs);


#endif

