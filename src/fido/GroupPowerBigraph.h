// Written by Oliver Serang 2009
// see license for more information

#ifndef _GroupPowerBigraph_H
#define _GroupPowerBigraph_H

#include "BasicGroupBigraph.h"
#include "StringTable.h"
#include "Array.h"
#include "Random.h"
#include "Model.h"
#include "Scores.h" // from Percolator

using namespace std;


class GroupPowerBigraph
{
  
public:
  
 GroupPowerBigraph(Scores* fullset,double __alpha, double __beta, double __gamma, 
		   bool __groupProteins = false, bool __noseparate = false , bool __noprune = false) :
  gm(__alpha,__beta,__gamma),
  groupProteins(__groupProteins), 
  noseparate(__noseparate), 
  noprune(__noprune),
  LOG_MAX_ALLOWED_CONFIGURATIONS(18),
  PsmThreshold(0.0),
  PeptideThreshold(1e-3),
  ProteinThreshold(1e-3),
  PeptidePrior(0.1)
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
  double getLogNumberStates() const;
  pair<Array<Array<string> >, Array<double> > getDescendingProteinsAndWeights() const;
  void setAlphaBetaGamma(double alpha, double beta, double gamma);
  Array<std::string> getSeveredProteins();
  void setMaxAllowedConfigurations(double max_conf);
  double getMaxAllowedConfigurations();
  void setPsmThreshold(double psm_threshold);
  double getPsmThreshold();
  void setPeptideThreshold(double peptide_threshold);
  double getPeptideThreshold();
  void setProteinThreshold(double protein_threshold);
  double getProteinThreshold();
  void setPeptidePrior(double peptide_prior);
  double getPeptidePrior();
  
  //NOTE to clone object
  //GroupPowerBigraph *clone();
  
private:

  void initialize();
  void getGroupProtNames();
  void read(Scores* fullset);
  Array<BasicBigraph> iterativePartitionSubgraphs(BasicBigraph & bb, double newPeptideThreshold );
  
  Model gm;
  Numerical zeroChecker;
  int numberClones;
  bool groupProteins;
  bool noseparate;
  bool noprune;
  double ProteinThreshold;
  double PeptideThreshold;
  double PsmThreshold;
  double PeptidePrior;
  double LOG_MAX_ALLOWED_CONFIGURATIONS;
  Array<std::string> severedProteins;
  Array<double> probabilityR;
  Array<Array<std::string> > groupProtNames;
  Array<BasicGroupBigraph> subgraphs;
};

ostream & operator <<(ostream & os, pair<double,double> rhs);

#endif