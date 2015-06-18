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

/*
* GroupPowerBigraph represents a collection of bigraphs of PSMs on one side and 
*   protein groups on the other. The power set enumerates all possible combinations 
*   of present and absent proteins.
* 
*/
class GroupPowerBigraph {
  
 public:
  
  GroupPowerBigraph(double alpha, double beta, double gamma, 
      bool noClustering = false, bool noPartitioning = false , 
      bool noPruning = false, bool trivialGrouping = false) :
        params_(alpha, beta, gamma), noPartitioning_(noPartitioning), 
        noClustering_(noClustering), noPruning_(noPruning),
        LOG_MAX_ALLOWED_CONFIGURATIONS(18),
        psmThreshold_(0.0), peptideThreshold_(1e-3),
        proteinThreshold_(1e-3), peptidePrior_(0.1),
        trivialGrouping_(trivialGrouping) {}
  ~GroupPowerBigraph();
  
  Array<double> proteinProbs();
  void printProteinWeights() const;
  void getProteinProbsPercolator(std::multimap<double, std::vector<std::string> > &pepProteins) const;
  void getProteinProbsAndNames(std::vector<std::vector<std::string> > &names, std::vector<double> &probs) const;
  void getProteinNames(std::vector<std::vector<std::string> > &names) const;
  void getProteinProbs();
  Array<string> peptideNames() const;
  double getLogNumberStates() const;
  pair<Array<Array<string> >, Array<double> > getDescendingProteinsAndWeights() const;
  
  void setAlphaBetaGamma(double alpha, double beta, double gamma) {
    params_.setAlphaBetaGamma(alpha, beta, gamma);
  }
  
  Array<std::string> getSeveredProteins() { return severedProteins_; }
  
  void setMaxAllowedConfigurations(double max_conf) {
    LOG_MAX_ALLOWED_CONFIGURATIONS = max_conf;
  }
  double getMaxAllowedConfigurations() const { 
    return LOG_MAX_ALLOWED_CONFIGURATIONS; 
  }
  
  void setPsmThreshold(double t) { psmThreshold_ = t; }
  double getPsmThreshold() const { return psmThreshold_; }
  void setPeptideThreshold(double t) { peptideThreshold_ = t; }
  double getPeptideThreshold() const { return peptideThreshold_; }
  void setProteinThreshold(double t) { proteinThreshold_ = t; }
  double getProteinThreshold() const { return proteinThreshold_; }
  
  void setPeptidePrior(double p) { peptidePrior_ = p; }
  double getPeptidePrior() const { return peptidePrior_; }
  
  void setTrivialGrouping(bool b) { trivialGrouping_ = b; }
  bool getTrivialGrouping() const { return trivialGrouping_; }
  
  void setNoClustering(bool b) { noClustering_ = b; }
  bool getNoClustering() const { return noClustering_; }
  
  void setNoPruning(bool b) { noPruning_ = b; }
  bool getNoPruning() const { return noPruning_; }
  
  void setNoPartitioning(bool b) { noPartitioning_ = b; }
  bool getNoPartitioning() const { return noPartitioning_; }
  
  void setMultipleLabeledPeptides(bool b) { addPeptideDecoyLabel_ = b; }
  bool getMultipleLabeledPeptides() const { return addPeptideDecoyLabel_; }
  
  void read(Scores* fullset);
  void read(istream & is);
  //NOTE to clone object
  //GroupPowerBigraph *clone();
private:
  void initialize(BasicBigraph& basicBigraph);
  void getGroupProtNames();
  
  Array<BasicBigraph> iterativePartitionSubgraphs(BasicBigraph & bb, double newPeptideThreshold );
  
  /* struct with alpha, beta and gamma parameter (Fig 2 in Serang et al. 2010) */
  Model params_;
  /* turns off partitioning, clustering or pruning of the graph (Fig 3 in Serang et al. 2010) */
  bool noPartitioning_;
  bool noClustering_;
  bool noPruning_;
  /* adds an asterisk (*) to decoy peptide sequences (see BasicBigraph.h:41) */
  bool addPeptideDecoyLabel_;
  /* log2 of the maximum number of configurations to limit runtime and complexity */
  double LOG_MAX_ALLOWED_CONFIGURATIONS;
  /* prunes PSMs with probability below the threshold */
  double psmThreshold_;
  /* prunes peptides with probability below threshold */
  double peptideThreshold_;
  /* prunes proteins with best peptide probability below threshold */
  double proteinThreshold_;
  /* prior probability of a peptide being present ("E" in Serang et al. 2010) */
  double peptidePrior_;
  /* groups are either present or absent and cannot be partially present */
  bool trivialGrouping_;
  /* proteins that have no PSMs remaining after pruning */
  Array<std::string> severedProteins_;
  /* probabilities for each protein to be present ("R" in Serang et al. 2010) */
  Array<double> probsPresentProteins_;
  /* protein names grouped by groups from the clustering step */
  Array<Array<std::string> > groupProtNames_;
  /* subgraphs resulting from the partitioning and pruning steps */
  Array<BasicGroupBigraph> subgraphs_;
};

ostream & operator <<(ostream & os, pair<double,double> rhs);

#endif
