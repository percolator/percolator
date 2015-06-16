/*******************************************************************************
 Copyright 2006-2012 Lukas Käll <lukas.kall@scilifelab.se>

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

 *******************************************************************************/

#ifndef FIDOINTERFACE_H
#define FIDOINTERFACE_H
#include "ProteinProbEstimator.h"
#include "GroupPowerBigraph.h"

/*
* FidoInterface is a class that computes probabilities and statistics based
* on provided proteins from the set of scored peptides from percolator. It
* uses number of additional parameters to setup the calculations.
*
* Here are some usefull abbreviations:
* AUC - Area Under Curve
* FDR - False Discovery Rate
* FIDO - ?
* MSE - Mean Squared Error
* PSM - Peptide Spectrum Match
* ROC - Receiver Operating Characteristic
*
*/
class FidoInterface : public ProteinProbEstimator {
      
  /** FIDO PARAMETERS **/

  /** compute peptide level prior probability instead of using default = 0.1 **/
  const static bool kComputePriors = false;
  /** threshold to remove poor PSMs **/
  const static double kPsmThreshold;
  /** threshold determine if a graph can be split into smaller subgraphs at low confidence PSMs **/
  const static double kPeptideThreshold;
  /** default value for peptide prior probability used in fido to compute the peptide likehood **/
  const static double kPeptidePrior; 
  /** log2 number of maximum of tree configurations allowed in fido **/
  const static double LOG_MAX_ALLOWED_CONFIGURATIONS;
  /** allow the presence of peptides with the same sequence but different label (target/decoy) **/
  const static bool kAddPeptideDecoyLabel = false;

  /** GRID SEARCH PARAMETERS **/

  /** value that balances the objective function equation (kObjectiveLambda * rocR) - (1-kObjectiveLambda) * (fdr_mse) **/
  const static double kObjectiveLambda;
  /** number of false positives allowed while estimating the ROC curve score **/
  /** if kUpdateRocN is true the N value will be estimated automatically according to the number of FP found at a certain threshold (maxN = 500) **/
  const static unsigned kDefaultRocN = 50u;
  const static bool kUpdateRocN = true;
  /** activate the optimization of the parameters to see the best boundaries**/
  const static bool kOptimizeParams = false;

 public:
  FidoInterface(double alpha = -1, double beta = -1, double gamma = -1, 
    bool noPartitioning = false, bool noClustering = false, 
    bool noPruning = true, 
    unsigned gridSearchDepth = 0u, double gridSearchThreshold = 0.0, 
    double proteinThreshold = 0.01, double mse_threshold = 0.1, 
    bool tiesAsOneProtein = true, bool usePi0 = false, 
    bool outputEmpirQVal = false, std::string decoyPattern = "random", 
    bool trivialGrouping = true);
  virtual ~FidoInterface();
  
  void run();
  void computeProbabilities(const std::string& fname = "");
  
  std::ostream& printParametersXML(std::ostream &os);
  string printCopyright();

 private:
  /** FIDO PARAMETERS **/
  
  /* bigraph used by fido to represent connections between PSMs and proteins */
  GroupPowerBigraph* proteinGraph_;
  /* alpha, beta and gamma parameter (Fig 2 in Serang et al. 2010) */
  double alpha_, beta_, gamma_;
  /* turns off optimization steps (Fig 3 in Serang et al. 2010) */
  bool noPartitioning_, noClustering_, noPruning_;
  /* turns off pruning of edges to proteins with only low confident PSMs */
  double proteinThreshold_;
  /* protein groups are either present or absent and cannot be partially present */
  bool trivialGrouping_;
  
  /** GRID SEARCH PARAMETERS **/
  
  /* sets the density of the grid search: 0 is the coarsest and 2 the finest */
  unsigned int gridSearchDepth_;
  /* determines if grid search is activated by the user */
  bool doGridSearch_;
  /* uses strict thresholds to create a sparse graph for the grid search */ 
  double gridSearchThreshold_;
  /* threshold in MSE estimation */
  double mseThreshold_;
  /* threshold for ROC AUC estimation */
  mutable unsigned int rocN_;
  
  /** fido extra functions to do the grid search for parameters alpha,beta and gamma **/
  void getROC_AUC(const std::vector<std::vector<string> > &names,
       const std::vector<double> &probabilities, double &auc);
  
  void getEstimated_and_Empirical_FDR(const std::vector<std::vector<string> > &names,
          const std::vector<double> &probabilities,
          std::vector<double> &empq,
          std::vector<double> &estq);
  
  void getFDR_MSE(const std::vector<double> &estFDR, const std::vector<double> &empFDR,double &mse);
  
  void gridSearch();
  void gridSearchOptimize();
  void gridSearch(std::vector<double>& alpha_search, 
                  std::vector<double>& beta_search, 
                  std::vector<double>& gamma_search);
  double calcObjective(double alpha, double beta, double gamma);  
  
};

/* FDRCalculator is a helper class for FidoInterface to facilitate the two 
*    options for FDR calculation (tiesAsOneProtein on/off).
*
*/
class FDRCalculator {
 public:
  FDRCalculator(double targetDecoyRatio, double pi0, bool countDecoyQvalue) : 
      fpCount_(0.0), tpCount_(0.0), totalFDR_(0.0), 
      previousEmpQ_(0.0), previousEstQ_(0.0),
      pi0_(pi0), targetDecoyRatio_(targetDecoyRatio),
      countDecoyQvalue_(countDecoyQvalue), rocN_(50u) {}
  
  unsigned getRocN() const { return rocN_; }
  double getPreviousEstQ() const { return previousEstQ_; }
  
  void calcFDRs(double fpChange, double tpChange, double prob,
               std::vector<double>& empq, std::vector<double>& estq);
 private:
  double fpCount_, tpCount_, totalFDR_;
  double pi0_, targetDecoyRatio_;
  double previousEmpQ_, previousEstQ_;
  bool countDecoyQvalue_;
  unsigned int rocN_;  
};

#endif // FIDOINTERFACE_H
