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
* PSM - ?
* ROC - Receiver Operating Characteristic
*/
class FidoInterface : public ProteinProbEstimator
{
      
    /** FIDO PARAMETERS **/

    /** compute peptide level prior probability instead of using default = 0.1 **/
    const static bool computePriors = false;
    /** threshold used for fido to remove poor PSMs **/
    const static double psmThreshold;
    const static double reduced_psmThreshold;
    /** threshold used for fido to classify a peptide as very low confidence **/
    const static double peptideThreshold;
    const static double reduced_peptideThreshold;
    /** threshold used for fido to classify a protein as very low confidence and prune it **/
    const static double proteinThreshold; 
    const static double reduced_proteinThreshold;
    /** default value for peptide prior probability used in fido to compute the peptide likehood **/
    const static double peptidePrior; 
    /** number of maximum of tree configurations allowed in fido **/
    const static double max_allow_configurations;
    /** allow the presence of peptides with the same sequence but different label (target/decoy) **/
    const static bool allow_multiple_labeled_peptides = false;
    
    /** GRID SEARCH PARAMETERS **/
    
    /** value that balances the objective function equation (lambda * rocR) - (1-lambda) * (fdr_mse) **/
    const static double lambda;
    
    /** number of false positives allowed while estiaming the ROC curve score **/
    /** if updateRocN is true the N value will be estimated automatically according to the number of FP found at a certain threshold **/
    const static unsigned default_rocN = 50;
    const static bool updateRocN = true;
    
    /** activate the optimization of the parameters to see the best boundaries**/
    const static bool optimize = false;
   

  public:

    FidoInterface(double __alpha = -1,double __beta = -1,double __gamma = -1,bool __nogroupProteins = false, 
       bool __noseparate = false, bool __noprune = false, unsigned __depth = 0,bool __reduceTree = false, 
       bool __truncate = true, double mse_threshold = 0.1,bool tiesAsOneProtein = false, bool usePi0 = false, 
       bool outputEmpirQVal = false, std::string decoyPattern = "random",bool trivialGrouping = false);

    virtual ~FidoInterface();
    
    void run();
    
    void computeProbabilities();
    void computeProbabilitiesFromFile(ifstream &fin);
    
    double getGamma(){return gamma;};
    double getBeta(){return beta;};
    double getAlpha(){return alpha;};
    
    string printCopyright();

  private:

    /** fido extra functions to do the grid search for parameters alpha,betha and gamma **/

    void getROC_AUC(const std::vector<std::vector<string> > &names,
         const std::vector<double> &probabilities, double &auc);
    
    void getEstimated_and_Empirical_FDR(const std::vector<std::vector<string> > &names,
            const std::vector<double> &probabilities,
            std::vector<double> &empq,
            std::vector<double> &estq);
    
    void getFDR_MSE(const std::vector<double> &estFDR, const std::vector<double> &empFDR,double &mse);
    
    void gridSearch();
    void gridSearchOptimize();
    
    double gamma;
    double alpha;
    double beta;
    double threshold;
    bool reduceTree;
    bool truncate;
    unsigned int depth;
    mutable unsigned int rocN;
    bool nogroupProteins;
    bool noseparate;
    bool noprune;
    bool dogridSearch;
    bool trivialGrouping;
    GroupPowerBigraph* proteinGraph;
};

#endif // FIDOINTERFACE_H
