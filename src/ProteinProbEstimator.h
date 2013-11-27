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
 
#ifndef PROTEINPROBESTIMATOR_H_
#define PROTEINPROBESTIMATOR_H_

#include "PosteriorEstimator.h"
#include "Globals.h"
#include <boost/algorithm/string.hpp>
#include <functional>
#include <numeric>
#include <iterator>
#include "ProteinFDRestimator.h"
#include "Protein.h"
#include <vector>
#include <boost/assign/list_of.hpp>
#include <math.h>
#include <cmath>
#include "percolator_in.hxx"
#include "Scores.h"

/** set of helper functions to sort data structures and some operations overloaded **/

struct IntCmpProb {
    bool operator()(const std::pair<const std::string,Protein*> &lhs, const std::pair<const std::string,Protein*> &rhs) {
        return 
	   (  (lhs.second->getPEP() < rhs.second->getPEP())
	   || ( (lhs.second->getPEP() == rhs.second->getPEP()) && (lhs.second->getQ() < rhs.second->getQ()) )
	   || ( (lhs.second->getPEP() == rhs.second->getPEP()) && (lhs.second->getQ() == rhs.second->getQ())
	      && (lhs.second->getName() < rhs.second->getName()) )  
      );
    }
};
  
inline double myminfunc(double a, double b) 
{
  return a > b ? b : a;
}

struct RetrieveKey
{
  template <typename T>
  typename T::first_type operator()(T keyValuePair) const
  {
    return keyValuePair.first;
  }
};
    
struct RetrieveValue
{
  template <typename T>
  typename T::second_type operator()(T keyValuePair) const
  {
    return keyValuePair.second;
  }
};

/*
* ProteinProbEstimator is a class that computes probabilities and statistics based
* on provided proteins from the set of scored peptides from percolator.
*
* Here are some usefull abbreviations:
* Mayu - a software package for the analysis of (large) mass
*           spectrometry-based shotgun proteomics data sets.
*/
 
class ProteinProbEstimator {
  
  public:
    
    /** PROTEIN FDR ESTIMATOR PARAMETERS **/
    
    /* Default configuration (changeable by functions)
     * decoy prefix = random
     * number of bins = 10
     * target decoy ratio = 1.0
     * binning mode = equal deepth
     * correct identical sequences = true
     */
    
    const static bool correct_identical_sequences = true;
    const static bool binning_equal_deepth = true;
    const static double target_decoy_ratio = 1.0;
    const static unsigned number_bins = 10;
    
    
    /** GENERAL PARAMETERS **/
    
    /** compute global protein FDR using MAYU based implementation **/
    const static bool mayufdr = false;
    /** threshold used to estimate the protein FDR(pi0) **/
    const static double psmThresholdMayu = 0.90;
    /** whether to use the decoy prefix(faster) to check if a protein is decoy or not **/
    const static bool useDecoyPrefix = true;
    /** whether to count decoy proteins when estimated q values or not **/
    const static bool countDecoyQvalue = true;
    /** protein prior probability used to estimate the peptide prior probabilities **/
    const static double prior_protein = 0.5;
    
    
    /******************************************************************************************************************/
    
    
    ProteinProbEstimator(bool tiesAsOneProtein = false, bool usePi0 = false, 
			  bool outputEmpirQVal = false, std::string decoyPattern = "random");
    
    virtual ~ProteinProbEstimator();
    
    /** reads the proteins from the set of scored peptides from percolator **/
    bool initialize(Scores* fullset);
    
    /** start the protein probabilities tool**/
    virtual void run() = 0;
    
    /** initialize the estimation of the protein probabilities **/
    virtual void computeProbabilities() = 0;
    
    /** print out the tab delimited list of proteins to std::cerr or the screen */
    void printOut(const std::string &proteinFN, 
		   const std::string &proteinDecoyFN);
    
    /** initialize the estimation of the q values and p values **/
    void computeStatistics();
    
    /** MAYUS method for estimation of Protein FDR **/
    void computeFDR();
    
    /** write the list of proteins to the output file **/
    void writeOutputToXML(string xmlOutputFN, bool outputDecoys);

    /** Return the number of proteins whose q value is less or equal than the threshold given**/
    unsigned getQvaluesBelowLevel(double level);
    unsigned getQvaluesBelowLevelDecoy(double level);
   
    /** populate the list of proteins**/
    void setTargetandDecoysNames();
    
    /** return the data structure for the proteins **/
    std::map<const std::string,Protein*> getProteins();
    
    /** add proteins read from the database **/
    void addProteinDb(const percolatorInNs::protein &protein);
    
    /** print copyright of the author**/
    virtual string printCopyright() = 0;
    
    /**some getters **/
    
    double getPi0(){return pi0;};
    double getFDR(){return fdr;};
    
    /** Hack for fido :( **/
    virtual double getGamma() = 0;
    virtual double getBeta() = 0;
    virtual double getAlpha() = 0;
    
  protected:
   
    /** functions to count number of target and decoy proteins **/
    unsigned countTargets(const std::vector<std::string> &proteinList);
    unsigned countDecoys(const std::vector<std::string> &proteinList);
    bool isTarget(const std::string& proteinName);
    bool isDecoy(const std::string& proteinName);
    
     /** print a tab delimited list of proteins probabilities in a file or stdout**/
    void print(ostream& myout, bool decoy=false);
    
    /** function that extracts a list of proteins from the peptides that have a qvalue lower than psmThresholdMayu
     * this function is used to estimate the protein FDR**/
    void getTPandPFfromPeptides(double threshold, std::set<std::string> &numberTP, 
				  std::set<std::string> &numberFP);
       
    /** estimate prior probabilities for peptide level **/
    double estimatePriors();
    
    /** this function generates a vector of pair protein pep and label **/
    void getCombinedList(std::vector<std::pair<double , bool> > &combined);
    
     /** update the proteins with the computed qvalues and pvalues**/
    void updateProteinProbabilities();
    
    /** compute estimated qvalues from the PEP**/
    void estimateQValues();
    
    /** compute pvalues from the scored target/decoy proteins**/
    void estimatePValues(); 
    
    /** compute empirical qvalues from the target/decoy proteins**/
    void estimateQValuesEmp();
    
    /** compute pi0 from the set of pvalues**/
    double estimatePi0(const unsigned int numBoot = 100);
    
    /** variables **/
    std::set<string> truePosSet, falsePosSet;
    ProteinFDRestimator *fastReader;
    std::map<const std::string,Protein*> proteins;    
    std::multimap<double,std::vector<std::string> > pepProteins;
    std::map<std::string,std::pair<std::string,double> > targetProteins;
    std::map<std::string,std::pair<std::string,double> > decoyProteins;
    std::vector<double> qvalues;
    std::vector<double> qvaluesEmp;
    std::vector<double> pvalues;
    Scores* peptideScores;
    bool tiesAsOneProtein;
    bool usePi0;
    bool outputEmpirQVal;
    double pi0;
    double fdr;
    unsigned int numberDecoyProteins;
    unsigned int numberTargetProteins;
    std::string decoyPattern;
};

#endif /* PROTEINPROBESTIMATOR_H_ */
