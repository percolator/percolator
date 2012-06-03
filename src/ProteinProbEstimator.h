/*******************************************************************************
 Copyright 2006-2009 Lukas Käll <lukas.kall@cbr.su.se>

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

#include "GroupPowerBigraph.h"
#include "Globals.h"
#include <boost/algorithm/string.hpp>
#include <functional>
#include <numeric>
#include <iterator>
#include "ProteinFDRestimator.h"
#include <vector>
#include <boost/assign/list_of.hpp>

/** data container that will store all the information of the proteins, peptides, type, qvalues,pvalues etc..**/

class Protein {
  
  public:
    
    struct Peptide{
      Peptide(std::string __name,bool __isdecoy,double __pep,double __q,double __empq)
      {
	name = __name;
	isdecoy = __isdecoy;
	pep =__pep;
	q = __q;
	empq = __empq;
      }
      double pep;
      double q;
      double empq;
      std::string name;
      bool isdecoy;
    };
    
    Protein() {
      q = qemp = pep = p = 0.0;
      isDecoy = false;
      name = "";
    }
    
    Protein(std::string namenew,double qnew, double qempnew, double pepnew, double pnew, bool isdecoy_new, Peptide *__peptide)
	      :name(namenew),q(qnew),qemp(qempnew),pep(pepnew),p(pnew),isDecoy(isdecoy_new)
	    {
	      if(__peptide)
		peptides.push_back(__peptide);
	    }
    
    ~Protein() {
      for(unsigned i = 0; i < peptides.size(); i++)
	delete peptides[i];
    }
    
    std::string getName()
    {
      return name;
    }
    
    std::string getName() const
    {
      return name;
    }
    
    double getQ() {
      return q;
    }
    
    double getQ() const {
      return q;
    }
    
    double getQemp() {
      return qemp;
    }
    
    double getQemp() const {
      return qemp;
    }
    
    double getPEP() {
      return pep;
    }
    
    double getPEP() const {
      return pep;
    }
    
    double getP() {
      return p; 
    }
    
    double getP() const {
      return p; 
    }
    
    bool getIsDecoy() {
      return isDecoy;
    }
    
    bool getIsDecoy() const {
      return isDecoy;
    }
    
    std::vector<Peptide*> getPeptides() {
      return peptides;
    }
    
    std::vector<Peptide*> getPeptides() const {
      return peptides;
    }
    
    void setName(std::string namenew)
    {
      name = namenew;
    }
    
    void setQ(double qnew) {
      q = qnew;
    }
    
    void setQemp(double qempnew) {
      qemp = qempnew;
    }
    
    void setIsDecoy(bool isdecoynew) {
      isDecoy = isdecoynew;
    }
    
    void setPEP(double pepnew) {
      pep = pepnew;
    }
    
    void setP(double pnew) {
      p = pnew;
    }
    
    void setPeptide(std::string peptide,bool isdecoy,double pep,double q,double empq) {
      peptides.push_back(new Peptide(peptide,isdecoy,pep,q,empq));
    }
    
    void setPeptide(Peptide *__peptide)
    {
      peptides.push_back(__peptide);
    }
    
    void setPeptides(std::vector<Peptide*> peptidesnew)
    {
       peptides = std::vector<Peptide*>(peptidesnew);
    }

   
  private:
    
    std::string name;
    double q, qemp, pep, p, pi0;
    string id;
    bool isDecoy;
    std::vector<Peptide*> peptides;
    
};

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

    inline bool operator<(const Protein& a,const Protein& b)  {
      return a.getPEP() > b.getPEP();
    }
    
    inline bool operator>(const Protein& a,const Protein& b)  {
      return a.getPEP() < b.getPEP();
    }
    
    inline bool operator!=(const Protein& a,const Protein& b) {
      return !boost::iequals(a.getName(),b.getName());
    }
    
    inline bool operator==(const Protein& a,const Protein& b) {
      return boost::iequals(a.getName(),b.getName());
    }
    
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

    
/** Interface with Fido **/

class ProteinProbEstimator {
  
  public:
    /** threshold used to estimate the protein FDR **/
    const static double psmThresholdMayu = 0.05;
    /** threshold to compute the N of the roc curve function **/
    const static double thresholdRoc = 0.05;
    
    ProteinProbEstimator(double alpha = -1, double beta = -1, double gamma = -1, bool tiesAsOneProtein = false,
			 bool usePi0 = false, bool outputEmpirQVal = false, bool groupProteins = false, 
			 bool noseparate = false, bool noprune = false, bool dogridSearch = true, unsigned depth = 3,
			 double lambda = 0.15, double threshold = 0.05, unsigned rocN = 0, std::string targetDB = "", 
			 std::string decoyDB = "", std::string decoyPattern = "random", bool mayufdr = false, 
			 bool conservative = false, bool outputDecoys = false, bool tabDelimitedOut = false, 
			 std::string proteinFN = "");
    
    virtual ~ProteinProbEstimator();
    
    /** reads the proteins from the set of scored peptides from percolator **/
    bool initialize(Scores* fullset);
    /** initialize the estimation of the protein probabilities and the grid search **/
    void run();
    /** write the list of proteins to the output file **/
    void writeOutputToXML(string xmlOutputFN);

    /** Return the number of proteins whose q value is less or equal than the threshold given**/
    unsigned getQvaluesBelowLevel(double level);
    unsigned getQvaluesBelowLevelDecoy(double level);
   
    /** populate the list of proteins**/
    void setTargetandDecoysNames();
    /** return the data structure for the proteins **/
    std::map<const std::string,Protein*> getProteins();
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
    /** print a tab delimited list of proteins probabilities in a file or stdout**/
    void print(ostream& myout);
    
    static string printCopyright();
	
     /**setters and getters for variables **/
    void setTiesAsOneProtein(bool tiesAsOneProtein);
    void setUsePio(bool usePi0);
    void setOutputEmpirQval(bool outputEmpirQVal);
    void setGroupProteins(bool groupProteins);
    void setPruneProteins(bool noprune);
    void setSeparateProteins(bool noseparate);
    void setGridSearch(bool dogridSearch);
    void setDepth(unsigned depth);
    void setLambda(double lambda);
    void setThreshold(double threshold);
    void setROCN(double rocn);
    void setTargetDb(std::string targetDB);
    void setDecoyDb(std::string decoyDB);
    void setMayusFDR(bool mayufdr);
    void setFDR(double fdr);
    void setOutputDecoys(bool outputDecoys);
    void setTabDelimitedOutput(bool tabDelimitedOut);
    void setProteinFN(std::string proteinFN);
    bool getTiesAsOneProtein();
    bool getUsePio();
    bool getOutputEmpirQval();
    bool getGroupProteins();
    bool getPruneProteins();
    bool getSeparateProteins();
    bool getMayuFdr();
    bool getDepth();
    bool getGridSearch();
    bool getOutputDecoys();
    bool getTabDelimitedOutput();
    std::string getProteinFN();
    std::string getDecoyPatter();
    std::string getDecoyDB();
    std::string getTargetDB();
    unsigned getROCN();
    double getThreshold();
    double getLambda();
    double getFDR();
    double getPi0();
    double getAlpha();
    double getBeta();
    double getGamma();
    
  private:
    
     /** fido extra functions to do the grid search for parameters alpha,betha and gamma **/
    double getROC_N(const std::vector<unsigned> &fpArray, const std::vector<unsigned> &tpArray, int N);
    void getEstimated_and_Empirical_FDR(const std::vector<std::vector<string> > &names,
					   const std::vector<double> &probabilities,
					   std::vector<double> &empq,
					   std::vector<double> &estq);
    double getFDR_divergence(const std::vector<double> &estFDR, const std::vector<double> &empFDR, double THRESH);
    void getROC(const std::vector<std::vector<string> > &names,std::vector<unsigned> &numberFP,std::vector<unsigned> &numberTP);
    void gridSearch(double alpha = -1, double gamma = -1, double  beta = -1);
    
    /** functions to count number of target and decoy proteins **/
    unsigned countTargets(const std::vector<std::string> &proteinList);
    unsigned countDecoys(const std::vector<std::string> &proteinList);
    bool isTarget(const std::string& proteinName);
    bool isDecoy(const std::string& proteinName);
    
    /** function that extracts a list of proteins from the peptides that have a qvalue lower than psmThresholdMayu
     * this function is used to estimate the protein FDR**/
    void getTPandPFfromPeptides(double threshold, std::set<std::string> &numberTP, std::set<std::string> &numberFP);
   
    /** variables **/
    std::set<string> truePosSet, falsePosSet;
    GroupPowerBigraph* proteinGraph;
    ProteinFDRestimator *fastReader;
    std::map<const std::string,Protein*> proteins;    
    std::multimap<double,std::vector<std::string> > pepProteins;
    std::vector<double> qvalues;
    std::vector<double> qvaluesEmp;
    std::vector<double> pvalues;
    Scores* peptideScores;
    bool tiesAsOneProtein;
    bool usePi0;
    bool outputEmpirQVal;
    bool groupProteins;
    bool noseparate;
    bool noprune;
    bool dogridSearch;
    bool mayufdr;
    bool updateRocN;
    bool conservative;
    bool tabDelimitedOut;
    bool outputDecoys;
    double pi0;
    double fdr;
    unsigned int numberDecoyProteins;
    unsigned int numberTargetProteins;
    unsigned int depth;
    double gamma;
    double alpha;
    double beta;
    double lambda;
    double threshold;
    unsigned rocN;
    std::string targetDB;
    std::string decoyDB;
    std::string decoyPattern;
    std::string proteinFN;
    
};

#endif /* PROTEINPROBESTIMATOR_H_ */
