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
#include "converters/MSToolkit/MSToolkitTypes.h"
#include <iterator>
#include "FastaProteinReader.h"
#include <vector>


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

/** To be used to sort the array of proteins **/

    struct ProbOrder : public binary_function<Protein, Protein, bool> {
      bool
      operator()(const Protein& __x, const Protein& __y) const 
      {
	return (__x.getPEP() < __y.getPEP() );
      }
    };
  
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
    
    struct IntCmp {
    bool operator()(const std::pair<const double,std::vector<std::string> > &lhs, 
		    const std::pair<const double,std::vector<std::string> > &rhs) {
        return lhs.first < rhs.first;
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
    
    inline bool mycomp_pair(const std::pair<const std::string,double>& a, 
			    const std::pair<const std::string,double>& b) {
      return a.second < b.second;
    }
    
    inline double myminfunc(double a, double b) 
    {
	return a > b ? b : a;
    }

    
    struct pair_min {
      pair_min() {}
      std::pair<double,std::vector<std::string> > operator()(std::pair<double,std::vector<std::string> > & sum, 
							     std::pair<double,std::vector<std::string> > & i) {
        return pair<double,std::vector<std::string> >(myminfunc(sum.first,i.first), i.second);
      }
    };


    struct mul_x {
      mul_x(double x) : x(x) {}
      std::pair<double,std::vector<std::string> > 
      operator()(std::pair<double,std::vector<std::string> > y) 
      { return std::pair<double,std::vector<std::string> >(y.first * x,y.second ); }
      private:
	double x;
    };

    inline std::map<const double,std::vector<std::string> >::const_iterator 
    MapSearchByValue(const std::map<double,std::vector<std::string> > & SearchMap, 
		     const std::string & SearchVal)
    {
      std::map<double,std::vector<std::string> >::const_iterator iRet = SearchMap.end();
      for (std::map<double,std::vector<std::string> >::const_iterator iTer = SearchMap.begin(); 
	   iTer != SearchMap.end(); iTer ++)
      {
        if (std::find(iTer->second.begin(),iTer->second.end(),SearchVal) != iTer->second.end())
        {
            iRet = iTer;
            break;
        }
      }
      return iRet;
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


class ProteinProbEstimator {
  
  public:

    const static double default_gamma = 0.5; 
    const static double default_alpha = 0.1; 
    const static double default_beta = 0.01;
    const static double psmThresholdMayu = 0.05;
    const static double thresholdRoc = 0.05;
    
    ProteinProbEstimator(double alpha, double beta, double gamma, bool tiesAsOneProtein = false,
			 bool usePi0 = false, bool outputEmpirQVal = false, bool groupProteins = false, 
			 bool noseparate = false, bool noprune = false, bool dogridSearch = true, unsigned deepness = 3,
			 double lambda = 0.15, double threshold = 0.05, unsigned rocN = 0, std::string targetDB = "", 
			 std::string decoyDB = "", std::string decoyPattern = "random", bool mayufdr = false, bool conservative = false);
    
    ~ProteinProbEstimator();
    
    bool initialize(Scores* fullset);
    void run();
    void writeOutputToXML(string xmlOutputFN);
    static string printCopyright();
    
    /**setters and getters for constants **/
    void setTiesAsOneProtein(bool tiesAsOneProtein);
    void setUsePio(bool usePi0);
    void setOutputEmpirQval(bool outputEmpirQVal);
    void setGroupProteins(bool groupProteins);
    void setPruneProteins(bool noprune);
    void setSeparateProteins(bool noseparate);
    void setGridSearch(bool dogridSearch);
    void setDeepness(unsigned deepness);
    void setLambda(double lambda);
    void setThreshold(double threshold);
    void setROCN(double rocn);
    void setTargetDb(std::string targetDB);
    void setDecoyDb(std::string decoyDB);
    void setMayusFDR(bool mayufdr);
    void setFDR(double fdr);
    bool getTiesAsOneProtein();
    bool getUsePio();
    bool getOutputEmpirQval();
    bool getGroupProteins();
    bool getPruneProteins();
    bool getSeparateProteins();
    bool getMayuFdr();
    bool getDeepness();
    bool getGridSearch();
    std::string getDecoyPatter();
    std::string getDecoyDB();
    std::string getTargetDB();
    unsigned getROCN();
    double getThreshold();
    double getLambda();
    double getFDR();
    /** Return the scores whose q value is less or equal than the threshold given**/
    unsigned getQvaluesBelowLevel(double level);
    unsigned getQvaluesBelowLevelDecoy(double level);
    
    void setTargetandDecoysNames();
    std::map<const std::string,Protein*> getProteins();
    void updateProteinProbabilities();
    void estimateQValues();
    void estimatePValues(); 
    void estimateQValuesEmp();
    double estimatePi0(const unsigned int numBoot = 100);
    double getPi0();
    double getAlpha();
    double getBeta();
    double getGamma();
    
  private:
    
     /** fido extra functions to do the grid search for parameters alpha,betha and gamma **/
    double getROC_N(const std::vector<int> & fpArray, const std::vector<int> & tpArray, int N);
    std::pair<std::vector<double>, std::vector<double> > getEstimated_and_Empirical_FDR(std::vector<std::vector<string> > names, 
								std::vector<double> probabilities);
    double getFDR_divergence(const std::vector<double> estFDR, const std::vector<double> empFDR, double THRESH);
    std::pair<std::vector<int>, std::vector<int> > getROC(std::vector<std::vector<string> > names);
    void gridSearch(double alpha = -1, double gamma = -1, double  beta = -1);
    int countTargets(std::vector<std::string> proteinList);
    int countDecoys(std::vector<std::string> proteinList);
    std::pair<std::set<std::string>,std::set<std::string> > getTPandPFfromPeptides(double threshold);
    double estimatePi0Bin(unsigned FP,unsigned TP);
    
    std::set<string> truePosSet, falsePosSet;
    GroupPowerBigraph* proteinGraph;
    FastaProteinReader *fastReader;
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
    double pi0;
    unsigned int numberDecoyProteins;
    unsigned int numberTargetProteins;
    unsigned int deepness;
    double gamma;
    double alpha;
    double beta;
    double lambda;
    double threshold;
    unsigned rocN;
    std::string targetDB;
    std::string decoyDB;
    std::string decoyPattern;
    
};

#endif /* PROTEINPROBESTIMATOR_H_ */
