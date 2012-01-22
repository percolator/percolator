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

class Protein {
  
  public:
    
    Protein() {
      q = qemp = pep = p = 0.0;
      isDecoy = false;
      name = "";
    }
    
    Protein(std::string namenew,double qnew, double qempnew, double pepnew, double pnew, bool isdecoy_new)
	      :name(namenew),q(qnew),qemp(qempnew),pep(pepnew),p(pnew),isDecoy(isdecoy_new)
	    {}
    
    ~Protein() {
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
    
    std::vector<std::string> getPeptides() {
      return peptides;
    }
    
    std::vector<std::string> getPeptides() const {
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
    
    void setPeptide(std::string peptide) {
      peptides.push_back(peptide);
    }
    
    void setPeptides(std::vector<std::string> peptidesnew)
    {
       peptides = std::vector<std::string>(peptidesnew);
    }

   
  private:
    
    std::string name;
    double q, qemp, pep, p, pi0;
    string id;
    bool isDecoy;
    std::vector<std::string> peptides;
    
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
    bool operator()(const std::pair<const std::string,Protein> &lhs, const std::pair<const std::string,Protein> &rhs) {
        return lhs.second.getPEP() < rhs.second.getPEP();
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
    
    inline double mymin(double a, double b) 
    {
	return a > b ? b : a;
    }

    
    struct pair_min {
      pair_min() {}
      std::pair<double,std::vector<std::string> > operator()(std::pair<double,std::vector<std::string> > & sum, 
							     std::pair<double,std::vector<std::string> > & i) {
        return pair<double,std::vector<std::string> >(mymin(sum.first,i.first), i.second);
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

    const static double default_gamma = 0.5; //0.01;
    const static double default_alpha = 0.1; //0.01;
    const static double default_beta = 0.01;
    
    ProteinProbEstimator(double alpha, double beta, double gamma, bool tiesAsOneProtein = false,
			 bool usePi0 = false, bool outputEmpirQVal = false, bool groupProteins = true, 
			 bool noseparate = false, bool noprune = false, bool dogridSearch = true);
    
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
    bool getTiesAsOneProtein();
    bool getUsePio();
    bool getOutputEmpirQval();
    bool getGroupProteins();
    bool getPruneProteins();
    bool getSeparateProteins();
    /** Return the scores whose q value is less or equal than the threshold given**/
    unsigned getQvaluesBelowLevel(double level);
    unsigned getQvaluesBelowLevelDecoy(double level);
    void setTargetandDecoysNames();
    std::map<const std::string,Protein> getProteins();
    void updateProteinProbabilities();
    void estimateQValues();
    void estimatePValues(); 
    void estimateQValuesEmp();
    double estimatePi0(const unsigned int numBoot = 100);
    double getPi0();
    double getAlpha();
    double getBeta();
    double getGamma();
    
    const static bool outputPEPs;
    
  private:
    
     /** fido extra functions to do the grid search for parameters alpha,betha and gamma **/
    double getROC_N(const std::vector<int> & fpArray, const std::vector<int> & tpArray, int N);
    pair<std::vector<double>, std::vector<double> > getEstimated_and_Empirical_FDR(std::vector<std::vector<string> > names, 
								   std::vector<double> probabilities);
    double getFDR_divergence(const std::vector<double> estFDR, const std::vector<double> empFDR, double THRESH);
    pair<std::vector<int>, std::vector<int> > getROC(std::vector<std::vector<string> > names);
    void gridSearch();
    

    int countTargets(std::vector<std::string> proteinList);
    int countDecoys(std::vector<std::string> proteinList);
    
    std::set<string> truePosSet, falsePosSet;
    GroupPowerBigraph* proteinGraph;
    /**proteins are mapped by the name, estQ and empQ are mapped by the q values**/
    std::map<const std::string,Protein> proteins;    
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
    double pi0;
    unsigned int numberDecoyProteins;
    unsigned int numberTargetProteins;
    double gamma;
    double alpha;
    double beta;
    
};

#endif /* PROTEINPROBESTIMATOR_H_ */
