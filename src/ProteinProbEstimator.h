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
    
    std::vector<ScoreHolderPeptide> getPeptide() {
      return peptide;
    }
    
    std::vector<ScoreHolderPeptide> getPeptide() const {
      return peptide;
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
    
    void setPeptide(ScoreHolderPeptide newpeptides) {
      peptide.push_back(newpeptides);; 
    }
    

   
  private:
    
    std::string name;
    double q, qemp, pep, p, pi0;
    string id;
    bool isDecoy;
    std::vector<ScoreHolderPeptide> peptide;
    
};

/** To be used to sort the array of proteins **/

    struct ProbOrder : public binary_function<Protein, Protein, bool> {
      bool
      operator()(const Protein& __x, const Protein& __y) const 
      {
	return (__x.getPEP() > __y.getPEP() );
      }
    };
  
    struct IntCmpProb {
    bool operator()(const std::pair<const std::string,Protein> &lhs, const std::pair<const std::string,Protein> &rhs) {
        return lhs.second.getPEP() < rhs.second.getPEP();
      }
    };
    
    struct IntCmp {
    bool operator()(const std::pair<const std::string,double> &lhs, const std::pair<const std::string,double> &rhs) {
        return lhs.second< rhs.second;
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
    
    inline std::pair<const std::string,double> mymin_pair(const std::pair<const std::string,double>& a, 
							   const std::pair<const std::string,double>& b) {
      return a.second > b.second ? b : a;
    }
    
    inline bool mycomp_pair(const std::pair<const std::string,double>& a, 
			    const std::pair<const std::string,double>& b) {
      return a.second < b.second;
    }
    
    struct pair_sum {
      
      std::pair<std::string,double> min;
      inline double mymin(double a, double b) {
	return a > b ? b : a;
      }
      std::pair<std::string,double> operator()(const std::pair<std::string,double> & i) {
        return std::pair<std::string,double>(i.first, mymin(min.second,i.second));
      }
    };
    
    inline void partial_sum(std::map<const std::string,double> a /*double (*op)(double,double)*/)
    {
       pair_sum pair_sum_func;
       std::for_each(a.rbegin(),a.rend(),pair_sum_func);
    }
    
    struct mul_x {
      mul_x(double x) : x(x) {}
      std::pair<const std::string,double> operator()(std::pair< const std::string,double> y) 
		    { return std::pair<const std::string,double>(y.first,y.second * x); }
      private:
	double x;
    };

class ProteinProbEstimator {
  public:
    double gamma;
    double alpha;
    double beta;
    const static double default_gamma = 0.5;
    const static double default_alpha = 0.1;
    const static double default_beta = 0.01;
    
    ProteinProbEstimator(double alpha, double beta, double gamma, bool tiesAsOneProtein = false,
			 bool usePi0 = false, bool outputEmpirQVal = false, bool groupProteins = false, 
			 bool noseparate = true, bool noprune = false);
    virtual ~ProteinProbEstimator();
    bool initialize(Scores* fullset);
    void run(bool startGridSearch = true);
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
    void setTargetandDecoysNames();
    std::map<const std::string,Protein> getProteins();
    void updateProteinProbabilities();
    void estimateQValues();
    void estimateQValuesEmp(double pi0);
    double getPi0();
    
    const static bool logScaleSearch;
    const static bool outputPEPs;
    
  private:
    
     /** fido extra functions to do the grid search for parameters alpha,betha and gamma **/
    int matchCount( const set<string> & positiveNames, const Array<string> & atThreshold );
    double getROC_N(const Array<int> & fpArray, const Array<int> & tpArray, int N);
    pair<Array<double>, Array<double> > getEstimated_and_Empirical_FDR(Array<Array<string> > names, 
								   Array<double> probabilities, 
								   const set<string> & falsePosSet, 
								   const set<string> & truePosSet);
    double getFDR_divergence(const Array<double> estFDR, const Array<double> empFDR, double THRESH);
    pair<Array<int>, Array<int> > getROC(Array<Array<string> > names, Array<double> probabilities, 
					 const set<string> & falsePosSet, const set<string> & truePosSet);
    void gridSearch();
    
    std::set<string> truePosSet, falsePosSet;
    GroupPowerBigraph* proteinGraph;
    /**assuming the uniqueness of the protein names**/
    std::map<const std::string,Protein> proteins;
    std::map<const std::string,double> estQ;
    std::map<const std::string,double> empQ;
    Scores* peptideScores;
    bool tiesAsOneProtein;
    bool usePi0;
    bool outputEmpirQVal;
    bool groupProteins;
    bool noseparate;
    bool noprune;
    double pi0;
    unsigned int numberDecoyProteins;
    unsigned int numberTargetProteins;
    
};

#endif /* PROTEINPROBESTIMATOR_H_ */
