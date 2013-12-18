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

#ifndef PROTEIN_H
#define PROTEIN_H
#include <vector>
#include <string>
using namespace std;

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
    
    Protein();
    
    Protein(std::string namenew,double qnew, double qempnew, double pepnew, 
	    double pnew, bool isdecoy_new, Peptide *__peptide);
    
    ~Protein();
    
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

    
    /*inline bool operator<(const Protein& a,const Protein& b)  {
      return a.getPEP() > b.getPEP();
    }
    
    inline bool operator>(const Protein& a,const Protein& b)  {
      return a.getPEP() < b.getPEP();
    }
    
    inline bool operator!=(const Protein& a,const Protein& b) {
      return !(a.getName() == b.getName());
    }
    
    inline bool operator==(const Protein& a,const Protein& b) {
      return (a.getName() == b.getName());
    }*/
   
  private:
    
    std::string name;
    double q, qemp, pep, p, pi0;
    std::string id;
    bool isDecoy;
    std::vector<Peptide*> peptides;
    
};


#endif // PROTEIN_H
