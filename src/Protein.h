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

/** data container that will store all the information of the proteins, 
    peptides, type, qvalues, pvalues etc..**/
class Protein {
 public:
  struct Peptide {
    Peptide(std::string __name, bool __isdecoy, double __pep, double __q, 
            double __empq) : name(__name), isdecoy(__isdecoy), pep(__pep),
            q(__q), empq(__empq) {}
    std::string name;
    bool isdecoy;
    double pep, q, empq;
  };
  
  Protein() : name(""), q(0.0), qemp(0.0), pep(0.0), p(0.0), groupId_(-1), 
              isDecoy(false) {}
  Protein(std::string namenew,double qnew, double qempnew, double pepnew, 
    double pnew, bool isdecoy_new, Peptide *__peptide, int groupId);
  ~Protein();
  
  void setName(std::string namenew) { name = namenew; }
  std::string getName() const { return name; }
  
  void setQ(double qnew) { q = qnew; }
  double getQ() const { return q; }
  
  void setQemp(double qempnew) { qemp = qempnew; }
  double getQemp() const { return qemp; }
  
  void setPEP(double pepnew) { pep = pepnew; }
  double getPEP() const { return pep; }
  
  void setP(double pnew) { p = pnew; }
  double getP() const { return p; }
  
  void setIsDecoy(bool isdecoynew) { isDecoy = isdecoynew; }
  bool getIsDecoy() const { return isDecoy; }
  
  void setGroupId(int groupId) { groupId_ = groupId; }
  int getGroupId() const { return groupId_; }
  
  void setPeptide(std::string peptide,bool isdecoy,double pep,double q,double empq) {
    peptides.push_back(new Peptide(peptide,isdecoy,pep,q,empq));
  }
  void setPeptide(Peptide *__peptide) {
    peptides.push_back(__peptide);
  }
  void setPeptides(std::vector<Peptide*> peptidesnew) {
     peptides = std::vector<Peptide*>(peptidesnew);
  }
  std::vector<Peptide*> getPeptides() { return peptides; }
  std::vector<Peptide*> getPeptides() const { return peptides; }
  
  /*
  inline bool operator<(const Protein& a,const Protein& b)  {
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
  }
  */
 
 private:
  std::string name;
  double q, qemp, pep, p, pi0;
  int groupId_;
  bool isDecoy;
  std::vector<Peptide*> peptides;
  
};


#endif // PROTEIN_H
