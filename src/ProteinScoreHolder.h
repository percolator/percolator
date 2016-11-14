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
class ProteinScoreHolder {
 public:
  struct Peptide {
    Peptide(std::string _name, bool _isdecoy, double _p, double _pep, double _q, 
            double _score) : name(_name), isdecoy(_isdecoy), 
            p(_p), pep(_pep), q(_q), score(_score) {}
    std::string name;
    bool isdecoy;
    double p, pep, q, score;
  };
  
  ProteinScoreHolder() : name_(""), q_(0.0), qemp_(0.0), pep_(0.0), p_(0.0), score_(0.0),
              groupId_(-1), isDecoy_(false) {}
  ProteinScoreHolder(std::string name, bool isdecoy, Peptide *peptide, int groupId);
  ~ProteinScoreHolder();
  
  inline void setName(std::string name) { name_ = name; }
  inline std::string getName() const { return name_; }
  
  inline void setQ(double q) { q_ = q; }
  inline double getQ() const { return q_; }
  
  inline void setQemp(double qemp) { qemp_ = qemp; }
  inline double getQemp() const { return qemp_; }
  
  inline void setPEP(double pep) { pep_ = pep; }
  inline double getPEP() const { return pep_; }
  
  inline void setP(double p) { p_ = p; }
  inline double getP() const { return p_; }
  
  inline void setScore(double score) { score_ = score; }
  inline double getScore() const { return score_; }
  
  inline void setIsDecoy(bool isdecoy) { isDecoy_ = isdecoy; }
  inline bool getIsDecoy() const { return isDecoy_; }
  
  inline void setGroupId(int groupId) { groupId_ = groupId; }
  inline int getGroupId() const { return groupId_; }
  
  void addPeptide(std::string peptide, bool isdecoy, double p, double pep, 
                  double q, double empq) {
    peptides_.push_back(new Peptide(peptide, isdecoy, p, pep, q, empq));
  }
  void addPeptide(Peptide *peptide) {
    peptides_.push_back(peptide);
  }
  void setPeptides(std::vector<Peptide*> peptides) {
     peptides_ = std::vector<Peptide*>(peptides);
  }
  std::vector<Peptide*> getPeptides() { return peptides_; }
  std::vector<Peptide*> getPeptides() const { return peptides_; }
  
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
  std::string name_;
  double q_, qemp_, pep_, p_, score_;
  int groupId_;
  bool isDecoy_;
  std::vector<Peptide*> peptides_;
  
};


#endif // PROTEIN_H
