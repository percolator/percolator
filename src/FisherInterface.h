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

#ifndef FISHERINTERFACE_H
#define FISHERINTERFACE_H

#include <cmath>
#include <functional>
#include <cfloat>
//#include <boost/math/special_functions/gamma.hpp>

#include "ProteinProbEstimator.h"
#include "PosteriorEstimator.h"
#include "FisherCaller.h"
#include "Enzyme.h"
#include "PseudoRandom.h"

/* 
* Protein inference methods:
* 1. Fisher's method for combining p-values
* 2. Product of Posterior error probabilities (PEPs)
* 3. Best peptide approach
*/
enum ProteinInferenceMethod {
  FISHER, PEPPROD, BESTPEPT
};

/*
* FisherInterface is a class that computes probabilities and statistics based
* on provided proteins from the set of scored peptides from percolator. It
* uses number of additional parameters to setup the calculations.
*
* Here are some usefull abbreviations:
* FDR - False Discovery Rate
* PSM - Peptide Spectrum Match
*
*/
class FisherInterface : public ProteinProbEstimator {

 public:
  FisherInterface(const std::string& fastaDatabase, double pvalueCutoff,
    bool reportFragmentProteins, bool reportDuplicateProteins, 
    bool trivialGrouping, double absenceRatio, 
    bool outputEmpirQval, std::string& decoyPattern);
  virtual ~FisherInterface();
  
  bool initialize(Scores* fullset);
  void run();
  void computeProbabilities(const std::string& fname = "");
  
  std::ostream& printParametersXML(std::ostream &os);
  string printCopyright();

 private:
  void pickedProteinStrategy(
    std::vector<std::pair<std::string,Protein*> >& protIdProtPairs);
  void pickedProteinStrategySubstring(
    std::vector<std::pair<std::string,Protein*> >& protIdProtPairs);
  void estimatePEPs(
    std::vector<std::pair<std::string,Protein*> >& protIdProtPairs,
    std::vector<double>& peps);
  
  /** FISHER PARAMETERS **/
  ProteinInferenceMethod protInferenceMethod_;
  std::string fastaProteinFN_;
  bool reportFragmentProteins_, reportDuplicateProteins_;
  FisherCaller fisherCaller_;
  double maxPeptidePval_;
  
};

#endif // FISHERINTERFACE_H
