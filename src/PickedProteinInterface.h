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

#ifndef PICKED_PROTEIN_INTERFACE_H
#define PICKED_PROTEIN_INTERFACE_H

#include <cmath>
#include <functional>
#include <cfloat>
//#include <boost/math/special_functions/gamma.hpp>

#include "ProteinProbEstimator.h"
#include "PosteriorEstimator.h"
#include "PickedProteinCaller.h"
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
* PickedProteinInterface is a class that computes probabilities and statistics based
* on provided proteins from the set of scored peptides from percolator. It
* uses number of additional parameters to setup the calculations.
*
* Here are some usefull abbreviations:
* FDR - False Discovery Rate
* PSM - Peptide Spectrum Match
*
*/
class PickedProteinInterface : public ProteinProbEstimator {

 public:
  PickedProteinInterface(const std::string& fastaDatabase, double pvalueCutoff,
    bool reportFragmentProteins, bool reportDuplicateProteins, 
    bool trivialGrouping, double absenceRatio, 
    bool outputEmpirQval, std::string& decoyPattern);
  virtual ~PickedProteinInterface();
  
  bool initialize(Scores& fullset);
  void run() {}
  void computeProbabilities(const std::string& fname = "");
  
  std::ostream& printParametersXML(std::ostream &os);
  string printCopyright();

 private:
  void groupProteins(Scores& peptideScores);
  
  void pickedProteinStrategy();
  bool pickedProteinCheckId(std::string& proteinId, bool isDecoy,
    std::set<std::string>& targetProts, std::set<std::string>& decoyProts);
  bool pickedProteinCheck(std::string& proteinName, bool isDecoy, 
    std::set<std::string>& targetProts, std::set<std::string>& decoyProts);
  void estimatePEPs();
  
  /** PICKED_PROTEIN PARAMETERS **/
  ProteinInferenceMethod protInferenceMethod_;
  std::string fastaProteinFN_;
  bool reportFragmentProteins_, reportDuplicateProteins_;
  PickedProteinCaller fisherCaller_;
  double maxPeptidePval_;
  
};

#endif // PICKED_PROTEININTERFACE_H
