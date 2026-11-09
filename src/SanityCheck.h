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
#ifndef SANITYCHECK_H_
#define SANITYCHECK_H_

#include <iostream>
#include <string>
#include <vector>

class Scores;
class Normalizer;

/*
 * SanityCheck sets the initial direction of the SVM, either from an input
 * weight file, a default direction given by the user or in the input file,
 * or the single most discriminative feature, and checks that the trained
 * direction is not worse than the initial one.
 */
class SanityCheck {
 public:
  void readWeights(std::istream& weightStream, std::vector<double>& w);
  void checkAndSetDefaultDir();
  int getInitDirection(std::vector<Scores>& testset,
                       std::vector<Scores>& trainset, const Normalizer* pNorm,
                       std::vector<std::vector<double> >& w, double test_fdr,
                       double initial_train_fdr);
  bool validateDirection(std::vector<Scores>& testset,
                         std::vector<Scores>& trainset,
                         std::vector<std::vector<double> >& w);
  void resetDirection(std::vector<Scores>& trainset,
                      std::vector<std::vector<double> >& w);

  void setInitWeightFN(const std::string& fn) { initWeightFN_ = fn; }
  void setInitDefaultDirName(const std::string& dirName) {
    initDefaultDirName_ = dirName;
  }
  void setOverrule(bool orl) { overRule_ = orl; }
  void addDefaultWeights(const std::vector<double>& defaultWeights) {
    defaultWeights_ = defaultWeights;
  }
  const std::vector<double>& getDefaultWeights() const {
    return defaultWeights_;
  }
  bool concatenatedSearch() const { return concatenatedSearch_; }
  void setConcatenatedSearch(const bool s) { concatenatedSearch_ = s; }

 protected:
  void getDefaultDirection(std::vector<Scores>& trainset,
                           std::vector<std::vector<double> >& w);

  int initPositives_ = 0;
  double test_fdr_ = 0.01, initial_train_fdr_ = 0.01;
  bool overRule_ = false;
  std::string initWeightFN_;
  std::string initDefaultDirName_;
  // Default Direction, 0=do not use, positive integer = feature number,
  // negative integer = lower score better
  int initDefaultDir_ = 0;
  // initial direction given in the input file
  std::vector<double> defaultWeights_;
  // input from concatenated search, i.e. not separate target and decoy searches
  bool concatenatedSearch_ = true;
};

#endif /*SANITYCHECK_H_*/
