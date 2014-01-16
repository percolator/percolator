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
 
/*****************************************************************************
  
    Implementation of method to estimate protein FDR as described in :
    
    http://prottools.ethz.ch/muellelu/web/LukasReiter/Mayu/
 
    Mayu
    Lukas Reiter
    Manfred Claassen

    Mayu Software
    Lukas Reiter - Hengartner Laboratory
    lukas.reiter@molbio.uzh.ch
    Institute of Molecular Biology
    Winterthurerstrasse 190
    University of Zürich - Irchel
    CH-8057 Zürich
    +++++++++++++++++++++++++++++++++++++++++
    Located at:
    Institute for Molecular Systems Biology
    Aebersold Laboratory
    Wolfgang-Pauli-Str. 16
    ETH Hönggerberg, HPT C 75
    CH-8093 Zürich
    Tel: +41 44 633 39 45

****************************************************************************************/

#ifndef PROTEINFDRESTIMATOR_H_
#define PROTEINFDRESTIMATOR_H_

#include "Globals.h"
#include <stdio.h>
#include <functional>
#include <numeric>
#include <iterator>
#include <vector>
#include <map>
#include <set>
#include <assert.h>
#include <math.h>
#include <cmath>


template <typename T>
void FreeAll( T & t ) {
    T tmp;
    t.swap( tmp );
}

class ProteinFDRestimator
{
  
public:
  
  ProteinFDRestimator(std::string decoy_prefix = "random", unsigned nbins = 10, 
		       double targetDecoyRatio = 1.0, bool binequalDeepth = true);
  virtual ~ProteinFDRestimator();

  /** return the number of proteins in bin i **/
  unsigned getBinProteins(unsigned bin);
  
  /** return the number of proteins in bin i that are in the list of proteins given **/
  unsigned countProteins(unsigned bin,const std::set<std::string> &proteins);
  
  /** estimate and return the global FDR for a given set of target and decoy proteins **/
  double estimateFDR(const std::set<std::string> &target, const std::set<std::string> &decoy);

  /**This function populates the proteins, proteins with same sequence only the alphabetical ordered first keeps the sequence
   * the rest of the sequences are set to null. This will keep only 1 protein when there is a degenerated peptide */
  void correctIdenticalSequences(const std::map<std::string,std::pair<std::string,double> > &targetProteins,
				   const std::map<std::string,std::pair<std::string,double> > &decoyProteins);
  
  /** SETTERS AND GETTERS **/
  
  void setNumberBins(unsigned nbins);
  void setEqualDeepthBinning(bool equal_deepth);
  void setDecoyPrefix(std::string prefix);
  void setTargetDecoyRatio(double ratio);
  unsigned getNumberBins();
  bool getEqualDeppthBinning();
  double getTargetDecoyRatio();
  std::string getDecoyPrefix();
  
private:
  
  //NOTE this is a version of the code to work with threads
  /** estimates the FDR of bin i 
  void estimateFDRthread(unsigned i,const std::set<std::string> &target, const std::set<std::string> &decoy);**/

  /**bins proteins according to the lenght**/
  void binProteinsEqualDeepth();
  void binProteinsEqualWidth();

  /**group proteins according to genes in order to estimate their lenght, proteins of the same gene group which has a tryptic peptide that has
   already been counted wont count that already counted tryptic peptide to estimate its lenght **/
  void groupProteinsGene();
  
  /** estimate the expected value of the hypergeometric distributions for N,TP and FP **/
  double estimatePi0HG(unsigned N,unsigned TP,unsigned FP);
  
  /** variables **/
  std::string decoy_prefix;
  unsigned nbins;
  double targetDecoyRatio;
  //NOTE this is a version of the code to work with threads
  //double fptol;
  //std::set<std::string> *target;
  //std::set<std::string> *decoy;
  bool binequalDeepth;
  std::map<unsigned,std::set<std::string> > binnedProteins;
  std::multimap<double,std::string> groupedProteins;
  std::vector<double> lengths; 

};
#endif /* PROTEINFDRESTIMATOR_H_ */
