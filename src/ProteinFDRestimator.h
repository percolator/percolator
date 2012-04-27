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

#define FASTA_MAXLINE 512	/* Requires FASTA file lines to be <512 characters */

template <typename T>
void FreeAll( T & t ) {
    T tmp;
    t.swap( tmp );
}

typedef struct fastafile_s {
  FILE *fp;
  char  buffer[FASTA_MAXLINE];
} FASTAFILE;

extern FASTAFILE *OpenFASTA(char *seqfile);
extern int        ReadFASTA(FASTAFILE *fp, char **ret_seq, char **ret_name, int *ret_L);
extern void       CloseFASTA(FASTAFILE *ffp);


class ProteinFDRestimator
{
  
public:
  
  ProteinFDRestimator(unsigned minpeplength = 4, unsigned minmaxx = 400, unsigned maxmass = 6000, 
		     std::string decoy_prefix = "random", double missed_cleavages = 0, unsigned nbins = 10, 
		     double targetDecoyRatio = 1.0, bool binequalDeepth = true);
  virtual ~ProteinFDRestimator();
  void parseDataBase(const char *seqfile,const char* seqfileDecoy);
  void parseDataBase(const char *seqfile);
  unsigned getNumberBins();
  unsigned getBinProteins(unsigned bin);
  unsigned countProteins(unsigned bin,const std::set<std::string> &proteins);
  double estimateFDR(const std::set<std::string> &target, const std::set<std::string> &decoy);
  void setDecoyPrefix(std::string);
  void setTargetDecoyRatio(double ratio);
  
private:
  
  double calculatePepMAss(std::string pepsequence, double charge = 2);
  unsigned calculateProtLength(std::string protsequence);
  void initMassMap(bool useAvgMass = false);
  void binProteinsEqualDeepth();
  void binProteinsEqualWidth();
  /**proteins with same sequence only the alphabetical ordered first keeps the sequence
   * the rest of the sequences are set to null. This will keep only 1 protein when there is a degenerated peptide */
  void correctIdenticalSequences(const std::map<std::string,std::string> &targetProteins,
						   const std::map<std::string,std::string> &decoyProteins);
  /**group proteins according to genes in order to estimate their lenght, proteins of the same gene group which has a tryptic peptide that has
   already been counted wont count that already counted tryptic peptide to estimate its lenght **/
  void groupProteinsGene();
  double estimatePi0HG(unsigned N,unsigned TP,unsigned FP);
  
  unsigned nbins;
  unsigned minpeplength;
  unsigned minmass;
  unsigned maxmass;
  unsigned missed_cleavages;
  double targetDecoyRatio;
  bool binequalDeepth;
  std::string decoy_prefix;
  std::map<unsigned,std::set<std::string> > binnedProteins;
  std::multimap<double,std::string> groupedProteins;
  std::map<char, double> massMap_;
  std::vector<double> lenghts; 
};
#endif /* PROTEINFDRESTIMATOR_H_ */