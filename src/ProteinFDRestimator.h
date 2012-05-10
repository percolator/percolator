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
#include <boost/thread.hpp> 
#define FASTA_MAXLINE 512	/* Requires FASTA file lines to be <512 characters */

template <typename T>
void FreeAll( T & t ) {
    T tmp;
    t.swap( tmp );
}

/** external functions in C to read fasta files **/

/* Simple API for FASTA file reading
 * for Bio5495/BME537 Computational Molecular Biology
 * SRE, Sun Sep  8 05:35:11 2002 [AA2721, transatlantic]
 * CVS $Id$
 */

typedef struct fastafile_s {
  FILE *fp;
  char  buffer[FASTA_MAXLINE];
} FASTAFILE;

extern FASTAFILE *OpenFASTA(char *seqfile);
extern int        ReadFASTA(FASTAFILE *fp, char **ret_seq, char **ret_name, int *ret_L);
extern void       CloseFASTA(FASTAFILE *ffp);
/** external functions in C to read fasta files **/

class ProteinFDRestimator
{
  
public:
  
  ProteinFDRestimator(unsigned minpeplength = 4, unsigned minmaxx = 400, unsigned maxmass = 6000, 
		     std::string decoy_prefix = "random", double missed_cleavages = 0, unsigned nbins = 10, 
		     double targetDecoyRatio = 1.0, bool binequalDeepth = true, unsigned maxSeqlength = 40);
  virtual ~ProteinFDRestimator();
  /** extract the protein and their sequences from the fasta files (target and decoy) given **/
  void parseDataBase(const char *seqfile,const char* seqfileDecoy);
  /** extract the protein and their sequences from the combined(target and decoy) fasta files  given **/
  void parseDataBase(const char *seqfile);
  /** return the variable that defines the number of bins **/
  unsigned getNumberBins();
  /** return the number of proteins in bin i **/
  unsigned getBinProteins(unsigned bin);
  /** return the number of proteins in bin i that are in the list of proteins given **/
  unsigned countProteins(unsigned bin,const std::set<std::string> &proteins);
  /** estimate and return the global FDR for a given set of target and decoy proteins **/
  double estimateFDR(const std::set<std::string> &target, const std::set<std::string> &decoy);
  /** define the decoy prefix used to identify the decoy proteins **/
  void setDecoyPrefix(std::string prefix);
  /** define the target decoy ratio used to adjust for an inequal number of target and decyo proteins in the database **/
  void setTargetDecoyRatio(double ratio);
  
private:
  /** estimates the FDR of bin i **/
  void estimateFDRthread(unsigned i,const std::set<std::string> &target, const std::set<std::string> &decoy);
  /** estimates the mass of a peptide sequence **/
  double calculatePepMAss(std::string pepsequence, double charge = 2);
  /**estimate the number of tryptic digested peptides of a protein sequence**/
  unsigned calculateProtLength(std::string protsequence);
  /**initialize the hash table of amino acid masses**/
  void initMassMap(bool useAvgMass = false);
  /**bins proteins according to the lenght**/
  void binProteinsEqualDeepth();
  void binProteinsEqualWidth();
  /**proteins with same sequence only the alphabetical ordered first keeps the sequence
   * the rest of the sequences are set to null. This will keep only 1 protein when there is a degenerated peptide */
  void correctIdenticalSequences(const std::map<std::string,std::string> &targetProteins,
						   const std::map<std::string,std::string> &decoyProteins);
  /**group proteins according to genes in order to estimate their lenght, proteins of the same gene group which has a tryptic peptide that has
   already been counted wont count that already counted tryptic peptide to estimate its lenght **/
  void groupProteinsGene();
  /** estimate the expected value of the hypergeometric distributions for N,TP and FP **/
  double estimatePi0HG(unsigned N,unsigned TP,unsigned FP);
  
  /** variables **/
  unsigned nbins;
  unsigned minpeplength;
  unsigned minmass;
  unsigned maxmass;
  unsigned missed_cleavages;
  double targetDecoyRatio;
  double fptol;
  bool binequalDeepth;
  unsigned maxSeqlength;
  std::string decoy_prefix;
  std::map<unsigned,std::set<std::string> > binnedProteins;
  std::multimap<double,std::string> groupedProteins;
  std::map<char, double> massMap_;
  std::vector<double> lenghts; 
  std::set<std::string> *target;
  std::set<std::string> *decoy;
};
#endif /* PROTEINFDRESTIMATOR_H_ */