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

#ifndef READER_H
#define READER_H

#include "Enzyme.h"
#include <iostream>
#include <fstream>
#include <numeric>
#include <map>
#include <xercesc/dom/DOM.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xsd/cxx/xml/string.hxx>
#include <sys/types.h>
#include <sys/stat.h>
#include <string>
#include <set>
#include <cmath>
#include <algorithm>
#include <vector>
#include "Globals.h"
#include "config.h"
#include "MassHandler.h"
#include "FeatureNames.h"
#include "percolator_in.hxx"
#include "parseoptions.h"
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp> 
#include "MSReader.h"
#include "Spectrum.h"
#include <assert.h>
#if defined (__WIN32__) || defined (__MINGW__) 
  #include <direct.h>
  #include <io.h>
  #include <stdio.h>
  #define  mkdir( D, M )   _mkdir( D )
  #include <fcntl.h>
  #include <errno.h>
#endif

#define FASTA_MAXLINE 512	/* Requires FASTA file lines to be <512 characters */

#if defined __LEVELDB__
  #include "FragSpectrumScanDatabaseLeveldb.h"
  typedef FragSpectrumScanDatabaseLeveldb serialize_scheme;
#elif defined __TOKYODB__ 
  #include "FragSpectrumScanDatabaseTokyodb.h"
  typedef FragSpectrumScanDatabaseTokyoDB serialize_scheme;
#else
  #include "FragSpectrumScanDatabaseBoostdb.h"
  typedef FragSpectrumScanDatabaseBoostdb serialize_scheme;
#endif
   
using namespace std;

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

class Protein
{
  
public:
  
  Protein() {
   length = 0;
   totalMass = 0.0;
   isDecoy = false;
   name = "";
   id = 0;
  }
    
  std::string name;
  std::string sequence;
  double totalMass;
  unsigned id;
  bool isDecoy;
  unsigned length;
  std::set<std::string> peptides;
};

class Reader
{

public:
  
  Reader(ParseOptions po);
  virtual ~Reader();
  
  void translateFileToXML(const std::string fn,bool isDecoy,unsigned int lineNumber_par);

  string getRidOfUnprintables(std::string inpString);
  
  virtual void read(const std::string fn,bool is_decoy,boost::shared_ptr<FragSpectrumScanDatabase> database){};
      
  virtual bool checkValidity(std::string file){};
  
  virtual void getMaxMinCharge(std::string fn){};
  
  virtual void addFeatureDescriptions(bool doEnzyme,const std::string& aaAlphabet,std::string fn){};
  
  virtual void readRetentionTime(std::string filename){};
	
  virtual void storeRetentionTime(boost::shared_ptr<FragSpectrumScanDatabase> database){};
  
  void push_backFeatureDescription(const char *str);

  void computeAAFrequencies(const string& pep,percolatorInNs::features::feature_sequence & f_seq);
  
  double calculatePepMAss(const std::string &pepsequence,double charge = 2);

  void initMassMap(bool useAvgMass);
  
  void init();
  
  void print(ofstream &xmlOutputStream);
  
  unsigned int peptideLength(const string& pep);
  
  unsigned int cntPTMs(const string& pep);
  
  double isPngasef(const string& peptide, bool isDecoy );
  
  void readProteins(std::string filenameTarget,std::string fileNamedecoy);
  
  void parseDataBase(const char* seqfile, bool isDecoy,bool isCombined, unsigned &proteins_counter);
  
  unsigned calculateProtLengthTrypsin(std::string protsequence,
				      std::set<std::string> &peptides,double &totalMass);
  
private:
  
   std::vector<char*> tmpDirs;
   std::vector<std::string> tmpFNs;

protected:
  
   static const std::string aaAlphabet;
   static const std::string ambiguousAA;
   static const std::string modifiedAA;
   std::vector<boost::shared_ptr<FragSpectrumScanDatabase> > databases;
   //NOTE as soon as I get the program working these 2 guys have to be smart pointers
   ::percolatorInNs::experiment::fragSpectrumScan_sequence fss;
   ::percolatorInNs::featureDescriptions f_seq;
   int maxCharge;
   int minCharge;
   ParseOptions po;
   std::map<char, double> massMap_;
   std::map<int, vector<double> > scan2rt;
   std::vector<Protein*> proteins;
};

#endif



