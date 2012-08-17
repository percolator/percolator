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
#include "FeatureNames.h"
#include "percolator_in.hxx"
#include "parseoptions.h"
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp> 
#include <boost/assign.hpp>
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

#define BUFFER_LEN 10240

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
  
  Reader(ParseOptions *po);
  virtual ~Reader();
  
  void translateFileToXML(const std::string &fn,bool isDecoy,
			  unsigned int lineNumber_par,bool isMeta = false);

  std::string getRidOfUnprintables(const std::string &inpString);
  
  virtual void read(const std::string &fn,bool is_decoy,
		    boost::shared_ptr<FragSpectrumScanDatabase> database) = 0;
      
  virtual bool checkValidity(const std::string &file) = 0;
  
  virtual bool checkIsMeta(const std::string &file) = 0;
  
  virtual void getMaxMinCharge(const std::string &fn, bool isDecoy) = 0;
  
  virtual void addFeatureDescriptions(bool doEnzyme) = 0;
  
  void readRetentionTime(const std::string &filename);
	
  void storeRetentionTime(boost::shared_ptr<FragSpectrumScanDatabase> database);
  
  void push_backFeatureDescription(const char *str);

  void computeAAFrequencies(const string& pep,
			    percolatorInNs::features::feature_sequence & f_seq);
  
  double calculatePepMAss(const std::string &pepsequence,double charge = 2);

  void initMassMap(bool useAvgMass);
  
  void init();
  
  void print(ofstream &xmlOutputStream);
  
  unsigned int peptideLength(const string& pep);
  
  unsigned int cntPTMs(const string& pep);
  
  double isPngasef(const string& peptide, bool isDecoy );
  
  void readProteins(const std::string &filenameTarget, const std::string &fileNamedecoy);
  
  void parseDataBase(const char* seqfile, bool isDecoy,bool isCombined, 
		     unsigned &proteins_counter);
  
  unsigned calculateProtLengthTrypsin(const std::string &protsequence,
				      std::set<std::string> &peptides,double &totalMass);
  
  void read_from_fasta(istream &buffer, std::string &name , std::string &seq); 
  
  double massDiff(double observedMass, double calculatedMass,unsigned int charge);
  
  bool checkPeptideFlanks(const std::string &pep);
  
private:
  
   std::vector<char*> tmpDirs;
   std::vector<std::string> tmpFNs;

protected:
  
   static const std::string aaAlphabet;
   static const std::string ambiguousAA;
   static const std::string modifiedAA;
   static const std::string additionalAA;
   static const std::map<unsigned,double> ptmMass;
   
   std::vector<boost::shared_ptr<FragSpectrumScanDatabase> > databases;
   //NOTE I should make these two guys pointers
   ::percolatorInNs::experiment::fragSpectrumScan_sequence fss;
   ::percolatorInNs::featureDescriptions f_seq;
   int maxCharge;
   int minCharge;
   ParseOptions *po;
   std::map<char, double> massMap_;
   std::map<int, vector<double> > scan2rt;
   std::vector<Protein*> proteins;
};

#endif