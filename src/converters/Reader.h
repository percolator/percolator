#ifndef READER_H
#define READER_H

#include <string>
#include <set>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <vector>
#include "FragSpectrumScanDatabase.h"
#include "Globals.h"
#include "Enzyme.h"
#include "MassHandler.h"
#include "DataSet.h"
#include "FeatureNames.h"
#include "percolator_in.hxx"
#include "parseoptions.h"
#include <boost/filesystem.hpp>

using namespace std;

class Reader
{

public:
  
  Reader();
  virtual ~Reader();
  
  enum parseType { justSearchMaxMinCharge, fullParsing };
  
  void translateSqtFileToXML(const std::string fn,
    ::percolatorInNs::featureDescriptions & fds,
     ::percolatorInNs::experiment::fragSpectrumScan_sequence& fsss,
      bool isDecoy, const ParseOptions& po, int* maxCharge, int* minCharge,
      parseType pType, vector<FragSpectrumScanDatabase*>& databases,
      unsigned int lineNumber/*, std::vector<char*>& DBtDirs,std::vector<std::string>& DBTmpFNs*/);

  string getRidOfUnprintables(string inpString);
  
  void setFile(std::string file);
  
  virtual void read(const std::string fn,
    ::percolatorInNs::featureDescriptions& fds,
     ::percolatorInNs::experiment::fragSpectrumScan_sequence& fsss,
      bool is_decoy, const ParseOptions& po, int* maxCharge, int* minCharge,
      parseType t, FragSpectrumScanDatabase* database){};

  void addFeatureDescriptions(percolatorInNs::featureDescriptions & fe_des,
    int minC, int maxC, bool doEnzyme, bool calcPTMs, bool doPNGaseF,
    const std::string& aaAlphabet, bool calcQuadratic);
  
  void push_backFeatureDescription(
    percolatorInNs::featureDescriptions::featureDescription_sequence&,
    const char*);

  void computeAAFrequencies(const string& pep,
    percolatorInNs::features::feature_sequence & f_seq);

private:
  
   std::vector<char*> tmpDirs;
   std::vector<std::string> tmpFNs;

protected:
  
   std::string aaAlphabet; 
   std::string ambiguousAA;
   std::string modifiedAA;
   
};


#endif



