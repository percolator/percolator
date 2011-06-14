#ifndef SQTREADER_H
#define SQTREADER_H

#include <string>
#include <set>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <vector>
#include "FragSpectrumScanDatabase.h"
#include "Globals.h"
#include "SqtReader.h"
#include "Enzyme.h"
#include "MassHandler.h"
#include "DataSet.h"
#include "FeatureNames.h"
#include "percolator_in.hxx"
using namespace std;

class ParseOptions {
  public:
    ParseOptions(): calcQuadraticFeatures(false),
    calcAAFrequencies(false),
    calcPTMs(false),
    isotopeMass(false),
    pngasef(false),
    calcDOC(false),
    hitsPerSpectrum(1),
    reversedFeaturePattern(""){};
    bool calcQuadraticFeatures;
    bool calcAAFrequencies;
    bool calcPTMs;
    bool calcDOC;
    bool isotopeMass;
    int hitsPerSpectrum;
    bool pngasef;
    std::string reversedFeaturePattern;
    std::map<char,int> ptmScheme;
  private:
    string getRidOfUnprintables(string inpString);
};

namespace SqtReader {
enum parseType { justSearchMaxMinCharge, fullParsing };


void translateSqtFileToXML(const std::string fn,
    ::percolatorInNs::featureDescriptions & fds,
     ::percolatorInNs::experiment::fragSpectrumScan_sequence& fsss,
      bool isDecoy, const ParseOptions& po, int* maxCharge, int* minCharge,
      parseType pType, vector<FragSpectrumScanDatabase*>& databases,
      unsigned int lineNumber, std::vector<char*>& tokyoCabinetDirs,
      std::vector<std::string>& tokyoCabinetTmpFNs);
void readSQT(const std::string fn,
    ::percolatorInNs::featureDescriptions& fds,
     ::percolatorInNs::experiment::fragSpectrumScan_sequence& fsss,
      bool is_decoy, const ParseOptions& po, int* maxCharge, int* minCharge,
      parseType t, FragSpectrumScanDatabase* database);
void addFeatureDescriptions(percolatorInNs::featureDescriptions & fe_des,
    int minC, int maxC, bool doEnzyme, bool calcPTMs, bool doPNGaseF,
    const std::string& aaAlphabet, bool calcQuadratic);
void readSectionS( std::string record,
    ::percolatorInNs::experiment::fragSpectrumScan_sequence& fsss,
     std::set<int> & theMs, bool is_decoy, const ParseOptions& po,
     int minCharge, int maxCharge, std::string psmId,
     FragSpectrumScanDatabase* database);
void readPSM(bool is_decoy, const std::string &in, int match,
    const ParseOptions& po,
    ::percolatorInNs::experiment::fragSpectrumScan_sequence & fsss,
     int minCharge, int maxCharge, std::string psmId,
     FragSpectrumScanDatabase* database );
void push_backFeatureDescription(
    percolatorInNs::featureDescriptions::featureDescription_sequence&,
    const char*);
void computeAAFrequencies(const string& pep,
    percolatorInNs::features::feature_sequence & f_seq);
}

#endif



