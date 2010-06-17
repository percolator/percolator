#ifndef SQTREADER_H
#define SQTREADER_H

#include <string>
#include "percolator_in.hxx"
using namespace std;
#include "FragSpectrumScanDatabase.h"
#include "Globals.h"
#include "SqtReader.h"

class ParseOptions {
public:
	  ParseOptions(): calcQuadraticFeatures(false),
	                  calcAAFrequencies(false),
	                  calcPTMs(false),
	                  isotopeMass(false),
	                  pngasef(false),
	                  reversedFeaturePattern(""){};

	  bool calcQuadraticFeatures;
	  bool calcAAFrequencies;
	  bool calcPTMs;
	  bool calcDOC;
	  bool isotopeMass;
	  int hitsPerSpectrum;
	  bool pngasef;
	  string reversedFeaturePattern;
}

class SqtReader {
public:
  enum parseType { justSearchMaxMinCharge, fullParsing };
  static void translateSqtFileToXML(const std::string fn,::percolatorInNs::featureDescriptions & fds, ::percolatorInNs::experiment::fragSpectrumScan_sequence  & fsss, bool is_decoy, ParseOptions & po, int * maxCharge,  int * minCharge, parseType t, FragSpectrumScanDatabase & database  );
  static void readSQT(const std::string fn,::percolatorInNs::featureDescriptions & fds, ::percolatorInNs::experiment::fragSpectrumScan_sequence  & fsss,  bool is_decoy, ParseOptions & po, int * maxCharge,  int * minCharge, parseType t, FragSpectrumScanDatabase & database  );
  static void addFeatureDescriptions( percolatorInNs::featureDescriptions & fe_des, int minC, int maxC, bool doEnzyme,
                                  bool calcPTMs, bool doPNGaseF,
			     const std::string& aaAlphabet,
			     bool calcQuadratic);
  static void readSectionS( std::string record ,  ::percolatorInNs::experiment::fragSpectrumScan_sequence & fsss, std::set<int> & theMs, bool is_decoy, ParseOptions & po, int minCharge, int maxCharge, std::string psmId, FragSpectrumScanDatabase & database  );
  static void readPSM( bool is_decoy, const std::string &in  ,  int match, ParseOptions & po ,  ::percolatorInNs::experiment::fragSpectrumScan_sequence  & fsss,   int minCharge, int maxCharge,  std::string psmId, FragSpectrumScanDatabase & database  );
  static void push_backFeatureDescription( percolatorInNs::featureDescriptions::featureDescription_sequence  & fd_sequence, const char *);
  static void computeAAFrequencies(const string& pep,   percolatorInNs::features::feature_sequence & f_seq );


}

#endif



