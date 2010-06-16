#ifndef SQTREADER_H
#define SQTREADER_H

#include <string>
#include "percolator_in.hxx"
#include "FragSpectrumScanDatabase.h"

namespace SqtReader {

  //  enum target_decoy_type { target, decoy  };

  enum parseType { justSearchMaxMinCharge, fullParsing };

  //void  TranslateSqtFileToXML(const std::string fn, const int label, ::percolatorInNs::target_decoys::target_decoy_sequence & tds , const std::string & wild, const bool match);

  //  int getScan( unsigned int scanNr, ::percolatorInNs::experiment::fragSpectrumScan_sequence  & fsss, double masscharge );


  void  translateSqtFileToXML(const std::string fn,::percolatorInNs::featureDescriptions & fds, ::percolatorInNs::experiment::fragSpectrumScan_sequence  & fsss, std::string & wild, bool is_decoy, bool calcQuadraticFeatures, bool calcAAFrequencies, bool calcPTMs, int * maxCharge,  int * minCharge, parseType t, FragSpectrumScanDatabase & database  );

  void readSQT(const std::string fn,::percolatorInNs::featureDescriptions & fds, ::percolatorInNs::experiment::fragSpectrumScan_sequence  & fsss, std::string & wild,  bool is_decoy,   bool calcQuadraticFeatures,    bool calcAAFrequencies ,     bool calcPTMs , int * maxCharge,  int * minCharge, parseType t, FragSpectrumScanDatabase & database  );

void addFeatureDescriptions( percolatorInNs::featureDescriptions & fe_des, int minC, int maxC, bool doEnzyme,
                                  bool calcPTMs, bool doPNGaseF,
			     const std::string& aaAlphabet,
			     bool calcQuadratic);

 void  readSectionS( std::string record ,  ::percolatorInNs::experiment::fragSpectrumScan_sequence & fsss, std::set<int> & theMs, bool is_decoy, bool calcPTMs, bool pngasef, bool calcAAFrequencies, int minCharge, int maxCharge, std::string psmId, FragSpectrumScanDatabase & database  );

 void readPSM( bool is_decoy, const std::string &in  ,  int match, bool calcPTMs, bool pngasef, bool calcAAFrequencies ,  ::percolatorInNs::experiment::fragSpectrumScan_sequence  & fsss,   int minCharge, int maxCharge,  std::string psmId, FragSpectrumScanDatabase & database  );



 void push_backFeatureDescription( percolatorInNs::featureDescriptions::featureDescription_sequence  & fd_sequence, const char *);
}

#endif



