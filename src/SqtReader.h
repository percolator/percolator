#ifndef SQTREADER_H
#define SQTREADER_H

#include <string>
#include "percolator-xml.hxx"


namespace SqtReader {



  //  enum target_decoy_type { target, decoy  };

  enum parseType { justSearchMaxMinCharge, fullParsing };

  //void  TranslateSqtFileToXML(const std::string fn, const int label, ::percolatorInNs::target_decoys::target_decoy_sequence & tds , const std::string & wild, const bool match);

  int getScan( unsigned int scanNr, ::percolatorInNs::experiment::frag_spectrum_scan_sequence  & fsss, double masscharge );


  void  translateSqtFileToXML(const std::string fn,::percolatorInNs::feature_descriptions & fds, ::percolatorInNs::experiment::frag_spectrum_scan_sequence  & fsss, std::string & wild,  ::percolatorInNs::type target_decoy_type, bool calcQuadraticFeatures, bool calcAAFrequencies, bool calcPTMs, int * maxCharge,  int * minCharge, parseType t );

  void readSQT(const std::string fn,::percolatorInNs::feature_descriptions & fds, ::percolatorInNs::experiment::frag_spectrum_scan_sequence  & fsss, std::string & wild, ::percolatorInNs::type target_decoy_type,   bool calcQuadraticFeatures,    bool calcAAFrequencies ,     bool calcPTMs , int * maxCharge,  int * minCharge, parseType t );

void addFeatureDescriptions( percolatorInNs::feature_descriptions & fe_des, int minC, int maxC, bool doEnzyme,
                                  bool calcPTMs, bool doPNGaseF,
			     const std::string& aaAlphabet,
			     bool calcQuadratic);

 void  readSectionS( std::string record ,  ::percolatorInNs::experiment::frag_spectrum_scan_sequence & fsss, std::set<int> & theMs, percolatorInNs::type & psmType, bool calcPTMs, bool pngasef, bool calcAAFrequencies, int minCharge, int maxCharge, std::string psmId );

 void readPSM(   percolatorInNs::type & psmType, const std::string &in  ,  int match, bool calcPTMs, bool pngasef, bool calcAAFrequencies ,  ::percolatorInNs::experiment::frag_spectrum_scan_sequence  & fsss,   int minCharge, int maxCharge,  std::string psmId );


 double isPngasef(const std::string& peptide, percolatorInNs::type & psmType );
}

#endif



