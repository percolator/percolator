#include "Reader.h"
#include <sys/types.h>
#include <sys/stat.h>

const std::string Reader::aaAlphabet("ACDEFGHIKLMNPQRSTVWY");
const std::string Reader::ambiguousAA("BZJX");
const std::string Reader::modifiedAA("#@*");

Reader::Reader(ParseOptions __po)
:po(__po)
{  
  tmpDirs = std::vector<char*>();
  tmpFNs = std::vector<std::string>();
  maxCharge = -1;
  minCharge = 10000;
  initMassMap(po.monoisotopic);
}

Reader::~Reader()
{
  //NOTE do this with boost?
  //NOTE remove files as well
  for(int i=0; i<tmpDirs.size(); i++)
    rmdir(tmpDirs[i]);
}



void Reader::push_backFeatureDescription( percolatorInNs::featureDescriptions::featureDescription_sequence  & fd_sequence , const char * str) {

  std::auto_ptr< ::percolatorInNs::featureDescription > f_p( new ::percolatorInNs::featureDescription(str));
  assert(f_p.get());
  fd_sequence.push_back(f_p);
  return;
}


/**
 * remove non ASCII characters from a string
 */
string Reader::getRidOfUnprintables(string inpString) {
  string outputs = "";
  for (int jj = 0; jj < inpString.size(); jj++) {
    signed char ch = inpString[jj];
    if (((int)ch) >= 32 && ((int)ch) <= 128) {
      outputs += ch;
    }
  }
  return outputs;
}

void Reader::computeAAFrequencies(const string& pep,   percolatorInNs::features::feature_sequence & f_seq ) {
  // Overall amino acid composition features
  assert(pep.size() >= 5);
  string::size_type aaSize = aaAlphabet.size();

  std::vector< double > doubleV;
  for ( int m = 0  ; m < aaSize ; m++ )  {
    doubleV.push_back(0.0);
  }
  int len = 0;
  for (string::const_iterator it = pep.begin() + 2; it != pep.end() - 2; it++) {
    string::size_type pos = aaAlphabet.find(*it);
    if (pos != string::npos) doubleV[pos]++;
    len++;
  }
  assert(len>0);
  for ( int m = 0  ; m < aaSize ; m++ )  {
    doubleV[m] /= len;
  }
  std::copy(doubleV.begin(), doubleV.end(), std::back_inserter(f_seq));
}

double Reader::calculatePepMAss(std::string pepsequence,double charge)
{
  double mass  =  0.0;
  if (pepsequence.length () > po.peptidelength) {
    
    for(unsigned i=0; i<pepsequence.length();i++)
    {
      if(isalpha(pepsequence[i])){
        mass += massMap_[pepsequence[i]];
      }
    }
    
    mass = (mass + massMap_['o'] + (charge * massMap_['h'])); 
  }
  return mass; 
}

void Reader::initMassMap(bool useAvgMass)
{
  if (useAvgMass) /*avg masses*/
    {
      massMap_['h']=  1.00794;  
      massMap_['o']= 15.9994;   
      massMap_['G']= 57.05192;
      massMap_['A']= 71.07880;
      massMap_['S']= 87.07820;
      massMap_['P']= 97.11668;
      massMap_['V']= 99.13256;
      massMap_['T']=101.10508;
      massMap_['C']=103.13880;
      massMap_['L']=113.15944;
      massMap_['I']=113.15944;
      massMap_['X']=113.15944;
      massMap_['N']=114.10384;
      massMap_['O']=114.14720;
      massMap_['B']=114.59622;
      massMap_['D']=115.08860;
      massMap_['Q']=128.13072;
      massMap_['K']=128.17408;
      massMap_['Z']=128.62310;
      massMap_['E']=129.11548;
      massMap_['M']=131.19256;
      massMap_['H']=137.14108;
      massMap_['F']=147.17656;
      massMap_['R']=156.18748;
      massMap_['Y']=163.17596;
      massMap_['W']=186.21320;
    }
  else /* monoisotopic masses */
    {
      massMap_['h']=  1.0078250;
      massMap_['o']= 15.9949146;
      massMap_['A']= 71.0371136;
      massMap_['C']=103.0091854;
      massMap_['D']=115.0269428;
      massMap_['E']=129.0425928;
      massMap_['F']=147.0684136;
      massMap_['G']= 57.0214636;
      massMap_['H']=137.0589116;
      massMap_['I']=113.0840636;
      massMap_['K']=128.0949626;
      massMap_['L']=113.0840636;
      massMap_['M']=131.0404854;
      massMap_['N']=114.0429272;
      massMap_['P']= 97.0527636;
      massMap_['Q']=128.0585772;
      massMap_['R']=156.1011106;
      massMap_['S']= 87.0320282;
      massMap_['T']=101.0476782;
      massMap_['U']=149.90419;
      massMap_['V']= 99.0684136;
      massMap_['W']=186.07931;
      massMap_['Y']=163.06333;
    }
}