#ifndef DESCRIPTIONOFCORRECT_H_
#define DESCRIPTIONOFCORRECT_H_
#include<string>
#include<vector>
using namespace std;
#include "PSMDescription.h"

class DescriptionOfCorrect
{
public:
  DescriptionOfCorrect();
  virtual ~DescriptionOfCorrect();
  static void fillFeaturesAllIndex(const string& peptide, double *features);
  static double isoElectricPoint(const string& peptide);
  static void setIsotopeMass(bool on) {doIsotopeMass=on;}
  void clear() {psms.clear();}
  void registerCorrect(PSMDescription* pPSM) {psms.push_back(pPSM);}
  void trainCorrect();
  void setFeatures(PSMDescription* pPSM);
  static size_t totalNumRTFeatures() {return 12 + aaAlphabet.size();}
  void print_10features();
  double estimateRT(double * features);
protected:
  void trainRetention();
  static inline double indexSum(const float *index, const string& peptide);
  static inline double indexN(const float *index, const string& peptide);
  static inline double indexC(const float *index, const string& peptide);
  static inline double indexNC(const float *index, const string& peptide);
  inline double deltadeltaMass(double dm);
  static double* fillAAFeatures(const string& pep, double *feat);
  static double* fillFeaturesIndex(const string& peptide, const float *index, double *features);
//  double kyteDolittle(string peptide);
  double avgPI,avgDM;
  size_t numRTFeat;
  
  vector<PSMDescription *> psms; 
  vector<double> rtW; 
  static float krokhin_index['Z'-'A'+1],hessa_index['Z'-'A'+1],kytedoolittle_index['Z'-'A'+1];
  static string aaAlphabet,isoAlphabet;
  static float pKiso[7]; 
  static float pKN,pKC;
  static bool doIsotopeMass;
};

#endif /*DESCRIPTIONOFCORRECT_H_*/
