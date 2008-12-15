#ifndef DESCRIPTIONOFCORRECT_H_
#define DESCRIPTIONOFCORRECT_H_
#include<string>
#include<iostream>
#include<vector>
using namespace std;
#include "PSMDescription.h"

struct svm_model;

class DescriptionOfCorrect
{
public:
  DescriptionOfCorrect();
  virtual ~DescriptionOfCorrect();
  static void fillFeaturesAllIndex(const string& peptide, double *features);
  static double isoElectricPoint(const string& peptide);
  static void setIsotopeMass(bool on) {doIsotopeMass=on;}
  static void setKlammer(bool on) {doKlammer=on;}
  void clear() {psms.clear();}
  void registerCorrect(PSMDescription* pPSM) {psms.push_back(pPSM);}
  void trainCorrect();
  void setFeatures(PSMDescription* pPSM);
  static size_t totalNumRTFeatures() {return (doKlammer?64:minimumNumRTFeatures() + aaAlphabet.size());}
  static size_t minimumNumRTFeatures() {return 3*9+2;}
  void print_10features();
  double estimateRT(double * features);
  void copyModel(svm_model* from);
  svm_model* getModel() {return model;}
  void copyDOCparameters(DescriptionOfCorrect& other) {
//    avgPI = other.avgPI; avgDM = other.avgDM; rtW = other.rtW; numRTFeat = other.numRTFeat;
    avgPI = other.avgPI; avgDM = other.avgDM; copyModel(other.getModel()); numRTFeat = other.numRTFeat;
  }

protected:
  void trainRetention();
  static inline double indexSum(const float *index, const string& peptide);
  static inline double indexAvg(const float *index, const string& peptide);
  static inline double indexN(const float *index, const string& peptide);
  static inline double indexC(const float *index, const string& peptide);
  static inline double indexNC(const float *index, const string& peptide);
  static inline double* indexPartialSum(const float* index, const string& peptide, const size_t window, double *features);
  static double indexNearestNeigbour(const float* index, const string& peptide);
  inline double deltadeltaMass(double dm);
  static double* fillAAFeatures(const string& pep, double *feat);
  static double* fillFeaturesIndex(const string& peptide, const float *index, double *features);
  double avgPI,avgDM;
  size_t numRTFeat;

  vector<PSMDescription *> psms; 
//  vector<double> rtW; 
  
  svm_model *model;
  static float krokhin_index['Z'-'A'+1],hessa_index['Z'-'A'+1],kytedoolittle_index['Z'-'A'+1], aa_weights['Z'-'A'+1];
  static string aaAlphabet,isoAlphabet,ptmAlphabet;
  static float pKiso[7]; 
  static float pKN,pKC;
  static bool doIsotopeMass,doKlammer;
};

#endif /*DESCRIPTIONOFCORRECT_H_*/
