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
  double getAvgDeltaMass() {return avgDM;}
  double getAvgPI() {return avgPI;}
  static void calcRegressionFeature(PSMDescription &psm);
  static void fillFeaturesAllIndex(const string& peptide, double *features);
  static double isoElectricPoint(const string& peptide);
  static void setKlammer(bool on) {doKlammer=on;}
  static void setDocType(const unsigned int dt) {docFeatures=dt;}
  void clear() {psms.clear();}
  void registerCorrect(PSMDescription* pPSM) {psms.push_back(pPSM);}
  void trainCorrect();
  void setFeatures(PSMDescription* pPSM);
  static size_t totalNumRTFeatures() {return (doKlammer?64:minimumNumRTFeatures() + aaAlphabet.size());}
  static size_t minimumNumRTFeatures() {return 3*10+2;}
  static size_t numDOCFeatures() {return 4;}
  void print_10features();
  double estimateRT(double * features);
  void copyModel(svm_model* from);
  svm_model* getModel() {return model;}
  void copyDOCparameters(DescriptionOfCorrect& other) {
//    avgPI = other.avgPI; avgDM = other.avgDM; rtW = other.rtW; numRTFeat = other.numRTFeat;
    avgPI = other.avgPI; avgDM = other.avgDM; copyModel(other.getModel()); numRTFeat = other.numRTFeat;
  }

protected:
  void trainRetention(vector<PSMDescription *>& trainset, const double C, const double gamma);
  double testRetention(vector<PSMDescription *>& testset);
  static inline double indexSum(const float *index, const string& peptide);
  static inline double indexAvg(const float *index, const string& peptide);
  static inline double indexN(const float *index, const string& peptide);
  static inline double indexC(const float *index, const string& peptide);
  static inline double indexNC(const float *index, const string& peptide);
  static inline double* indexPartialSum(const float* index, const string& peptide, const size_t window, double *features);
  static double indexNearestNeigbourPos(const float* index, const string& peptide);
  static double indexNearestNeigbourNeg(const float* index, const string& peptide);
  static double* fillAAFeatures(const string& pep, double *feat);
  static double* fillFeaturesIndex(const string& peptide, const float *index, double *features);
  double avgPI,avgDM;
  double c,gamma;
  size_t numRTFeat;

  vector<PSMDescription *> psms;
//  vector<double> rtW;

  svm_model *model;
  static float krokhin_index['Z'-'A'+1],krokhin100_index['Z'-'A'+1],krokhinC2_index['Z'-'A'+1],krokhinTFA_index['Z'-'A'+1],
               hessa_index['Z'-'A'+1],kytedoolittle_index['Z'-'A'+1], aa_weights['Z'-'A'+1];
  static string aaAlphabet,isoAlphabet,ptmAlphabet;
  static unsigned int docFeatures;
  static float pKiso[7];
  static float pKN,pKC;
  static bool doKlammer;
};

#endif /*DESCRIPTIONOFCORRECT_H_*/
