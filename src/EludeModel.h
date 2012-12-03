/*******************************************************************************
 Copyright 2006-2012 Lukas Käll <lukas.kall@scilifelab.se>

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

 *****************************************************************************/

#ifndef RTMODEL_H_
#define RTMODEL_H_

#include <vector>
#include <string>
#include "PSMDescription.h"
#include "Normalizer.h"
#include "svm.h"

using namespace std;

// number of feature groups
#define NO_FEATURE_GROUPS 14
// types of grid
typedef enum {
  NO_GRID, NORMAL_GRID, FINE_GRID
} GridType;
// types of evaluation
typedef enum {
  SIMPLE_EVAL, K_FOLD_CV
} EvaluationType;
// parameter values for a grid
struct gridInformation {
    vector<double> gridGamma;
    vector<double> gridC;
    vector<double> gridEpsilon;
};

class RTModel {
  public:
    RTModel();
    ~RTModel();
    // functions to calculate retention features
    static double* amphipathicityHelix(const float* index,
                                       const string& peptide,
                                       double* features);
    static double indexSum(const float* index, const string& peptide);
    static double indexAvg(const float* index, const string& peptide);
    static double indexNearestNeigbourPos(const float* index,
                                          const string& peptide);
    static double indexNearestNeigbourNeg(const float* index,
                                          const string& peptide);
    static inline double indexN(const float* index, const string& peptide);
    static inline double indexC(const float* index, const string& peptide);
    static inline double
        indexNC(const float* index, const string& peptide);
    static double* indexPartialSum(const float* index,
                                   const string& peptide,
                                   const size_t win, double* features);
    static double* fillAAFeatures(const string& alphabet,
                                  const string& pep, double* feat);
    static double* fillPTMFeatures(const string& pep, double* feat);
    static double bulkinessSum(const string& peptide);
    static double* hydrophobicMoment(const float* index,
                                     const string& peptide,
                                     const double angle, const int window,
                                     double* features);
    static void fillFeaturesAllIndex(const string& pep, double* features);
    static double* fillFeaturesIndex(const string& peptide,
                                     const float* index, double* features);
    void calcRetentionFeatures(PSMDescription& psm);
    void calcRetentionFeatures(vector<PSMDescription> &psms);
    // calculate the difference in hydrophobicity between 2 peptides
    static double calcDiffHydrophobicities(const string& parent,
                                           const string& child);
    // train a SVR
    void trainSVM(vector<PSMDescription> & psms);
    void trainRetention(vector<PSMDescription>& trainset);
    void trainRetention(vector<PSMDescription>& trainset, const double C,
                        const double gamma, const double epsilon,
                        int noPsms);
    bool isModelNull() {
      if (model == NULL) {
        return true;
      } else {
        return false;
      }
    }
    // k-fold evaluation
    double computeKfoldCV(const vector<PSMDescription> & psms,
                          const double gamma, const double epsilon,
                          const double c);
    // use simple evaluation
    double computeSimpleEvaluation(const vector<PSMDescription> & psms,
                                   const double gamma,
                                   const double epsilon, const double c);
    // estima rt using a trained model
    double testRetention(vector<PSMDescription>& testset);
    double estimateRT(double* features);
    // load, save, copy and destroy the svr model
    void loadSVRModel(string modelFile, Normalizer* theNormalizer);
    void saveSVRModel(string modelFile, Normalizer* theNormalizer);
    void copyModel(svm_model* from);
    void destroyModel() {
      svm_destroy_model(model);
    }
    // get functions
    svm_model* getModel() {
      return model;
    }
    size_t getRTFeat() {
      return noFeaturesToCalc;
    }
    int getNoFeaturesToCalc() {
      return noFeaturesToCalc;
    }
    int getSelectFeatures() {
      return selected_features;
    }
    int getSelect(int sel_features, int max, size_t* finalNumFeatures);
    string getGridType();
    string getEvaluationType();
    static double getNoPtms(string pep);
    static size_t minimumNumRTFeatures() {
      return 3 * 17 + 1 + 1 + 1 + 2;
    }
    static size_t totalNumRTFeatures();
    // set functions
    void setNumRtFeat(const size_t nRtFeat) {
      noFeaturesToCalc = nRtFeat;
    }
    static void setDoKlammer(const bool switchKlammer);
    void setSelectFeatures(const int sf);
    void setCalibrationFile(const string calFile) {
      calibrationFile = calFile;
    }
    void setK(const int newk) {
      k = newk;
    }
    void setEvaluationType(const EvaluationType e) {
      eType = e;
    }
    void setGridType(const GridType& g);
    // print functions
    void printFeaturesInUse(ostringstream& oss);
    static string aaAlphabet, isoAlphabet;
    // print the inhouse index to the given file
    void printInhouseIndex(string& filename);
    // EXPERIMENTAL
    void trainIndexSVRNoCCalibration(vector<PSMDescription> & psms,
                                     const double C);
    void computeHydrophobicityIndex(vector<PSMDescription> & psms);
    void trainIndexRetention(vector<PSMDescription>& trainset,
                             const double C, const double epsilon);
    void copyIndexModel(svm_model* from);
    double computeKfoldCVIndex(const vector<PSMDescription> & psms,
                               const double epsilon, const double c);
    void trainIndexSVR(vector<PSMDescription> & psms);
    double estimateIndexRT(double* features);
    double testIndexRetention(vector<PSMDescription>& testset);
    void printInhouseIndex();
    bool calculateIndex();
    // EXPERIMENTAL 2
    static double noHydrophobicAA(const string& peptide);
    static double noConsecHydrophobicAA(const string& peptide);
    static double noPolarAA(const string& peptide);
    static double noConsecPolarAA(const string& peptide);
    /*static double noSmallAA(const string& peptide);
     static double noConsecAliphatic(const string& peptide);
     static double noBBranchedAA(const string& peptide);
     static double noConsecRepeats(const string& peptides, const char& letter);*/
    double* fillHydrophobicFeatures(const string& peptide,
                                    double* features);
    double* fillPolarFeatures(const string& peptide, double* features);
    static double indexSumSquaredDiff(const float* index,
                                      const string& peptide);
    void printPsms(string& filename, vector<PSMDescription> & psms);
    void setinhouseIndexAlphabet(const string& alph) {
      inhouseIndexAlphabet = alph;
    }

  protected:
    // EXPERIMENTAL
    float our_index['Z' - 'A' + 1];
    static float Luna_120_index['Z' - 'A' + 1];
    static float TFA_index['Z' - 'A' + 1];
    svm_model* index_model;
    double c_index, eps_index;
    static string inhouseIndexAlphabet;

    // number of features to be calculated (depends on the selected features)
    size_t noFeaturesToCalc;
    // svr model
    svm_model* model;
    // parameters for the SVR
    double c, gamma, epsilon;
    // type of grid
    GridType gType;
    // type of evaluation for parameter calibration
    EvaluationType eType;
    // k for cross validation
    size_t k;
    // file to save the steps of the calibration
    string calibrationFile;
    // save calibration data?
    bool saveCalibration;
    struct gridInformation grids;
    // features used (bit-wise operation are used to decode, for ex 2 means krokhin100_index)
    int selected_features;
    // number of points and step for the fine grid
    size_t noPointsFineGrid;
    double stepFineGrid;
    // symbols for post-translational modifications
    static string ptmAlphabet;
    // true if klammer features are to be used
    static bool doKlammer;
    // indices, bulkiness
    static float krokhin_index['Z' - 'A' + 1], krokhin100_index['Z' - 'A'
        + 1], krokhinC2_index['Z' - 'A' + 1], krokhinTFA_index['Z' - 'A'
        + 1], hessa_index['Z' - 'A' + 1], kytedoolittle_index['Z' - 'A'
        + 1], aa_weights['Z' - 'A' + 1], bulkiness['Z' - 'A' + 1],
        cornette_index['Z' - 'A' + 1], eisenberg_index['Z' - 'A' + 1],
        meek_index['Z' - 'A' + 1];
    // the groups of features to be used
    static string feature_groups[NO_FEATURE_GROUPS];
    // how many features are in each group?
    static int no_features_per_group[NO_FEATURE_GROUPS];
};

#endif /* RTMODEL_H_ */

