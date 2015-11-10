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

 *******************************************************************************/
#ifndef DESCRIPTIONOFCORRECT_H_
#define DESCRIPTIONOFCORRECT_H_
#include <functional>
#include <string>
#include <vector>
using namespace std;
#include "EludeModel.h"


class DescriptionOfCorrect {
  public:
    DescriptionOfCorrect();
    virtual ~DescriptionOfCorrect();
    double getAvgDeltaMass() {
      return avgDM;
    }
    double getAvgPI() {
      return avgPI;
    }
    static void calcRegressionFeature(PSMDescription& psm);
    static double isoElectricPoint(const string& peptide);
    static void setKlammer(bool on) {
      RTModel::setDoKlammer(on);
    }
    static void setDocType(const unsigned int dt) {
      docFeatures = dt;
    }
    void clear() {
      psms.clear();
    }
    void registerCorrect(PSMDescription& psm) {
      psms.push_back(psm);
    }
    void trainCorrect();
    void setFeatures(PSMDescription& psm);
    void setFeaturesNormalized(PSMDescription& psm, Normalizer* pNorm);
    //static size_t totalNumRTFeatures() {return (doKlammer?64:minimumNumRTFeatures() + 20);}
    //static size_t minimumNumRTFeatures() {return 3*10+1+3;}
    void print_10features();
    svm_model* getModel() {
      return rtModel.getModel();
    }
    size_t getRTFeat() {
      return rtModel.getRTFeat();
    }
    static size_t numDOCFeatures() {
      return 4;
    }
    void copyDOCparameters(DescriptionOfCorrect& other) {
      //    avgPI = other.avgPI; avgDM = other.avgDM; rtW = other.rtW; numRTFeat = other.numRTFeat;
      avgPI = other.avgPI;
      avgDM = other.avgDM;
      rtModel.copyModel(other.getModel());
      rtModel.setNumRtFeat(other.getRTFeat());
    }
    double estimateRT(double* features) {
      return rtModel.estimateRT(features);
    }

  protected:
    double avgPI, avgDM;
    vector<PSMDescription> psms;
    //  vector<double> rtW;
    double c, gamma, epsilon;
    RTModel rtModel;
    static string isoAlphabet;
    static float pKiso[7];
    static float pKN, pKC;
    static unsigned int docFeatures;
};

#endif /*DESCRIPTIONOFCORRECT_H_*/
