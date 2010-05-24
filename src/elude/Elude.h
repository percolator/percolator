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

 *****************************************************************************/
#ifndef RTPREDICTOR_H_
#define RTPREDICTOR_H_

#include <string>
#include <utility>
#include "PSMDescription.h"
#include "Normalizer.h"
#include "EludeModel.h"
#include "LTSRegression.h"
#include "Enzyme.h"

//#define VERSION "1.0"

class RTPredictor {
  public:
    RTPredictor();
    ~RTPredictor();
    // print introductory message
    string greeter();
    // command line options
    bool parseOptions(int argc, char** argv);
    // print information about the current run
    void printRunInformation();
    // functions to load and process input data
    void loadTrainFile();
    void processTrainData();
    void loadTestFile();
    void processTestData();
    // remove redundant peptides; the peptide with lowest rt is kept
    void removeRedundantPeptides(vector<PSMDescription> & psms,
                                 const char* dataset);
    // initialize the table storing the retention features
    void initFeaturesTable(const unsigned int numberRecords, vector<
        PSMDescription> & psms, double* retentionFeatures,
                           size_t noFeatures = -1);
    // remove CID fragments
    void removeAberrantPeptides(vector<PSMDescription> & psms);
    // check if child is CID fragments with parent as parental peptide
    bool isChildOf(PSMDescription& child, PSMDescription& parent);
    // check if psm has a parent peptide in the input data
    bool hasParent(PSMDescription psm);
    // push the decaying peptides at the end of the vector; it returns the index of the first decaying peptide
    int pushBackDecayingPeptides(vector<PSMDescription> & psms);
    // for a peptide given in the format X.Y.Z, return Y
    static string getMSPeptide(string& peptide);
    static bool fileExists(string& strFilename);
    void addToPairVector(vector<PSMDescription> psms, bool value, vector<
        pair<pair<PSMDescription, bool> , bool> > & psmPairs);
    // write CID fragments to a file
    void
        writeDecayToFile(
                         vector<pair<pair<PSMDescription, bool> , bool> > & psms);
    // remove CID fragments
    void
        removeSourceCIDs(
                         vector<pair<pair<PSMDescription, bool> , bool> > & psms);
    // remove non-tryptic peptides
    void removeNonEnzymatic(vector<PSMDescription> & psms, string text);
    // return the observed and predicted retention times as a pair of vectors
    pair<vector<double> , vector<double> >
        getRTs(vector<PSMDescription> & psms);
    // train the SVM regressor
    void trainSVRRegressor();
    // load a model
    void loadSVRModel();
    // get all the model files from the library
    static vector<string> getModelFiles();
    // load a model from the library
    void loadLibModel();
    // load the model that best fits the training data from the library
    void loadBestModel();
    // add the current model to the library
    void addModelLibrary();
    // implement least squares regression
    static void
        findLeastSquaresSolution(const vector<PSMDescription> & psms,
                                 double& a, double& b);
    // predict retention time
    void predictRTs();
    void estimateRetentionTime(vector<PSMDescription> & psms);
    static void unNormalizeRetentionTimes(vector<PSMDescription> & psms);
    // compute mean-squares, Pearson and Spearman correlations, window where 95% peptides are located
    double computeMS(vector<PSMDescription> & psms);
    double computePearsonCorrelation(vector<PSMDescription> & psms);
    double computeSpearmanCorrelation(vector<PSMDescription> & psms);
    double computeWindow(vector<PSMDescription> & psms);
    // writing functions
    void writeOutputFile(vector<PSMDescription> & psms);
    void
        writeRetentionTimeFile(const char*, vector<PSMDescription> & psms);
    void writeFeaturesFile(const char* file,
                           vector<PSMDescription> & psms,
                           bool unnormalized);
    // print a vector of psms
    void printPsms(vector<PSMDescription> & psms);
    // main function
    void run();
    bool operator()(PSMDescription& psm) {
      return theEnzyme->isEnzymatic(psm.getPeptide());
    }

  protected:
    // EXPERIMENTAL
    //double C;
    double a, b;
    // the difference in hydrophobicity used to detect CID fragments
    static float diff_hydrophobicity;
    // how many peptides should be included in the time window reported (by default 95%)
    static float fractionPeptides;
    // path to the library
    static string libPath;
    // true if the non tryptic peptides are to be removed
    bool removeNEnzymatic;
    // true if the test file includes the observed rt
    bool testIncludesRT;
    // linearly adjustment?
    bool linearAdjust;
    // remove redundant peptides from the test set?
    bool removeRedundant;
    // remove CID fragments from the test set? (peptides that are included in a longer peptide with a similar retention time)
    // this option should be set ONLY when the test set includes retention times
    bool removeDecaying;
    // true if the model is to be loaded from the library?
    bool autoModelSelection;
    // is the model to be added to the library
    bool addLibModel;
    // least trimmed squares regression
    LTSRegression lts;
    // information about svr model
    RTModel model;
    // file containing the training data
    string trainFile;
    // file including the test data
    string testFile;
    // output file for the predictions
    string outputFile;
    // file to save the created model
    string saveModelFile;
    // file to load the model from
    string loadModelFile;
    // file where the CID fragments should be stored
    string decayingPeptidesFile;
    // file where the hydrophobicity index should be saved
    string hydrophobicityFile;
    // train and test peptides
    vector<PSMDescription> trainPsms;
    vector<PSMDescription> testPsms;
    // pointer to the feature table of training and testing data
    double* trainRetentionFeatures;
    double* testRetentionFeatures;
    // normalizer
    Normalizer* theNormalizer;
    // enzyme
    Enzyme* theEnzyme;
};

#endif /* RTPREDICTOR_H_ */
