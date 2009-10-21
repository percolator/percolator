/*
 * RTModel.h
 *
 */
#ifndef RTMODEL_H_
#define RTMODEL_H_

#include <vector>
#include <string>
#include "PSMDescription.h"
#include "Normalizer.h"
#include "svm.h"

using namespace std;
#define NO_FEATURE_GROUPS 12

// types of grid
typedef enum {NO_GRID, NORMAL_GRID, FINE_GRID} GridType;
typedef enum {SIMPLE_EVAL, K_FOLD_CV} EvaluationType;
struct gridInformation
{
	vector<double> gridGamma;
	vector<double> gridC;
	vector<double> gridEpsilon;
};

class RTModel
{
	public:
		RTModel();
		~RTModel();
		void calcRetentionFeatures(PSMDescription &psm);
		void calcRetentionFeatures(vector<PSMDescription> &psms);
		static void fillFeaturesAllIndex(const string& pep, double *features);
		static double* fillFeaturesIndex(const string& peptide, const float *index, double *features);
		//static double* conformationalPreferences(const string& peptide, double *features);
		static double* amphipathicityHelix(const float *index, const string& peptide, double *features);
		static double indexSum(const float* index, const string& peptide);
		static double indexAvg(const float* index, const string& peptide);
		static double indexNearestNeigbourPos(const float* index, const string& peptide);
		static double indexNearestNeigbourNeg(const float* index, const string& peptide);
		static inline double indexN(const float *index, const string& peptide);
		static inline double indexC(const float *index, const string& peptide);
		static inline double indexNC(const float *index, const string& peptide);
		static double* indexPartialSum(const float* index, const string& peptide, const size_t win, double *features);
		static double* fillAAFeatures(const string& pep, double *feat);
		static double* fillPTMFeatures(const string& pep, double *feat);
		static double  bulkinessSum(const string& peptide);
		/*static double  hydrophobicMoment(const string& peptide, const double angle);*/
		static int noConsecKRDENQ(const string& peptide);
		static double getNoPtms(string pep);
		void copyModel(svm_model* from);
		void trainRetention(vector<PSMDescription>& trainset);
		void trainRetention(vector<PSMDescription>& trainset, const double C, const double gamma, const double epsilon, int noPsms);
		double testRetention(vector<PSMDescription>& testset);
		double estimateRT(double * features);
		double computeKfoldCV(const vector<PSMDescription> & psms, const double gamma, const double epsilon, const double c);
		double computeSimpleEvaluation(const vector<PSMDescription> & psms, const double gamma, const double epsilon, const double c);
		void trainSVM(vector<PSMDescription> & psms);
		void setGridType(const GridType & g);
		string getGridType();
		string getEvaluationType();
		void saveSVRModel(const char* fileName, Normalizer *p);
		void loadSVRModel(const char* fileName, Normalizer *p);
        static size_t totalNumRTFeatures();
        // features for the three hydrophobicity indices, 3 features for ptms, peptide size, bulkiness,
        // no of consec KRDNEQ
	    static size_t minimumNumRTFeatures() {return 3*14 + 3 + 1 + 1 + 1;}
		size_t getRTFeat() {return numRTFeat;}
		void setNumRtFeat(const size_t nRtFeat) {numRTFeat = nRtFeat;}
		static void setDoKlammer(const bool switchKlammer);
		svm_model* getModel() {return model;}
		static string aaAlphabet,isoAlphabet;
		bool isModelNull() {if (model == NULL) return true; else return false;}
		void setCalibrationFile(const string calFile) {calibrationFile = calFile;}
		void setK(const int newk) {k = newk;}
		void setEvaluationType(const EvaluationType e) {eType = e;}
		void loadSVRModel(string modelFile, Normalizer * theNormalizer);
		void saveSVRModel(string modelFile, Normalizer * theNormalizer);
		void setSelectFeatures(const int sf);
		int getNoFeaturesToCalc() {return noFeaturesToCalc;}
		int getSelectFeatures() {return selected_features;}
		void printFeaturesInUse(ostringstream & oss);
		int getSelect(int sel_features, int max, size_t *finalNumFeatures);
		void destroyModel() {svm_destroy_model(model); }
	protected:
		// remove redundant peptides from the test set?
		bool removeRedundant;
		// number of features used to generate the model
		size_t numRTFeat;
		// number of features to be calculated (depends on the selected features)
		size_t noFeaturesToCalc;
		// svm model
		svm_model *model;
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
		static float krokhin_index['Z'-'A'+1],krokhin100_index['Z'-'A'+1],krokhinC2_index['Z'-'A'+1],krokhinTFA_index['Z'-'A'+1],
		             hessa_index['Z'-'A'+1],kytedoolittle_index['Z'-'A'+1], aa_weights['Z'-'A'+1], bulkiness['Z'-'A' + 1];
		// conformational preferences of aa
		//static float alpha_helix['Z'-'A'+1], beta_sheet['Z'-'A'+1];
		// the groups of features to be used
		static string feature_groups[NO_FEATURE_GROUPS];
		// how many features are in each group?
		static int no_features_per_group[NO_FEATURE_GROUPS];
		// features to be used (bit-wise operation are used to decode, till ex 2 means use just krokhin100_index, 3 means use
		int selected_features;
		// krokhin100_index and krokhin_index etc)
		static string ptmAlphabet;
		static bool doKlammer;
		double stepFineGrid;
		size_t noPointsFineGrid;
};

#endif /* RTMODEL_H_ */



