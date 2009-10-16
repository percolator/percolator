
#ifndef RTPREDICTOR_H_
#define RTPREDICTOR_H_

#include <string>
#include "PSMDescription.h"
#include "Normalizer.h"
#include "RTModel.h"

#define VERSION "1.0"

class RTPredictor
{
	public:
	RTPredictor();
	~RTPredictor();
	string greeter();
	bool parseOptions(int argc, char **argv);
	void loadInputFile(const string & fileName, vector<PSMDescription> & psms, bool includesRT);
	void printPsms(vector<PSMDescription> & psms);
	void printRunInformation();
	void removeRedundantPeptides(vector<PSMDescription> & psms);
	void run();
	void initFeaturesTable(const unsigned int numberRecords, vector<PSMDescription> & psms, double * retentionFeatures,
						    size_t noFeatures = -1);
	void estimateRetentionTime(vector<PSMDescription> & psms);
	double computeMS(vector<PSMDescription> & psms);
	double computeCorrelation(vector<PSMDescription> & psms);
	void writeOutputFile(vector<PSMDescription> & psms);
	void loadBestModel();
	static vector<string> getModelFiles();
	// find the a and b that fit best the observed retention time using the least squared method
	static void findLeastSquaresSolution(const vector<PSMDescription> & psms, double & a, double & b);
	static void unNormalizeRetentionTimes(vector<PSMDescription> & psms);

	protected:
		// path to the library
		static string libPath;
		// linearly adjustment?
		bool linearAdjust;
		// the coefficients a and b of the line that fits best the data (l = ax + b)
		double a, b;
		// should the model be searched in the library?
		bool autoModelSelection;
		// is the model to be added to the library
		bool addLibModel;
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
		// log file
		string logFile;
		// train and test peptides
		vector<PSMDescription> trainPsms;
		vector<PSMDescription> testPsms;
		// pointer to the feature table of training and testing data
		double *trainRetentionFeatures;
		double *testRetentionFeatures;
		// normalizer
		Normalizer *theNormalizer;
		// remove redundant peptides from the test set?
		bool removeRedundant;
};

#endif /* RTPREDICTOR_H_ */
