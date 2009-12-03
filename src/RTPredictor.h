
#ifndef RTPREDICTOR_H_
#define RTPREDICTOR_H_

#include <string>
#include "PSMDescription.h"
#include "Normalizer.h"
#include "RTModel.h"
#include <utility>

#define VERSION "1.0"

class RTPredictor
{
	public:
		RTPredictor();
		~RTPredictor();
		string greeter();
		bool parseOptions(int argc, char **argv);
		void loadTrainFile();
		void loadTestFile();
		void printPsms(vector<PSMDescription> & psms);
		void printRunInformation();
		void removeRedundantPeptides(vector<PSMDescription> & psms);
		void removeAberrantPeptides(vector<PSMDescription> & psms);
		void run();
		void initFeaturesTable(const unsigned int numberRecords, vector<PSMDescription> & psms, double * retentionFeatures,
								size_t noFeatures = -1);
		void estimateRetentionTime(vector<PSMDescription> & psms);
		double computeMS(vector<PSMDescription> & psms);
		double computePearsonCorrelation(vector<PSMDescription> & psms);
		double computeSpearmanCorrelation(vector<PSMDescription> & psms);
		void writeOutputFile(vector<PSMDescription> & psms);
		void writeRetentionTimeFile(const char *, vector<PSMDescription> & psms);
		void loadBestModel();
		static vector<string> getModelFiles();
		// find the a and b that fit best the observed retention time using the least squared method
		static void findLeastSquaresSolution(const vector<PSMDescription> & psms, double & a, double & b);
		static void unNormalizeRetentionTimes(vector<PSMDescription> & psms);
		void writeFeaturesFile(const char* file, vector<PSMDescription> & psms, bool unnormalized);
		// check is child is a peptide included in parent and with the same rt
		static bool isChildOf(PSMDescription & child, PSMDescription & parent);
		// check if psm has a parent peptide in the training or test set
		bool hasParent(PSMDescription  psm);
		// push the decaying peptides at the end of the vector; it returns the index of the first decaying peptide
		int pushBackDecayingPeptides(vector <PSMDescription> & psms);
		static bool isTryptic(string & pep);
		static string getMSPeptide(string & peptide);
		// write decay peptides starting at index index to the outDecay
		void addToPairVector(vector<PSMDescription> psms, bool value, vector< pair<pair<PSMDescription, bool>,bool> > & psmPairs);
		void writeDecayToFile(vector< pair <pair<PSMDescription, bool>, bool> > & psms);
		// return indices of decay peptides
		void removeDecays(vector< pair< pair<PSMDescription, bool>, bool> > & psms);
		void removeNonTryptic(vector<PSMDescription> & psms);

	protected:
	    // the difference inhydrophobicity used to detect aberrant peptides
	    static float diff_hydrophobicity;
		// path to the library
		static string libPath;
		// true if the non tryptic peptides are to be removed
		bool removeNTryptic;
		// true if the test file includes the observed rt
		bool testIncludesRT;
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
		// file including the decaying peptides
		string decayingPeptidesFile;
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
		// remove aberrant peptides from the test set? (peptides that are included in a longer peptide with a similar retention time)
		// this option should be set ONLY when the test set includes retention times
		bool removeDecaying;
};

#endif /* RTPREDICTOR_H_ */
