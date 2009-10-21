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
#include <sys/types.h>
#include <dirent.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <cstring>
#include <math.h>
#include <assert.h>
#include "Option.h"
#include "Globals.h"
#include "RTPredictor.h"
#include "PSMDescription.h"
#include <sstream>
#include <ctime>

using namespace std;

// path to the library
string RTPredictor::libPath = "../lib/";

RTPredictor::RTPredictor() : trainFile(""), testFile(""),outputFile(""),
			 saveModelFile(""), loadModelFile(""), trainRetentionFeatures(NULL), testRetentionFeatures(NULL), theNormalizer(NULL),
			 logFile(""), removeRedundant(false), autoModelSelection(false), a(1.0), b(0.0), linearAdjust(true), addLibModel(false)
{
	RTModel::setDoKlammer(false);
	Normalizer::setType(Normalizer::UNI);
}

RTPredictor::~RTPredictor()
{
	if (trainRetentionFeatures)
	{
		// delete the retention feature table
		delete[] trainRetentionFeatures;
		trainRetentionFeatures = NULL;
	}

	if (testRetentionFeatures)
	{
		// delete the retention feature table
		delete[] testRetentionFeatures;
		testRetentionFeatures = NULL;
	}

	if (theNormalizer)
		delete theNormalizer;
}

// introductory message
string RTPredictor::greeter()
{
	ostringstream oss;

  	oss << "rtPredictor version " << VERSION << ", ";
  	oss << "Build Date " << __DATE__ << " " << __TIME__ << endl;
  	oss << "Distributed under MIT License" << endl;
  	oss << "Written by Lukas Käll (lukas.kall@cbr.su.se) in the" << endl;
  	oss << "Center for Biomembrane Reasearch." << endl;
  	oss << "Dept. of Biochemistry, Stockholm University, Stockholm." << endl;
  	oss << "Usage:" << endl;
  	oss << "   rtPredictor [options]" << endl << endl;

  	return oss.str();
}

// parse the command line
bool RTPredictor::parseOptions(int argc, char **argv)
{
	ostringstream intro;
  	intro << greeter() << endl;
  	CommandLineParser cmd(intro.str());

  	// define available options
  	cmd.defineOption("v", "verbose", "Set verbosity of output: 0=no processing info, 5=all, default is 2", "level");
  	cmd.defineOption("n", "no_grid", "Specifies that no calibration of parameters should be carried out", "", TRUE_IF_SET);
  	cmd.defineOption("f", "fine_grid", "Specifies that the calibration of parameters will be carried out using a coarse grid followed by a fine grid", "", TRUE_IF_SET);
  	cmd.defineOption("p", "plain_evaluation", "The performance of a model is estimated by training the model on 3/4 of training data and testing on the 1/4 left.", "", TRUE_IF_SET);
  	cmd.defineOption("t", "train", "Specifies the file including the data to train the model", "filename");
  	cmd.defineOption("l", "log", "Specifies the log file", "filename");
  	cmd.defineOption("c", "calibration", "Files storing information about calibration steps", "filename");
  	cmd.defineOption("e", "evaluate", "Specifies the file including the test data", "filename");
  	cmd.defineOption("s", "save-model", "Specifies the file in which the model will be saved", "filename");
  	cmd.defineOption("m", "load-model", "Specifies a file including a SVM model to be loaded", "filename");
  	cmd.defineOption("o", "out", "Create an output file including information regarding the current run", "filename");
  	cmd.defineOption("u", "unique", "Remove all redundant peptides from the test set", "", TRUE_IF_SET);
	cmd.defineOption("k", "k_fold", "Specify the number of folds for cross validation", "value");
	cmd.defineOption("r", "rt_feat", "Specify a number between 2^0 and 2^9 to select the rt feature groups to be used ", "value");
	cmd.defineOption("a", "auto", "Specifies that the model that matches the best the training data will be used for the test set ", "", TRUE_IF_SET);
	cmd.defineOption("d", "append", "Append model to library ", "", TRUE_IF_SET);
	cmd.defineOption("j", "no_linear_adjust", "Specifies that the model will not be linearly adjusted ", "", TRUE_IF_SET);

	// parse command line
	cmd.parseArgs(argc, argv);

	// process options
	if (cmd.optionSet("c"))
  		model.setCalibrationFile(cmd.options["c"]);

	if (cmd.optionSet("l"))
  		logFile = cmd.options["l"];

  	if (cmd.optionSet("v"))
  		Globals::getInstance()->setVerbose(cmd.getInt("v",0,10));

  	if (cmd.optionSet("t"))
  	{
  		trainFile = cmd.options["t"];
  	}

  	if (cmd.optionSet("e"))
  	{
  		testFile = cmd.options["e"];
  	}

  	if (cmd.optionSet("s"))
  		saveModelFile = cmd.options["s"];

  	if (cmd.optionSet("m"))
  	  	loadModelFile = cmd.options["m"];

	if (cmd.optionSet("o"))
		outputFile = cmd.options["o"];

	if (cmd.optionSet("n"))
		model.setGridType(NO_GRID);

	if (cmd.optionSet("f"))
		model.setGridType(FINE_GRID);

	if (cmd.optionSet("p"))
		model.setEvaluationType(SIMPLE_EVAL);

	if (cmd.optionSet("u"))
		removeRedundant = true;

	if (cmd.optionSet("k"))
		model.setK(cmd.getInt("k",1,100));

	if (cmd.optionSet("r"))
		model.setSelectFeatures(cmd.getInt("r",1,1000));

	if (cmd.optionSet("a"))
		autoModelSelection = true;

	if (cmd.optionSet("d"))
	{
		addLibModel = true;
	}

	if (cmd.optionSet("j"))
		linearAdjust = false;

  	return true;
}

// load the input file fileName and store the information in a vector of PSMs
// if includeRT = false, then only the peptides are included in the file, otherwise both peptide and retention time are included
void RTPredictor::loadInputFile(const string & fileName, vector<PSMDescription> & psms, bool includesRT)
{
	// retention time
	double retTime;
	string peptide;
	// opens the file for reading
	ifstream in(fileName.c_str(), ios::in);
	// number of peptides in the file

	cout << "Loading file " << fileName << endl;

	// check if the file was open successfully
	if (!in)
	{
		cerr << "Unable to open " << fileName << endl;
		cerr << "Please check the file location and try again " << endl;
		cerr << "Execution aborted "<< endl;
		exit(-1);
	}

	if (includesRT)
	{
		// read the observed retention time and the peptide sequence
		while (in >> peptide >> retTime)
			// add a new psm
			psms.push_back( PSMDescription(peptide, retTime) );
	}
	else
	{
		// read the observed retention time and the peptide sequence
		while (in >> peptide)
			// add a new psm
			psms.push_back( PSMDescription(peptide, 0.0) );
	}

	cout << psms.size() << " peptides were loaded." << endl;
	cout << "Done." << endl << endl;
}

// print the PSMs
void RTPredictor::printPsms(vector<PSMDescription> & psms)
{
	cout << endl << "Printing..." << endl;

	vector<PSMDescription>::iterator it;
	for(it = psms.begin(); it != psms.end(); it ++)
		cout << (*it) << endl;

	cout << "Done." << endl << endl;
}

void RTPredictor::printRunInformation()
{
	ostringstream oss;

	oss << endl << "----------------------------------------" << endl;
	oss <<  "Running rtPredictor " << VERSION << " with the following options: " << endl;
	if (!trainFile.empty())
		oss << " *Training data: " << trainFile << endl;

	if (!testFile.empty())
		oss << " *Test data: " << testFile << endl;

	if (!loadModelFile.empty())
		oss << " *Loading model from file: " << loadModelFile << endl;

	if (!saveModelFile.empty())
		oss << " *Saving model to file: " << saveModelFile << endl;

	if (autoModelSelection)
	{
		oss << " *The SVR model will be selected automatically" << endl;
		if (linearAdjust)
			oss << " *Linear fitting will be carried out " << endl;
		else
			oss << " *No linear fitting will be carried out " << endl;
	}
	if (!outputFile.empty())
		oss << " *Output file: " << outputFile << endl;

	if (removeRedundant)
		oss << " *Redundant peptides will be removed from the test set" << endl;
	else
		oss << " *Redundant peptides are allowed in the test set" << endl;

	if (!trainFile.empty() && !autoModelSelection)
	{
		oss << " *Calibration of parameters: " << model.getGridType() << endl;
		oss << " *Evaluation type: " << model.getEvaluationType() << endl;
	}
	oss << "-----------------------------------------" << endl << endl;


	cout << oss.str();
}

// remove all the duplicate peptides
void RTPredictor::removeRedundantPeptides(vector<PSMDescription> & psms)
{
	cout << "Removing redundant peptides... " << endl;
	int initial, final;

	initial = psms.size();

	// sort first so we keep the peptide with the lowest retention time
	sort(psms.begin(), psms.end(),less<PSMDescription>());
	psms.resize(distance(psms.begin(), unique(psms.begin(), psms.end())));

	cout << (initial - psms.size()) << " peptides were removed " << endl;
	cout << "Dataset includes now " << psms.size() << " peptides." << endl;
	cout << "Done." << endl<< endl;

}

// read all the files from the libPath
vector<string> RTPredictor::getModelFiles()
{
	vector<string> modelFiles;
	DIR *dp;
	struct dirent *dirp;
	string name;
	size_t extension;

	if ((dp = opendir(libPath.c_str())) == NULL)
	{
		cout << "Unable to open " << libPath << endl;
		cout << "Execution aborted. " << endl;
		exit(-1);
	}

	while ((dirp = readdir(dp)) != NULL)
	{
		name = string(dirp->d_name);
		extension = name.rfind(".model");
		if (extension == (name.length() - 6))
			modelFiles.push_back(libPath + name);
	}

	closedir(dp);

	return modelFiles;
}


// load the best model out of the library
void RTPredictor::loadBestModel()
{
	// variables to store best correlation
	double bestCorr, corr;
	// name of the file including the best model
	string bestModelFile;
	vector<string> modelFiles;
	double *retFeatures;

	cout << "Checking available model files..." << endl;
	// if no training data available, we cannot choose a model
	if (trainPsms.empty())
	{
		cout << "No training data available. " << endl;
		cout << "No model will be loaded. " << endl;
		return;
	}

	// get all the model files from the library
	modelFiles = getModelFiles();
	// if no model files in the library, there is nothing to do at this step
	if (modelFiles.empty())
	{
		cout << "No model file available. " << endl;
		return;
	}
	else
	{
		cout << "The library includes " << modelFiles.size() << " models." << endl;
		cout << "Done." << endl << endl;
	}

	bestCorr = 0.0;
	bestModelFile ="";
	for(int i = 0; i < modelFiles.size(); i++)
	{
		// load the model and information about the normalizer
		cout << "Processing " << modelFiles[i] << "..." << endl << "-------" << endl ;

		cout << "Loading model " << modelFiles[i] << " ..." << endl;
		model.loadSVRModel(modelFiles[i], theNormalizer);
		cout << "Done." << endl << endl;

		// if the model was not loaded then just give a warning and go to the next one
		if (model.isModelNull())
		{
			cout << "Warning: Model file " << modelFiles[i] << " was not loaded correctly." << endl;
			cout << "This model will not be considered." << endl;
		}
		else
		{
			// normalize the retention times; the scaling parameters should have already been set when the model was loaded
			PSMDescription::normalizeRetentionTimes(trainPsms);

			// calculate the retention features
			model.calcRetentionFeatures(trainPsms);

			// normalize the retention features; the scaling parameters should have already been set when loading the model
	        cerr << "FIXME: Replace with new normalizer" << endl;
			// theNormalizer->normalizeSet(trainPsms);

			// estimate the retention time
			estimateRetentionTime(trainPsms);

			// compute correlation
			corr = computeCorrelation(trainPsms);

			if (corr > bestCorr)
			{
				bestCorr = corr;
				bestModelFile = modelFiles[i];
			}

			// unnormalize the retention time
			for(int i = 0; i < trainPsms.size(); i++)
				trainPsms[i].retentionTime = PSMDescription::unnormalize(trainPsms[i].retentionTime);
		}
	}

	if (bestModelFile != "")
	{
		cout << endl << "----------------------" << endl ;
		cout << "Best model: " << bestModelFile << " with correlation " << bestCorr << endl;
		cout << endl << "Initialize model..." << endl;
		model.loadSVRModel(bestModelFile, theNormalizer);
		// normalize the retention times; the scaling parameters should have already been set when the model was loaded
		PSMDescription::normalizeRetentionTimes(trainPsms);
		// calculate the retention features
		model.calcRetentionFeatures(trainPsms);
		// normalize the retention features; the scaling parameters should have already been set when loading the model
        cerr << "FIXME: Replace with new normalizer" << endl;
//		theNormalizer->normalizeSet(trainPsms);
		// estimate the retention time
		estimateRetentionTime(trainPsms);
	}
	else
		cout << "No model loaded " << endl;
}

// find the best line that fits the data (basically the coefficients a, b)
void RTPredictor::findLeastSquaresSolution(const vector<PSMDescription> & psms, double & a, double & b)
{
	// the x is the predicted value, while the y is the observed one
	double avgX, avgY, sumX = 0.0, sumY = 0.0, sumXSquared = 0.0, sumXY = 0.0;
	double ssxx, ssxy;
	int n = psms.size();
	double x, y;

	for(int i = 0; i < n; ++i)
	{
		x = psms[i].predictedTime;
		y = psms[i].retentionTime;
		//cout << x << ", " << y << endl;
		sumX += x;
		sumY += y;
		sumXSquared += x * x;
		sumXY += x * y;
	}

	avgX = sumX / (double)n;
	avgY = sumY / (double)n;
	ssxx = sumXSquared - (n * avgX * avgX);
	ssxy = sumXY - (n * avgX * avgY);

	a = ssxy / ssxx;
	b = avgY - (b * avgX);
}

// unnormalize both observed and predicted retention times for a set of peptides
void RTPredictor::unNormalizeRetentionTimes(vector<PSMDescription> & psms)
{
	for(int i = 0; i < psms.size(); ++i)
	{
		psms[i].retentionTime = PSMDescription::unnormalize(psms[i].retentionTime);
		psms[i].predictedTime = PSMDescription::unnormalize(psms[i].predictedTime);
	}
}

void RTPredictor::run()
{
	double trainingTime = 0.0, testingTime = 0.0;
	// it will be used to check whether the test psms were loaded
	bool loadedTestPsms = false;

	// redirect standard output to the log file (if such file was provided)
	if (!logFile.empty())
		freopen (logFile.c_str(), "w", stdout);

	// print information about the current run
	printRunInformation();

	// get a normalizer
	theNormalizer = Normalizer::getNormalizer();
	theNormalizer->resizeVecs(RTModel::totalNumRTFeatures());

	// if standalone application
	if (!autoModelSelection)
	{
		///////////////////////////////// GENERATE MODEL
		// generate a model (either by training if train data given, either by loading a model)
		// if training data is given, the model is generated by training a SVR
		if (!trainFile.empty())
		{
			cout << "************* Load and Process training data ***************" << endl;
			// read the input data
			loadInputFile(trainFile, trainPsms, true);

			// remove redundant peptides
			removeRedundantPeptides(trainPsms);

			// this has to be included here so I can remove from the training set the peptides that appear in the test set as well
			// this will not be done in case model selection is to be performed
			if (!testFile.empty())
			{
				loadInputFile(testFile, testPsms, true);
				loadedTestPsms = true;

				// remove redundant peptides if the user wanted so
				if (removeRedundant)
					removeRedundantPeptides(testPsms);

				// remove from the train set the peptides that are in the test
				cout << "Removing from the train set the peptides featuring in both train and test sets..." << endl;
				vector<PSMDescription>::iterator it;
				int initialSize = trainPsms.size();
				for(it = testPsms.begin(); it != testPsms.end(); ++it)
				{
					trainPsms.erase(remove(trainPsms.begin(), trainPsms.end(), (*it)), trainPsms.end());
				}

				cout << (initialSize - trainPsms.size())<< " peptides were removed." << endl;
				cout << "Final training set includes " << trainPsms.size() << " peptides" << endl << endl;
			}

			// normalize the retention times; since sub and div are static members, their values will be set correctly
			// for the test set
			PSMDescription::setPSMSet(trainPsms);
			PSMDescription::normalizeRetentionTimes(trainPsms);

			// allocate memory for the feature table
			initFeaturesTable(trainPsms.size(), trainPsms, trainRetentionFeatures);

			// compute the retention features
			model.calcRetentionFeatures(trainPsms);

			// normalize all the features; note that the same scaling parameters will be used to normalize the test data if any
	        cerr << "FIXMET: Replace with new normalizer" << endl;
	        vector<double*> tmp;
	        vector<double*> tRetFeat = PSMDescription::getRetFeatures(trainPsms);
	        theNormalizer->setSet(tmp, tRetFeat, (size_t)0, model.getNoFeaturesToCalc());
	        theNormalizer->normalizeSet(tmp, tRetFeat);
//			theNormalizer->setPsmSet(trainPsms, model.getNoFeaturesToCalc());
//			theNormalizer->normalizeSet(trainPsms);

			// build the model
			cout << "************* Build the SVR model ***************" << endl;
			clock_t start = std::clock();
			model.trainSVM(trainPsms);
			trainingTime = ((std::clock() - start) / (double)CLOCKS_PER_SEC);
			cout << "Training lasted: "  << trainingTime << " seconds " << endl << endl;
		}

		// load a model
		if (!loadModelFile.empty())
		{
			cout << "************* Load the SVR model ***************" << endl;
			// if a training file was also provided and a model was already generated, this step is skipped
			if ((!trainFile.empty()) && (!model.isModelNull()))
			{
				cerr << "WARNING: both options -m and -t generate a model.\nSince a model was already generated using data from "
					 << trainFile << " option -m will be ignored" << endl << "Please use option -m alone to load an existent model\n";
			}
			else
			{
				cout << "Loading model from file " << loadModelFile << endl;

				// load the svm model
				//theNormalizer->resizeVecs(RTModel::totalNumRTFeatures());
				//theNormalizer->setNumberRetentionFeatures(RTModel::totalNumRTFeatures());
				///////////////////////////////////////////////////////77

				model.loadSVRModel(loadModelFile, theNormalizer);
				//model.setNumRtFeat((*theNormalizer->getNumRetFeatures()));

				cout << "Done" << endl;
			}
		}
	}
	else
	{
		// AUTOMATIC MODEL SELECTION
		// If no training data provided the first model available will be used
		if (trainFile.empty())
		{
			// load first available model
			// to be implemented
		}
		else
		{
			cout << "************* Load and Process training data ***************" << endl;
			// load input data
			loadInputFile(trainFile, trainPsms, true);

			// remove redundant peptides
			removeRedundantPeptides(trainPsms);

			// initialize the table storing the retention features
			initFeaturesTable(trainPsms.size(), trainPsms, trainRetentionFeatures, RTModel::totalNumRTFeatures());

			cout << endl << "************* Load best model for your data ***************" << endl;
			// load the model that fits best the data
			loadBestModel();

			// find the best line that fits the data
			if (linearAdjust)
			{
				cout << "Compute linear fitting parameters..." << endl;
				findLeastSquaresSolution(trainPsms, a, b);
				cout << "prediction = " << a << " * x + " << b << endl;
				cout << "Done." << endl;
			}
		}

	}
	////////////////////////// TEST OR SAVE a model
	// TEST if testing data provided
	if(!testFile.empty())
	{
		cout << endl << "************* Testing ***************" << endl;
		// first check that a retention model was generated
		if (model.isModelNull())
		{
			cout << "Error: No SVM model available " << endl;
			cout << "Please train a model using the -t option or load an existing model from a file using the"
					" -m option and try again" << endl;
			cout << "Execution aborted" << endl;
			exit(-1);
		}

		cout << endl;
		// if the file was not already loaded
		if (!loadedTestPsms)
		{
			// read the input file
			loadInputFile(testFile, testPsms, true);

			// remove redundant peptides if the user wanted so
			if (removeRedundant)
				removeRedundantPeptides(testPsms);
		}

		// normalize the retention times; the scaling parameters should have already been set when the model was generated or loaded
		PSMDescription::normalizeRetentionTimes(testPsms);

		// initialize the table storing the retention features
		initFeaturesTable(testPsms.size(), testPsms, testRetentionFeatures);

		// calculate the retention features
		model.calcRetentionFeatures(testPsms);

		// normalize the retention features; the scaling parameters should have already been set when generating the model
        cerr << "FIXMET: Replace with new normalizer" << endl;
        vector<double*> tmp;
 	    vector<double*> tRetFeat = PSMDescription::getRetFeatures(testPsms);
		theNormalizer->normalizeSet(tmp, tRetFeat);

		// predict retention time
		//clock_t start = std::clock();
		estimateRetentionTime(testPsms);
		//unnormalizeRetentionTimes(testPsms);
		// if the model was select automatically, then the prediction will be ax + b
		if (linearAdjust)
		{
			double ms = computeMS(testPsms);
			cout << "Adjusting linearly..." << endl;
			for(int i = 0; i < testPsms.size(); i++)
				testPsms[i].predictedTime = a * testPsms[i].predictedTime + b;
			cout << "Done." << endl << endl;
		}
		//testingTime = ((std::clock() - start) / (double)CLOCKS_PER_SEC);
		//cout << "Testing lasted: "  << testingTime << " seconds " << endl << endl;

		double correlation, ms;
		// compute correlation and ms
		correlation = computeCorrelation(testPsms);
		ms = computeMS(testPsms);
		cout << "--------------------------------------\n";
		cout << "Pearson Correlation = " << correlation << endl;
		cout << "MS = " << ms << endl;
		cout << "--------------------------------------\n\n";

		// write the output file
		if (!outputFile.empty())
			writeOutputFile(testPsms);
	}

	// save the model
	if (!saveModelFile.empty())
	{
		cout << "************* Save the SVR model ***************" << endl;
		if (model.isModelNull())
		{
			cout << "WARNING: No model was generated; please use options -m and -t to generate a model and then save it. \n"
				<< "Option -s has no effect;" << endl;
		}
		else
		{
			// SAVE A MODEL in the given file
			cout << "Saving SVR model to file " << saveModelFile << endl;
			model.saveSVRModel(saveModelFile, theNormalizer);

			cout << "Done." << endl;
		}
	}

	// adding to library
	if (addLibModel)
	{
		cout << "************* Add model o library ***************" << endl;

		if (model.isModelNull())
		{
			cout << "No model available to add to library." << endl;
		}
		else
		{
			string libModelFile;
			if (!trainFile.empty())
			{
				size_t found;
				found = trainFile.find_last_of("/");
				if (found != string::npos)
					libModelFile = trainFile.substr(found + 1, trainFile.length());
				else
				{
					found = trainFile.find_last_of("\\");
					if (found != string::npos)
						libModelFile = trainFile.substr(found + 1, trainFile.length());
					else
						libModelFile = trainFile;
				}
			}
			else
			{
				time_t currentTime;
				struct tm * timeInfo;

				currentTime = time(NULL);
				timeInfo = localtime (&currentTime);
				libModelFile = asctime(timeInfo);
				libModelFile = libModelFile.substr(4, (libModelFile.length()-10));

				size_t found = libModelFile.find_first_of(" ");
				while(found != string::npos)
				{
					libModelFile.replace(found, 1 ,"_");
					found = libModelFile.find_first_of(" ");
				}
			}

			libModelFile = libPath + libModelFile + ".model" ;
			cout << "Model name: " << libModelFile << endl;
			model.saveSVRModel(libModelFile, theNormalizer);
			cout << "Done." << endl;
		}
	}

	if (!logFile.empty())
		fclose(stdout);
}


/*
void RTPredictor::run()
{
	double trainingTime = 0.0, testingTime = 0.0;
	// it will be used to check whether the test psms were loaded
	bool loadedTestPsms = false;

	// redirect standard output to the log file (if such file was provided)
	if (!logFile.empty())
		freopen (logFile.c_str(), "w", stdout);

	// print information about the current run
	printRunInformation();

	// get a normalizer
	theNormalizer = Normalizer::getNormalizer();
	theNormalizer->resizeVecs(RTModel::totalNumRTFeatures());

	///////////////////////////////// GENERATE MODEL
	// generate a model (either by training if train data given, either by loading a model)
	// if training data is given, the model is generated by training a SVR
	if (!trainFile.empty())
	{
		// read the input data
		loadInputFile(trainFile, trainPsms, true);

		// remove redundant peptides
		removeRedundantPeptides(trainPsms);

		// this has to be included here so I can remove from the training set the peptides that appear in the test set as well
		if (!testFile.empty())
		{
			loadInputFile(testFile, testPsms, true);
			loadedTestPsms = true;

			// remove redundant peptides if the user wanted so
			if (removeRedundant)
				removeRedundantPeptides(testPsms);

			// remove from the train set the peptides that are in the test
			cout << "Removing from the train set the peptides featuring in both train and test sets..." << endl;
			vector<PSMDescription>::iterator it;
			int initialSize = trainPsms.size();
			for(it = testPsms.begin(); it != testPsms.end(); ++it)
			{
				trainPsms.erase(remove(trainPsms.begin(), trainPsms.end(), (*it)), trainPsms.end());
			}

			cout << (initialSize - trainPsms.size())<< " peptides were removed." << endl;
			cout << "Final training set includes " << trainPsms.size() << " peptides" << endl << endl;
		}

		// normalize the retention times; since sub and div are static members, their values will be set correctly
		// for the test set
		PSMDescription::setPSMSet(trainPsms);
		PSMDescription::normalizeRetentionTimes(trainPsms);

		// allocate memory for the feature table
		initFeaturesTable(trainPsms.size(), trainPsms, trainRetentionFeatures);

		// compute the retention features
		model.calcRetentionFeatures(trainPsms);

		// normalize all the features; note that the same scaling parameters will be used to normalize the test data if any
		theNormalizer->setPsmSet(trainPsms, model.getNoFeaturesToCalc());
		theNormalizer->normalizeSet(trainPsms);

		// build the model
		clock_t start = std::clock();
		model.trainSVM(trainPsms);
		trainingTime = ((std::clock() - start) / (double)CLOCKS_PER_SEC);
		cout << "Training lasted: "  << trainingTime << " seconds " << endl << endl;
	}

	// load a model
	if (!loadModelFile.empty())
	{
		// if a training file was also provided and a model was already generated, this step is skipped
		if ((!trainFile.empty()) && (!model.isModelNull()))
		{
			cout << "WARNING: both options -m and -t generate a model.\nSince a model was already generated using data from "
				 << trainFile << " option -m will be ignored" << endl << "Please use option -m alone to load an existent model\n";
		}
		else
		{
			cout << "Loading model from file " << loadModelFile << endl;

			// load the svm model
			//theNormalizer->resizeVecs(RTModel::totalNumRTFeatures());
			//theNormalizer->setNumberRetentionFeatures(RTModel::totalNumRTFeatures());
			///////////////////////////////////////////////////////77

			model.loadSVRModel(loadModelFile, theNormalizer);
			//model.setNumRtFeat((*theNormalizer->getNumRetFeatures()));

			cout << "Done" << endl;
		}
	}

	////////////////////////// TEST OR SAVE a model
	// TEST if testing data provided
	if(!testFile.empty())
	{
		// first check that a retention model was generated
		if (model.isModelNull())
		{
			cout << "Error: No SVM model available " << endl;
			cout << "Please train a model using the -t option or load an existing model from a file using the"
					" -m option and try again" << endl;
			cout << "Execution aborted" << endl;
			exit(-1);
		}

		cout << endl;
		// if the file was not already loaded
		if (!loadedTestPsms)
		{
			// read the input file
			loadInputFile(testFile, testPsms, true);

			// remove redundant peptides if the user wanted so
			if (removeRedundant)
				removeRedundantPeptides(testPsms);
		}

		// normalize the retention times; the scaling parameters should have already been set when the model was generated or loaded
		PSMDescription::normalizeRetentionTimes(testPsms);

		// initialize the table storing the retention features
		initFeaturesTable(testPsms.size(), testPsms, testRetentionFeatures);

		// calculate the retention features
		model.calcRetentionFeatures(testPsms);

		// normalize the retention features; the scaling parameters should have already been set when generating the model
		theNormalizer->normalizeSet(testPsms);

		// predict retention time
		clock_t start = std::clock();
		estimateRetentionTime(testPsms);
		testingTime = ((std::clock() - start) / (double)CLOCKS_PER_SEC);
		cout << "Testing lasted: "  << testingTime << " seconds " << endl;

		double correlation, ms;
		// compute correlation and ms
		correlation = computeCorrelation(testPsms);
		ms = computeMS(testPsms);
		cout << "--------------------------------------\n";
		cout << "Pearson Correlation = " << correlation << endl;
		cout << "MS = " << ms << endl;
		cout << "--------------------------------------\n\n";

		// write the output file
		if (!outputFile.empty())
			writeOutputFile(testPsms);
	}

	// save the model
	if (!saveModelFile.empty())
	{
		if (model.isModelNull())
		{
			cout << "WARNING: No model was generated; please use options -m and -t to generate a model and then save it. \n"
				<< "Option -s has no effect;" << endl;
		}
		else
		{
			// SAVE A MODEL in the given file
			cout << "Saving SVR model to file " << saveModelFile << endl;
			model.saveSVRModel(saveModelFile, theNormalizer);

			cout << "Done." << endl;
		}
	}

	if (!logFile.empty())
		fclose(stdout);
}*/

// allocating memory for the feature table
void RTPredictor::initFeaturesTable(const unsigned int numberRecords, vector<PSMDescription> & psms, double * retentionFeatures,
									size_t noFeatures)
{
	// number of rt related features
 	size_t noRTFeatures;
 	// copy of the pointer
 	double *ptr;

 	if (noFeatures == -1)
 		noRTFeatures = model.getNoFeaturesToCalc();
 	else
 		noRTFeatures = noFeatures;

 	cout << "Initializing features table for " << noRTFeatures << " features and " << numberRecords << " records ..." << endl;

 	// allocate memory to all the feature table
	retentionFeatures = new double[noRTFeatures * numberRecords];
	// get the beginning of the memory block
	ptr = retentionFeatures;

	if (!retentionFeatures)
	{
		cout << "Unable to allocate memory for the features table" << endl;
		cout << "Execution aborted" << endl;
		exit(-1);
	}

	// resize the size of the vector
  	psms.resize(numberRecords);

 	// divide the big allocated memory block to all rtPredictors (each one gets a line)
  	for(int i = 0; i < numberRecords; i++, ptr += noRTFeatures)
      		psms[i].retentionFeatures = ptr;

	cout << "Done." << endl << endl;
}

// estimates the retention time and fills in the predicted retention time field in each PSMDescription
void RTPredictor::estimateRetentionTime(vector<PSMDescription> & psms)
{
	assert(model.getModel());
	vector<PSMDescription>::iterator it;
	double predictedRT;

	cout << "Predicting retention features..." << endl;
	for(it = psms.begin(); it != psms.end(); ++it)
	{
		predictedRT =  model.estimateRT(it->retentionFeatures);
		//cout << predictedRT << " ";
		it->predictedTime = predictedRT;
	}

	cout << "Done." << endl << endl;
}

// calculate rms between observed and predicted retention times
double RTPredictor::computeMS(vector<PSMDescription> & psms)
{
	double ms = 0.0, diff = 0.0;
	int noPsms;

	cout << "Computing ms..." << endl;

	noPsms = psms.size();

	for(int i = 0; i < noPsms; i++)
	{
    	diff = psms[i].predictedTime - psms[i].retentionTime;
		ms += diff * diff;
  	}

  	ms = ms / noPsms;

	cout << "Done." << endl << endl;

  	return ms;
}

// calculate the correlation rho
double RTPredictor::computeCorrelation(vector<PSMDescription> & psms)
{
	double sumObserved = 0.0, sumPredicted = 0.0;
	double meanObserved = 0.0, meanPredicted = 0.0;
	double stdevObserved = 1, stdevPredicted = 1;
	double devObserved, devPredicted;
	double numerator = 0.0;
	double corr;
	int noPsms = psms.size();

	cout << "Computing correlation..." << endl;

	// calculate means
	for(int i = 0; i < noPsms; i++)
	{
		sumObserved += psms[i].retentionTime;
		sumPredicted += psms[i].predictedTime;
	}

	meanObserved = sumObserved / noPsms;
	meanPredicted = sumPredicted / noPsms;

	sumObserved = 0.0;
	sumPredicted = 0.0;
	// calculate stdevs
	for(int i = 0; i < noPsms; i++)
	{
		devObserved = psms[i].retentionTime - meanObserved;
		devPredicted = psms[i].predictedTime - meanPredicted;
		numerator += devObserved * devPredicted;
		sumObserved += pow(devObserved, 2);
		sumPredicted += pow(devPredicted, 2);
	}

	stdevObserved = sqrt(sumObserved / (noPsms - 1));
	stdevPredicted = sqrt(sumPredicted / (noPsms - 1));

	corr = numerator / ((noPsms - 1) * stdevObserved * stdevPredicted);

	cout << "r = " << corr << endl;
	cout << "Done." << endl << endl;

	return corr;
}

void RTPredictor::writeOutputFile(vector<PSMDescription> & psms)
{
	ofstream out(outputFile.c_str());
	vector<PSMDescription>::iterator it;
	double predictedRT, observedRT;

	out << "# File generated by rtPredictor " << endl;
	out << "# " << __DATE__  << " , " << __TIME__  << endl;
	out << "Peptide\tpredicted_retention_time\t[observed_retention_time]" << endl;

	for(it = psms.begin(); it != psms.end(); ++it)
	{
		predictedRT = PSMDescription::unnormalize(it->predictedTime);
		observedRT = PSMDescription::unnormalize(it->retentionTime);
		out << (it->peptide) << "\t" << predictedRT << "\t" << observedRT << "\n";
	}

	out.close();
}
