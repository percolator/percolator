/*
 * @ Created by L. Moruz
 * Nov 23rd, 2009
 * Implementation of Fast LTS regression for p = 2 (y = ax + b) and (aproximately) n <= 600. If n > 600
 * then the SVM should be trained on this dataset (more reliable results)
 * The method is implemented as explained in "Computing LTS Regression for Large Data Sets" by
 * Peter J. Rousseeuw and Katrien van Driessen, Data Mining and Knowledge Discovery, Vol. 12, No. 1. (January 2006), pp. 29-45.
 */

#ifndef LTSREGRESSION_H_
#define LTSREGRESSION_H_

#include <utility>
#include <vector>

using namespace std;

struct dataPoint
{
	double x;
	double y;
	// absoluet value of residual
	double absr;
};

class LTSRegression
{
	public:
		LTSRegression();
		~LTSRegression();
		// set the data points used for regression
		void setData(vector<double> & x, vector<double> & y);
		// construct an initial random p-subset; the elements of data will be sorted, and the first h elements are returned
		vector<dataPoint> getInitialHSubset();
    	// fill the residuals for all the data points using the line ax + b, with a = par.first, b = par.second
		void fillResiduals(pair<double, double> par);
		// fit a line using the least-squares method using h; return the a and b of the model
		pair<double, double> fitLSLine(vector<dataPoint> h);
		// perform a C-step starting with h(build the regression line, compute residuals, sort data according to residuals)
		// it returns the new h
		vector<dataPoint> performCstep(vector<dataPoint> h);
		// main function to detect the parameters a,b
		void runLTS();
		// predict the y values of x
		vector<double> predict(vector<double> & x);
		// print the values of the data points
		void printDataPoints();
		// get the data points
		vector<dataPoint> getDataPoints() {return data;};
		// calculate the squares of the residuals
		double calculateQ();
		// get the first 3 best h-subsets
		vector< vector<dataPoint> > getBestHSubsets();
		// print vector
		void printVector(vector<dataPoint> v);

	protected:
		// the max number of initial sets H1 generated; be default we use 500 (as suggested in the article)
	    // if the number of possible subsets is smaller than 500, we generate them all
	    static int noSubsets;
	    // maximum difference to acheive convergence
	    static double epsilon;
	    // coverage (number of points used to generate the regression line); the default value will be 0.75*n
	    int h;
	    // the points
	    vector <dataPoint> data;
	    // the regression coefficients (y = ax + b => a, b are the coefficients)
	    pair<double, double> regCoefficients;
};

#endif /* LTSREGRESSION_H_ */
