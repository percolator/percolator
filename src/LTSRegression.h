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
/*
 * @ Created by L. Moruz
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
	// absolute value of residual
	double absr;
};

class LTSRegression
{
	public:
		LTSRegression();
		~LTSRegression();
		// set the data points used for regression
		void setData(vector<double> & x, vector<double> & y);
		// construct an initial random p-subset; data will be sorted, and the first h elements are returned
		vector<dataPoint> getInitialHSubset();
    	// fill the absolute values of residuals for all the data points using the line ax + b, with a = par.first, b = par.second
		void fillResiduals(pair<double, double> par);
		// fit a line using the least-squares method using h; return the a and b of the model
		pair<double, double> fitLSLine(vector<dataPoint> h);
		// perform a C-step starting with h (build the regression line, compute abs(residuals), sort data according to abs(residuals))
		// it returns the new h
		vector<dataPoint> performCstep(vector<dataPoint> h);
		// predict the y values of x
		double predict(double x) { return ((regCoefficients.first * x) + regCoefficients.second); }
		// calculate the squares of the residuals
		double calculateQ();
		// apply LTS regression
		void runLTS();
		// get functions
		vector<dataPoint> getDataPoints() { return data; };
		pair<double, double> getRegCoefficients() { return regCoefficients; }
		// printing functions
		void printVector(vector<dataPoint> v);
		void printDataPoints();

	protected:
		// the max number of initial sets H1 generated; be default we use 500 (as suggested in the article)
	    static int noSubsets;
	    // maximum difference to acheive convergence
	    static double epsilon;
	    // coverage (number of points used to generate the regression line); the default value will be 0.75*n
	    int h;
	    // data points
	    vector <dataPoint> data;
	    // the regression coefficients (y = ax + b => a, b are the coefficients)
	    pair<double, double> regCoefficients;
};

#endif /* LTSREGRESSION_H_ */
