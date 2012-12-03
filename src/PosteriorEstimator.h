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

#ifndef POSTERIORESTIMATOR_H_
#define POSTERIORESTIMATOR_H_
#include<vector>
#include<string>
#include<utility>
#include "LogisticRegression.h"
using namespace std;


class PosteriorEstimator {
  public:
    PosteriorEstimator(){};
    virtual ~PosteriorEstimator(){};
    bool parseOptions(int argc, char** argv);
    string greeter();
    int run();
    static void estimatePEP(vector<pair<double, bool> >& combined,
                            double pi0, vector<double>& peps,
			      bool include_negative = false);
    static void estimatePEPGeneralized(vector<pair<double, bool> >& combined,
					 vector<double>& peps,
					 bool include_negative = false);
    static void estimate(vector<pair<double, bool> >& combined,
                         LogisticRegression& lr);
    static void getPValues(const vector<pair<double, bool> >& combined,
                           vector<double>& p);
    static void getQValues(double pi0,
                           const vector<pair<double, bool> >& combined,
                           vector<double>& q);
    static void getQValuesFromP(double pi0, const vector<double>& p,
                                vector<double>& q);
    static void getQValuesFromPEP(const vector<double>& pep,
                                vector<double>& q);
    static double estimatePi0(vector<double>& p,
                              const unsigned int numBoot = 100);
    static void setReversed(bool status) {
		reversed = status;
    }
    static void setGeneralized(bool general) {
		competition = general;
		assert(!(general && pvalInput));
    }
    static void setNegative(bool negative) {
	        includeNegativesInResult = negative;
    }
protected:
    void finishStandalone(vector<pair<double, bool> >& combined,
                          const vector<double>& peps,
                          const vector<double>& p, double pi0);
    void finishStandaloneGeneralized(vector<pair<double, bool> >& combined,
                          const vector<double>& peps);
    static void binData(const vector<pair<double, bool> >& combined,
                        vector<double>& medians,
                        vector<unsigned int>& negatives, vector<
                            unsigned int> & sizes);
    string targetFile, decoyFile;
    static bool reversed, pvalInput, competition,includeNegativesInResult;
    string resultFileName;
};

#endif /*POSTERIORESTIMATOR_H_*/
