/*
 * ProteinProbEstimator.h
 *
 *  Created on: Feb 25, 2011
 *      Author: tomasoni
 */

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

 *******************************************************************************/

#ifndef PROTEINPROBESTIMATOR_H_
#define PROTEINPROBESTIMATOR_H_

#include "GroupPowerBigraph.h"

struct fidoOutput {
    fidoOutput(Array<double> peps_par, Array<Array<string> > protein_ids_par,
        Array<double> estimQvalues_par, Array<double> empirQvalues_par,
        unsigned int proteinsAtThr1_par, unsigned int proteinsAtThr2_par,
        unsigned int totTargets_par, unsigned int totDecoys_par,
        double pi_0_par, double alpha_par, double beta_par,
        bool wellFormed_par) {
      peps = peps_par;
      protein_ids = protein_ids_par;
      estimQvalues = estimQvalues_par;
      empirQvalues = empirQvalues_par;
      proteinsAtThr1 = proteinsAtThr1_par;
      proteinsAtThr2 = proteinsAtThr2_par;
      totTargets = totTargets_par;
      totDecoys = totDecoys_par;
      pi_0 = pi_0_par;
      alpha = alpha_par;
      beta = beta_par;
      wellFormed = wellFormed_par;
    }
    Array<double> peps;
    Array< Array<string> > protein_ids;
    Array<double> estimQvalues;
    Array<double> empirQvalues;
    const static double threshold1 = 0.015;
    const static double threshold2 = 0.1;
    unsigned int proteinsAtThr1;
    unsigned int proteinsAtThr2;
    unsigned int totTargets;
    unsigned int totDecoys;
    double pi_0;
    double alpha;
    double beta;
    bool wellFormed;
};

class ProteinProbEstimator {
  public:
    double gamma;
    double alpha;
    double beta;
    ProteinProbEstimator(double alpha, double beta);
    virtual ~ProteinProbEstimator();
    bool initialize(Scores* fullset);
    void setDefaultParameters();
    fidoOutput run(bool startGridSearch);
    void writeOutput(const fidoOutput& output);
    void writeOutputToXML(string xmlOutputFN, const fidoOutput& output);
    static string printCopyright();
    static void testGridRanges();
    void plotQValues(const fidoOutput& output);
    void plotRoc(const fidoOutput& output, int N);
    void printStatistics(const fidoOutput& output);
    map<string, vector<ScoreHolder*> > proteinsToPeptides;
    unsigned int numberDecoyProteins;
    unsigned int numberTargetProteins;
    const static bool gridSearchDebugPlotting;
    const static bool debugginMode;
    const static bool logScaleSearch;
    const static bool tiesAsOneProtein;
    const static bool usePi0;
  private:
    void gridSearchAlphaBeta();
    GroupPowerBigraph* proteinGraph;
    Scores* peptideScores;
};

#endif /* PROTEINPROBESTIMATOR_H_ */
