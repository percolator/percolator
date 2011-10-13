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
        unsigned int targetsAtThr1_par, unsigned int targetsAtThr2_par,
        unsigned int decoysAtThr2_par,
        unsigned int totTargets_par, unsigned int totDecoys_par,
        double pi_0_par, double alpha_par, double beta_par,
        bool wellFormed_par) {
      peps = peps_par;
      protein_ids = protein_ids_par;
      estimQvalues = estimQvalues_par;
      empirQvalues = empirQvalues_par;
      targetsAtThr1 = targetsAtThr1_par;
      targetsAtThr2 = targetsAtThr2_par;
      decoysAtThr2 = decoysAtThr2_par;
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
    const static double threshold1 = 0.01;
    const static double threshold2 = 0.05;
    unsigned int targetsAtThr1;
    unsigned int targetsAtThr2;
    unsigned int decoysAtThr2;
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
    const static double default_alpha = 0.1;
    const static double default_beta = 0.01;
    ProteinProbEstimator(double alpha, double beta, bool tiesAsOneProtein = false
			 ,bool usePi0 = false, bool outputEmpirQVal = false);
    virtual ~ProteinProbEstimator();
    bool initialize(Scores* fullset);
    void setDefaultParameters();
    fidoOutput run(bool startGridSearch);
    void writeOutputToStream(const fidoOutput& output, ostream& stream);
    void writeOutputToXML(const fidoOutput& output, string xmlOutputFN);
    static string printCopyright();
    static void testGridRanges();
    void plotQValues(const fidoOutput& output);
    void plotRoc(const fidoOutput& output, int N);
    void printStatistics(const fidoOutput& output);
    
    void setTiesAsOneProtein(bool tiesAsOneProtein);
    void setUsePio(bool usePi0);
    void setOutputEmpirQval(bool outputEmpirQVal);
    bool getTiesAsOneProtein();
    bool getUsePio();
    bool getOutputEmpirQval();
    
    map<string, vector<ScoreHolder*> > proteinsToPeptides;
    unsigned int numberDecoyProteins;
    unsigned int numberTargetProteins;
    const static bool debugginMode;
    const static bool logScaleSearch;
    const static bool outputPEPs;
    
  private:
    void gridSearchAlphaBeta();
    GroupPowerBigraph* proteinGraph;
    Scores* peptideScores;
    bool tiesAsOneProtein;
    bool usePi0;
    bool outputEmpirQVal;
};

#endif /* PROTEINPROBESTIMATOR_H_ */
