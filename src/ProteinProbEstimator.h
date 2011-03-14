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
    Array<double> peps;
    Array< Array<string> > protein_ids;
    fidoOutput(Array<double> peps_par, Array<Array<string> > protein_ids_par) {
      peps = peps_par;
      protein_ids = protein_ids_par;
    }
};

class ProteinProbEstimator {
  public:
    double gamma;
    double alpha;
    double beta;
    ProteinProbEstimator(double alpha, double beta);
    virtual ~ProteinProbEstimator() {
      delete proteinGraph;
    }
    bool initialize(Scores* fullset);
    fidoOutput calculateProteinProb(bool gridSearch);
    void writeOutput(fidoOutput output);
    void writeOutputToXML(string xmlOutputFN);
    map<string, vector<ScoreHolder*> > proteinsToPeptides;
  private:
    void gridSearchAlphaBeta();
    GroupPowerBigraph* proteinGraph;
    Scores* peptideScores;
};

#endif /* PROTEINPROBESTIMATOR_H_ */
