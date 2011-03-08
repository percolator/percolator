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

class ProteinProbEstimator {
  public:
    virtual ~ProteinProbEstimator() {
      delete instance;
      delete proteinGraph;
    }
    static ProteinProbEstimator* getInstance();
    int calculateProteinProb(Scores* fullset, bool gridSearch=false);
    void writeXML(ofstream& os);
  private:
    ProteinProbEstimator() {
      assert(instance==0);
      instance = this;
      proteinGraph = 0;
      gamma = 0.5;
      alpha = -1;
      beta = -1;
    }
    void gridSearchAlphaBeta();
    static ProteinProbEstimator* instance;
    GroupPowerBigraph* proteinGraph;
    map<string, vector<ScoreHolder*> > proteinsToPeptides;
    double gamma;
    double alpha;
    double beta;
};

#endif /* PROTEINPROBESTIMATOR_H_ */
