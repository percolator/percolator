/*
 * ProteinProbEstimatorDebugger.h
 *
 *  Created on: Apr 19, 2011
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

#ifndef PROTEINPROBESTIMATORDEBUGGER_H_
#define PROTEINPROBESTIMATORDEBUGGER_H_

#include "ProteinProbEstimatorDebugger.h"

///////////////////////////////////////////////////////////////////////////////
// DEBUGGING FUNCTIONS FOR ProteinProbEstimator
///////////////////////////////////////////////////////////////////////////////

struct DebugHelper{
    DebugHelper();
    ~DebugHelper();
    static void plotQValues(const fidoOutput&, ProteinProbEstimator*);
    static void testGridRanges();
};

/** plot to file the number of target proteins identified as a function of the
 *  q-value (empirical and estimated)
 */
void DebugHelper::plotQValues(const fidoOutput& output,
    ProteinProbEstimator* estimator){
  //files that will contain the plot points
  string estim = string(WRITABLE_DIR) + "estimQValuesPlot.dat";
  string empir = string(WRITABLE_DIR) + "empirQValuesPlot.dat";
  ofstream o_est(estim.c_str());
  ofstream o_emp(empir.c_str());
  int targetsCount, decoysCount = 0;
  // count targets and decoys as you go down the list of qvalues
  for(int k=0; k<output.qvalues.size(); k++){
    int targetsAtQValue = countTargets(output.protein_ids[k], estimator);
    targetsCount += targetsAtQValue;
    decoysCount += output.protein_ids[k].size() - targetsAtQValue;
    if(k % 10 == 0){
      // output estimated q-values
      o_est << fixed << setprecision(6) <<
          output.qvalues[k] <<"\t"<< targetsCount <<"\n";
      // output empirical q-values
      o_emp << fixed << setprecision(6) <<
          output.pi_0*decoysCount/targetsCount <<"\t"<< targetsCount <<"\n";
    }
  }
}

/** prints the grid
 */
void DebugHelper::testGridRanges(){
  Grid grid = Grid();
  grid.current_a = grid.getLower_a();
  for(; grid.current_a<=grid.getUpper_a(); grid.updateCurrent_a()){
    grid.current_b = grid.getLower_b();
    for(; grid.current_b<=grid.getUpper_b(); grid.updateCurrent_b()){
    }
  }
}

#endif /* PROTEINPROBESTIMATORDEBUGGER_H_ */
