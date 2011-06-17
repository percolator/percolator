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

#include "boost/lexical_cast.hpp"
using namespace std;

///////////////////////////////////////////////////////////////////////////////
// DEBUGGING FUNCTIONS FOR ProteinProbEstimator
///////////////////////////////////////////////////////////////////////////////

struct ProteinDebugger{
    ProteinDebugger();
    ~ProteinDebugger();
    static void runScript(const string& script);
    static void plotQValues(const fidoOutput&, ProteinProbEstimator*);
    static void plotRoc(const fidoOutput&, ProteinProbEstimator*, int N);
    static void testGridRanges();
};

/**
 * routine that sets permissions and runs a script
 */
void ProteinDebugger::runScript(const string& script){
  string allowX = "chmod +x " + script;
  system(allowX.c_str());
  system(script.c_str());
}

/** plot to file the number of target proteins identified as a function of the
 *  q-value (empirical and estimated)
 */
void ProteinDebugger::plotQValues(const fidoOutput& output,
    ProteinProbEstimator* estimator){
  double threshold = 0.1;
  ostringstream a; a << fixed << setprecision(5) << output.alpha;
  ostringstream b; b << fixed << setprecision(5) << output.beta;
  //files that will contain the plot points
  string estim = string(TEMP_DIR) +
      "a" + a.str() + "_" + "b" + b.str() + "_" + "estimQValues.dat";
  string empir = string(TEMP_DIR) +
      "a" + a.str() + "_" + "b" + b.str() + "_" + "empirQValues.dat";
  ofstream o_est(estim.c_str());
  ofstream o_emp(empir.c_str());
  int targetsCount(0);
  int decoysCount(0);
  double previousEmpirQval(0), currentEmpirQval(0);

  // count targets and decoys as you go down the list of qvalues
  for(int k=0; k<output.estimQvalues.size(); k++){
    targetsCount += ProteinHelper::countTargets(
        output.protein_ids[k], estimator);
    decoysCount += ProteinHelper::countDecoys(
        output.protein_ids[k], estimator);
    // output estimated q-values
    if(output.estimQvalues[k] <= threshold){ // only plot below q-value thresh
      o_est << fixed << setprecision(6) <<
          output.estimQvalues[k] <<"\t"<< targetsCount <<"\n";
    }
    // output empirical q-values
    if(ProteinProbEstimator::usePi0){
      currentEmpirQval = (double)decoysCount/targetsCount*output.pi_0;
    } else {
      currentEmpirQval = (double)decoysCount/targetsCount;
    }
    if(currentEmpirQval>1.0) currentEmpirQval=1.0;
    if(currentEmpirQval<previousEmpirQval) currentEmpirQval=previousEmpirQval;
    double stored = output.empirQvalues[k];
    assert(abs(currentEmpirQval-stored)<1e-10);
    if (currentEmpirQval <= threshold)
      o_emp << fixed << setprecision(6) <<
      currentEmpirQval <<"\t"<< targetsCount <<"\n";
    previousEmpirQval = currentEmpirQval;
  }
  o_emp.close();
  o_est.close();
  // create gnuplot script
  string script = string(TEMP_DIR) +
      "a" + a.str() + "_" + "b" + b.str() + "_" + "plotQValues.sh";
  string image = string(TEMP_DIR) +
        "a" + a.str() + "_" + "b" + b.str() + "_" + "qValues.png";
  ofstream o_script(script.c_str());
  o_script << "#!/usr/bin/gnuplot\n\n";
  o_script << "set key right bottom\n";
  o_script << "set title \"alpha=" + a.str() + ", beta=" + b.str() + "\"\n";
  o_script << "set ylabel \"numb target proteins\"\n";
  o_script << "plot \""+estim+"\" using 1:2 title \"estimated q-values\" with lines,\\\n";
  o_script << "\""+empir+"\" using 1:2 title \"empirical q-values\" with lines\n";
  o_script << "set term png\n";
  o_script << "set output '" + image + "'\n";
  o_script << "replot";
  o_script.close();
  // invoke gnuplot
  runScript(script);
}

/** plot to file the number of decoy proteins identified as a function of the
 *  q-value (estimated)
 */
void ProteinDebugger::plotRoc(const fidoOutput& output,
    ProteinProbEstimator* estimator, int N){
  ostringstream a; a << fixed << setprecision(5) << output.alpha;
  ostringstream b; b << fixed << setprecision(5) << output.beta;
  string roc = string(TEMP_DIR) +
      "a" + a.str() + "_" + "b" + b.str() + "_" + "rocPlot.dat";
  ofstream o_roc(roc.c_str());
  int targetsCount(0), decoysCount(0);
  // count decoys as you go down the list of qvalues
  for(int k=0; k<output.estimQvalues.size(); k++){
    if(decoysCount>N) break;
    int targetsAtQValue =
        ProteinHelper::countTargets(output.protein_ids[k], estimator);
    targetsCount += targetsAtQValue;
    decoysCount += output.protein_ids[k].size() - targetsAtQValue;
    if(k % 10 == 0){
      // output number of decoys at current q-value
      o_roc << fixed << setprecision(6) <<
          targetsCount <<"\t"<< decoysCount <<"\n";
    }
  }
  o_roc.close();
  // create gnuplot script
  string script = string(TEMP_DIR) +
      "a" + a.str() + "_" + "b" + b.str() + "_" + "plotRoc.sh";
  string image = string(TEMP_DIR) +
        "a" + a.str() + "_" + "b" + b.str() + "_" + "roc.png";
  ofstream o_script(script.c_str());
  o_script << "#!/usr/bin/gnuplot\n\n";
  o_script << "set key right bottom\n";
  o_script << "set title \"alpha=" + a.str() + ", beta=" + b.str() + "\"\n";
  o_script << "set ylabel \"# targets\"\n";
  o_script << "set xlabel \"# decoys\"\n";
  o_script << "plot \"" + roc + "\" using 2:1 title \"roc\" with lines\n";
  o_script << "set term png\n";
  o_script << "set output '" + image + "'\n";
  o_script << "replot";
  o_script.close();
  // invoke gnuplot
  runScript(script);
}

/** prints the grid
 */
void ProteinDebugger::testGridRanges(){
  Grid* grid;
  if(ProteinProbEstimator::logScaleSearch)
    grid = new Grid();
  else grid = new Grid(0.001,0.05,0.0001,0.005,0.001,0.0002);
  grid->current_a = grid->getLower_a();
  for(; grid->current_a<=grid->getUpper_a(); grid->updateCurrent_a()){
    grid->current_b = grid->getLower_b();
    for(; grid->current_b<=grid->getUpper_b(); grid->updateCurrent_b()){
    }
  }
}

#endif /* PROTEINPROBESTIMATORDEBUGGER_H_ */
