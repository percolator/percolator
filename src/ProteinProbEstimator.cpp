/*
 * ProteinProbEstimator.cpp
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

#include <iostream>
#include <iomanip>
#include <fstream>
#include <limits>
#include "ProteinProbEstimator.h"
#include "ProteinProbEstimatorHelper.h"
#include "Globals.h"


ProteinProbEstimator::ProteinProbEstimator(double alpha_par, double beta_par) {
  peptideScores = 0;
  proteinGraph = 0;
  gamma = 0.5;
  alpha = alpha_par; // 0.1;
  beta = beta_par; // 0.01;
}

/**
 * sets alpha and beta to default values (avoiding the need for grid search)
 */
void ProteinProbEstimator::setDefaultParameters(){
  alpha = 0.1;
  beta = 0.01;
}

bool ProteinProbEstimator::initialize(Scores* fullset){
  peptideScores = fullset;
  populateProteinsToPeptidesTable(fullset, this);
  bool scheduleGridSearch;
  if(alpha ==-1 || beta ==-1) scheduleGridSearch = true;
  else scheduleGridSearch = false;
  return scheduleGridSearch;
}

/**
 * choose the set of parameters that jointly maximizes the ROC50 score (the
 * average sensitivity when allowing between zero and 50 false positives) and
 * minimizes the mean squared error (MSE_FDR) from an ideally calibrated
 * probability.
 * minimize: (1 − λ) MSE_FDR − λ ROC50 with λ = 0.15
 * range [0.01,0.76] at resolution of 0.05 for α
 * range [0.00,0.80] at resolution 0.05 for β
 */
void ProteinProbEstimator::gridSearchAlphaBeta(){
  // a λ approaching 1.0 will shift the emphasis to the most accurate model,
  // and a λ approaching 0.0 will result in a more calibrated model
  double lambda = 0.15;
  // before the search starts, the best point seen so far is artificially set
  // to the ideal absolute maximum (we are minimizing!)
  gridPoint bestSoFar;
  bestSoFar.objectiveFnValue = numeric_limits<double>::max();
  double lower_a=0.01, upper_a=0.76;
  double lower_b=0.01, upper_b=0.81;
  // if a parameter had previously been set (from command line) exclude it from
  // the grid search
  if(alpha != -1) lower_a = upper_a = alpha;
  if(beta != -1) lower_b = upper_b = beta;

  if(VERB > 1) cerr << endl << "             ";
  for(double b=log(lower_b); b<=log(upper_b); b+=0.5){
    if(VERB > 1) cerr << "beta = " << fixed<<std::setprecision(3) << exp(b) << "  ";
  }
  if(VERB > 1) cerr << endl;
  for(double a=log(lower_a); a<=log(upper_a); a+=0.5){
    if(VERB > 1) cerr << "alpha=" << fixed<<std::setprecision(3) << exp(a);
    for(double b=log(lower_b); b<=log(upper_b); b+=0.5){
      gridPoint current = gridPoint(exp(a),exp(b));
      current.calculateObjectiveFn(lambda, this);
      if(VERB > 1) {
        if(isinf(current.objectiveFnValue)) cerr << "  _infinity_  ";
        else cerr << "  " << current.objectiveFnValue << " ";
      }
      if(current<bestSoFar) bestSoFar = current;
    }
    if(VERB > 1) cerr << endl;
  }
  // the search is concluded: set the parameters
  if(bestSoFar.objectiveFnValue == numeric_limits<double>::max()){
    cerr << "ERROR: it was not possible to estimate values for parameters alpha and beta.\n"
        << "Please invoke Percolator with -a and -b option to set them manually.";
  }
  alpha = bestSoFar.alpha;
  beta = bestSoFar.beta;
}

/**
 * Calculate protein level probabilities. By default the parameters alpha and
 * beta will be estimated by grid search. If the function is invoked with
 * gridSearch set to false, whatever values for alpha and beta were previously
 * set will be used. If no values were set, the dafault will be enforced.
 *
 * @param fullset set of unique peptides with scores computed by Percolator
 * @param gridSearch indicate whether the values of alpha and beta parameters
 * should be estimated by grid search
 */
fidoOutput ProteinProbEstimator::calculateProteinProb(bool gridSearch){
  srand(time(NULL)); cout.precision(8); cerr.precision(8);
  // by default, a grid search is executed to estimate the values of the
  // parameters alpha and beta

  if(gridSearch) {
    if(VERB > 1) cerr << "Estimating parameters for the model by grid search\n";
    gridSearchAlphaBeta();
    if(VERB > 1) {
      cerr << "The following parameters have been chosen;\n";
      cerr << "alpha = " << alpha << endl;
      cerr << "beta = " << beta << endl;
      cerr << "Protein level probabilities will now be calculated\n";
    }
  }
  // at this point the parameters alpha and beta must have been initialized:
  // either statically set through command line or set after the grid search
  // or temporarily set in one of the grid search's iteration steps
  assert(alpha != -1);
  assert(beta != -1);

  //GroupPowerBigraph::LOG_MAX_ALLOWED_CONFIGURATIONS = ;
  delete proteinGraph;
  proteinGraph = new GroupPowerBigraph ( RealRange(alpha, 1, alpha),
      RealRange(beta, 1, beta), gamma );
  proteinGraph->read(peptideScores);
  proteinGraph->getProteinProbs();

  fidoOutput output = buildOutput(proteinGraph);
  // uncomment the following line to print protein level probabilities to file
  //writeOutputToFile(output, "/tmp/fido/6_final_fido_output.txt");
  return output;
}

/**
 * output protein level probabilites results in xml format
 *
 * @param os stream to which the xml is directed
 * @param output object containing the output to be written to file
 */
void ProteinProbEstimator::writeOutputToXML(string xmlOutputFN,
    const fidoOutput& output){
  ofstream os;
  os.open(xmlOutputFN.data(), ios::app);
  // append PROTEINs
  os << "  <proteins>" << endl;
  // for each probability
  for (int k=0; k<output.size(); k++) {
    Array<string> protein_ids = output.protein_ids[k];
    // for each protein with a certain probability
    for(int k2=0; k2<protein_ids.size(); k2++) {
      string protein_id = protein_ids[k2];
      // check wether is a decoy...
      double probOfDecoy = isDecoyProbability(protein_id, this);
      // if it's not a decoy, output it. If it is, only output if the option
      // isOutXmlDecoys() is activated
      if (probOfDecoy<0.5 || (probOfDecoy>0.5 && Scores::isOutXmlDecoys())) {
        os << "    <protein p:protein_id=\"" << protein_id << "\"";
        if (Scores::isOutXmlDecoys()) {
          if(probOfDecoy>0.5) os << " p:decoy=\"true\"";
          else  os << " p:decoy=\"false\"";
        }
        os << ">" << endl;
        os << "      <pep>" << output.peps[k] << "</pep>" << endl;
        os << "      <q_value>" << output.qvalues[k] << "</q_value>" << endl;
        writeXML_writeAssociatedPeptides(protein_id, os, proteinsToPeptides);
        os << "    </protein>" << endl;
      }
    }
  }
  os << "  </proteins>" << endl << endl;
  os.close();
}

/**
 * writes protein weights to cerr.
 *
 * @param proteinGraph proteins and associated probabilities to be outputted
 */
void ProteinProbEstimator::writeOutput(const fidoOutput& output) {
  cerr << endl;
  int size = output.size();
  for (int k=0; k<size; k++) {
    if (Scores::isOutXmlDecoys())
      cerr << output.peps[k] << " " << output.protein_ids[k] << endl;
    else {
      // filter decoys
      Array<string> filtered;
      Array<string>::Iterator protIt = output.protein_ids[k].begin();
      // only keep proteins associated with a majority of non-decoy psms
      for(; protIt< output.protein_ids[k].end(); protIt++){
        if(isDecoyProbability(*protIt, this)<0.5)
          filtered.add(*protIt);
      }
      if(filtered.size()>0){
        cerr << output.peps[k] << " " << filtered << endl;
      }
    }
  }
}
