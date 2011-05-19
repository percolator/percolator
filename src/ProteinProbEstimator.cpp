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
#include <fstream>
#include "ProteinProbEstimator.h"
#include "ProteinProbEstimatorHelper.h"
#include "ProteinProbEstimatorDebugger.h"

/**
 * enable plot to debug q-value calculation and Roc50 calculation: methods
 * ProteinProbEstimator::plotQValues and ProteinProbEstimator::plotRoc will be
 * invoked
 */
const bool ProteinProbEstimator::gridSearchDebugPlotting = true;
/**
 * extra debug information is printed to file in ProteinProbEstimatorHelper.h
 * and BasicBigraph.cpp
 */
const bool ProteinProbEstimator::debugginMode = true;
/**
 * during grid search, should the grid values be sampled exponentially farther
 * away from each other? (or just linearly?)
 */
const bool ProteinProbEstimator::logScaleSearch=true;
/**
 * treat ties as if it were one protein
 */
const bool ProteinProbEstimator::tiesAsOneProtein = false;
/**
 * use pi_0 value when calculating empirical q-values
 */
const bool ProteinProbEstimator::usePi0 = true;


ProteinProbEstimator::ProteinProbEstimator(double alpha_par, double beta_par) {
  peptideScores = 0;
  proteinGraph = 0;
  gamma = 0.5;
  alpha = alpha_par;
  beta = beta_par;
  numberDecoyProteins = 0;
  numberTargetProteins = 0;
}

ProteinProbEstimator::~ProteinProbEstimator(){
  delete proteinGraph;
}

/**
 * sets alpha and beta to default values (avoiding the need for grid search)
 */
void ProteinProbEstimator::setDefaultParameters(){
  alpha = 0.1;
  beta = 0.01;
}

/**
 * initializes the ProteinProbEstimator by storing a pointer to percolator's
 * peptide level results and by scheduling a grid search in case one or both
 * of the parameters alpha and beta have not been set.
 */
bool ProteinProbEstimator::initialize(Scores* fullset){
  // percolator's peptide level statistics
  peptideScores = fullset;
  // populated a hash table to retrieve peptides given a protein name
  ProteinHelper::populateProteinsToPeptidesTable(fullset, this);
  bool scheduleGridSearch = false;
  if(alpha ==-1 || beta ==-1) scheduleGridSearch = true;
  return scheduleGridSearch;
}

/**
 * Calculate protein level probabilities.
 *
 * @param startGridSearch indicates whether the values of alpha and beta
 * parameters should be estimated by grid search. In Caller::run(): intended to
 * be initialized by the return value of ProteinProbEstimator::initialize(); in
 * GridPoint::calculateObjectiveFn: always set to false (avoid recursively
 * starting nested grid searches)
 */
fidoOutput ProteinProbEstimator::run(bool startGridSearch){
  srand(time(NULL)); cout.precision(8); cerr.precision(8);
  // by default, a grid search is executed to estimate the values of the
  // parameters alpha and beta
  if(startGridSearch) {
    if(VERB > 1) {
      cerr << "The parameters for the model will be estimated by grid search."
          << endl;
    }
    gridSearchAlphaBeta();
    if(VERB > 1) {
      cerr << "\nThe following parameters have been chosen;\n";
      cerr << "alpha = " << alpha << endl;
      cerr << "beta  = " << beta << endl;
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
  proteinGraph->read(peptideScores, ProteinProbEstimator::debugginMode);
  proteinGraph->getProteinProbs();

  fidoOutput output = ProteinHelper::buildOutput(proteinGraph, this);
  if(ProteinProbEstimator::debugginMode) {
    // print protein level probabilities to file
    ProteinHelper::writeOutputToFile(output, string(WRITABLE_DIR) +
        "6percolator_final_fido_output.txt");
  }
  return output;
}

/**
 * choose the set of parameters that jointly maximizes the ROC50 score (the
 * average sensitivity when allowing between zero and 50 false positives) and
 * minimizes the mean squared error (MSE_FDR) from an ideally calibrated
 * probability.
 * minimize: (1 − λ) MSE_FDR − λ ROC50 with λ = 0.35
 * range [0.01,0.76] at resolution of 0.05 for α
 * range [0.00,0.80] at resolution 0.05 for β
 */
void ProteinProbEstimator::gridSearchAlphaBeta(){
  // a λ approaching 1.0 will shift the emphasis to the most accurate model,
  // and a λ approaching 0.0 will result in a more calibrated model
  //Grid grid = Grid(0.001,0.05,0.0001,0.005,0.001,0.0002);
  Grid grid = Grid();
  // if a parameter had previously been set (from command line) limit the
  // search in the corresponding direction
  if(alpha != -1) grid.limitSearch(Grid::alpha,alpha);
  if(beta != -1) grid.limitSearch(Grid::beta,beta);
  // evaluate the objective function in the default location to compare
  // performances
  grid.compareAgainstDefault(this);
  // start the search
  grid.current_a = grid.getLower_a();
  for(; grid.current_a<=grid.getUpper_a(); grid.updateCurrent_a()){
    grid.current_b = grid.getLower_b();
    for(; grid.current_b<=grid.getUpper_b(); grid.updateCurrent_b()){
      grid.toCurrentPoint(); // create point
      grid.calculateObjectiveFn(this);
      grid.updateBest();
    }
  }
  if(VERB > 2){
    cerr << "\n\nThe search was completed; debugging details follow:\n"
        << "alpha\t\t" << "beta\t\t" << "MSE_FDR\t\t" << "ROC50\t\t"
        << "ObjFn\t\t" << "@0.015\t\t" << "@0.1\n"
        << grid.debugInfo.str();
  }
  // the search is concluded: set the parameters
  if(grid.wasSuccessful()){
    grid.setToBest(this);
  } else {
    cerr << "ERROR: it was not possible to estimate values for parameters alpha and beta.\n"
        << "Please invoke Percolator with -a and -b option to set them manually.";
  }
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
  for (int k=0; k<output.peps.size(); k++) {
    Array<string> protein_ids = output.protein_ids[k];
    // for each protein with a certain probability
    for(int k2=0; k2<protein_ids.size(); k2++) {
      string protein_id = protein_ids[k2];
      // check wether is a decoy...
      double probOfDecoy =
          ProteinHelper::isDecoyProbability(protein_id, this);
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
        os << "      <q_value>" << output.estimQvalues[k] << "</q_value>\n";
        ProteinHelper::writeXML_writeAssociatedPeptides(
            protein_id, os, proteinsToPeptides);
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
  cout << "PEP\t\t" << "est qvalues\t" << "emp qvalues\t" << "proteins\n";
  int size = output.peps.size();
  for (int k=0; k<size; k++) {
    if (Scores::isOutXmlDecoys())
      cout << scientific << setprecision(7) << output.peps[k] << "\t"
      << scientific << setprecision(7) << output.estimQvalues[k]<< "\t"
      << scientific << setprecision(7) << output.empirQvalues[k]<< "\t"
      << scientific << setprecision(7) << output.protein_ids[k] << endl;
    else {
      // filter decoys
      Array<string> filtered;
      Array<string>::Iterator protIt = output.protein_ids[k].begin();
      // only keep proteins associated with a majority of non-decoy psms
      for(; protIt< output.protein_ids[k].end(); protIt++){
        if(ProteinHelper::isDecoyProbability(*protIt, this)<0.5)
          filtered.add(*protIt);
      }
      if(filtered.size()>0){
        cout << scientific << setprecision(7) << output.peps[k] << "\t"
            << scientific << setprecision(7) << output.estimQvalues[k]<< "\t"
            << scientific << setprecision(7) << output.empirQvalues[k]<< "\t"
            << scientific << setprecision(7) << filtered << endl;
      }
    }
  }
}

string ProteinProbEstimator::printCopyright(){
  ostringstream oss;
  oss << "Copyright (c) 2008-9 University of Washington. All rights reserved.\n"
      << "Written by Oliver R. Serang (orserang@u.washington.edu) in the\n"
      << "Department of Genome Sciences at the University of Washington.\n";
  return oss.str();
}

void ProteinProbEstimator::testGridRanges(){
  ProteinDebugger::testGridRanges();
}

/** plot to file the number of target proteins identified as a function of the
 *  q-value (empirical and estimated)
 */
void ProteinProbEstimator::plotQValues(const fidoOutput& output){
  if(!ProteinProbEstimator::gridSearchDebugPlotting) return;
  assert(output.peps.size()!=0);
  ProteinDebugger::plotQValues(output,this);
}

/** plot to file the number of decoy proteins identified as a function of the
 *  q-value (estimated)
 */
void ProteinProbEstimator::plotRoc(const fidoOutput& output, int N){
  if(!ProteinProbEstimator::gridSearchDebugPlotting) return;
  assert(output.peps.size()!=0);
  ProteinDebugger::plotRoc(output,this,N);
}

void ProteinProbEstimator::printStatistics(const fidoOutput& output){
  ProteinHelper::printStatistics(output);
}
