/*
 * ProteinProbEstimatorHelper.h
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

#ifndef PROTEINPROBESTIMATORHELPER_H_
#define PROTEINPROBESTIMATORHELPER_H_

#include <math.h>
#include <string>
#include <set>
#include <limits>
#include <iomanip>
#include <boost/unordered_set.hpp>
#include "Vector.h"
#include "Globals.h"
using namespace std;


///////////////////////////////////////////////////////////////////////////////
// HELPER FUNCTIONS FOR ProteinProbEstimator METHODS
///////////////////////////////////////////////////////////////////////////////

struct ProteinHelper{
    ProteinHelper();
    ~ProteinHelper();
    static unsigned int countTargets(const Array<string>& protein_ids,
        ProteinProbEstimator* estimator);
    static unsigned int countDecoys(const Array<string>& protein_ids,
        ProteinProbEstimator* estimator);
    static fidoOutput buildOutput(
        GroupPowerBigraph* proteinGraph, ProteinProbEstimator* estimator);
    static void writeOutputToFile(fidoOutput output, string fileName);
    static void writeXML_writeAssociatedPeptides(
        string& protein_id, ofstream& os,
        map<string, vector<ScoreHolder*> >& proteinsToPeptides);
    static void populateProteinsToPeptidesTable(
        Scores* fullset, ProteinProbEstimator* thisEstimator);
    static double isDecoyProbability(
        string protein_id, ProteinProbEstimator* estimator);
    static void printStatistics(const fidoOutput& output);
};

/**
 * given a list of protein ids returns the number of targets
 *
 * @param protein_ids
 * @param estimator object containing the complete lists of targets and decoys
 * for the experiment
 */
unsigned int ProteinHelper::countTargets(const Array<string>& protein_ids,
    ProteinProbEstimator* estimator){
  unsigned int countTargets(0);
  // for each protein in the list
  for(int p=0; p<protein_ids.size(); p++){
    // if the protein is a target...
    if(isDecoyProbability(protein_ids[p], estimator)<0.5){
      if(ProteinProbEstimator::tiesAsOneProtein) {
        // ... return in case ties are being counted as one protein
        return 1;
      } else {
        // ... otherwise, increment the counter and keep going
        if(isDecoyProbability(protein_ids[p], estimator)<0.5) countTargets++;
      }
    }
  }
  return countTargets;
}

/**
 * given a list of protein ids returns the number of decoys
 *
 * @param protein_ids
 * @param estimator object containing the complete lists of targets and decoys
 * for the experiment
 */
unsigned int ProteinHelper::countDecoys(const Array<string>& protein_ids,
    ProteinProbEstimator* estimator){
  // in case ties are being counted as one protein...
  if(ProteinProbEstimator::tiesAsOneProtein) {
    // ... return as soon as a decoy is found
    unsigned int countDecoys(0);
    for(int p=0; p<protein_ids.size(); p++){
      if(isDecoyProbability(protein_ids[p], estimator)>=0.5) return 1;
    }
    return 0; // zero if none were found
  } else {
    // ... otherwise count targets and subtract that to total
    return protein_ids.size() - countTargets(protein_ids, estimator);
  }
}

/**
 * after calculating protein level probabilities, the output is stored in a
 * dedicated structure that can be printed out or evaluated during the grid
 * search
 *
 * @param proteinGraph graph with proteins and corresponding probabilities
 * calculated by fido
 * @return output results of fido encapsulated in a fidoOutput structure
 */
fidoOutput ProteinHelper::buildOutput(GroupPowerBigraph* proteinGraph,
    ProteinProbEstimator* estimator){
  // array containing the PPs (Posterior ProbabilkitieS)
  Array<double> pps = proteinGraph->probabilityR;
  assert(pps.size()!=0);
  // array containing the PEPs (Posterior Error Probabilities)
  Array<double> peps = Array<double>(pps.size());
  // arrays that (will) contain protein ids and corresponding indexes
  Array< Array<string> > protein_ids =  Array< Array<string> >(pps.size());
  Array<int> indices = pps.sort();
  // arrays that (will) contain estimated and empirical q-values
  Array<double> estimQvalues = Array<double>(pps.size());
  Array<double> empirQvalues = Array<double>(pps.size());
  double sumTargetPepSoFar(0);
  unsigned int targetProtSoFar(0), decoyProtSoFar(0),
      targetsAtThr1(0), targetsAtThr2(0), decoysAtThr2(0);
  bool wellFormed = true;

  // calculating protein_ids, peps, qvalues
  for (int k=0; k<pps.size(); k++) {
    protein_ids[k] = proteinGraph->groupProtNames[indices[k]];
    double pep = 1.0 - pps[k];
    if(pep < 0.0) peps[k] = 0.0;
    else peps[k] = pep;
    if(isnan(peps[k])|| isinf(peps[k])) wellFormed = false;
    // estimated q-value: average pep of the proteins scoring better than the
    // current one (only targets, current one included)
    unsigned int targetsAtCurrentPep = countTargets(protein_ids[k], estimator);
    sumTargetPepSoFar += (peps[k]*targetsAtCurrentPep);
    targetProtSoFar += targetsAtCurrentPep;
    unsigned int decoysAtCurrentPep = countDecoys(protein_ids[k], estimator);
    decoyProtSoFar += decoysAtCurrentPep;
    estimQvalues[k] = sumTargetPepSoFar/targetProtSoFar;
    if(isnan(estimQvalues[k])|| isinf(estimQvalues[k])) wellFormed = false;
    // empirical q-value: #decoys/#targets
    empirQvalues[k] = (double)decoyProtSoFar/targetProtSoFar;
    if(isnan(empirQvalues[k])|| isinf(empirQvalues[k])) wellFormed = false;
    // update protein count under q-value thresholds
    if(estimQvalues[k]<fidoOutput::threshold2){
      targetsAtThr2+=targetsAtCurrentPep;
      decoysAtThr2+=decoysAtCurrentPep;
      if(estimQvalues[k]<fidoOutput::threshold1) {
        targetsAtThr1+=targetsAtCurrentPep;
      }
    }
  }
  // appending proteins with pps=0 (peps=1)
  if(proteinGraph->severedProteins.size()!=0){
    peps.add(1);
    protein_ids.add(proteinGraph->severedProteins);
    unsigned int targetsAtCurrentPep =
        countTargets(proteinGraph->severedProteins, estimator);
    sumTargetPepSoFar += (1.0*targetsAtCurrentPep);
    targetProtSoFar += targetsAtCurrentPep;
    decoyProtSoFar += countDecoys(proteinGraph->severedProteins, estimator);
    estimQvalues.add(sumTargetPepSoFar/targetProtSoFar);
    empirQvalues.add((double)decoyProtSoFar/targetProtSoFar);
  }

  double pi_0 = 1.0;
  // pi_0 is set to the q-value of the targets with the highest highest q-value
  // and empirical q-values are multiplied by it
  if(ProteinProbEstimator::usePi0==true) {
    pi_0 = estimQvalues[estimQvalues.size()-1];
    for(int k=0; k<empirQvalues.size(); k++)
      empirQvalues[k] = pi_0 * empirQvalues[k];
  }
  return fidoOutput(peps, protein_ids, estimQvalues, empirQvalues,
      targetsAtThr1, targetsAtThr2, decoysAtThr2,
      targetProtSoFar, decoyProtSoFar,
      pi_0, estimator->alpha, estimator->beta, wellFormed);
}

/**
 * writes protein weights to file.
 *
 * @param proteinGraph proteins and associated probabilities to be outputted
 * @param fileName file that will store the outputted information
 */
void ProteinHelper::writeOutputToFile(fidoOutput output, string fileName) {
  ofstream of(fileName.c_str());
  int size = output.peps.size();
  of << "PEP\t\t" << "est qvalues\t" << "emp qvalues\t" << "proteins\n";
  for (int k=0; k<size; k++) {
    of << scientific << setprecision(7) << output.peps[k] << "\t"
        << scientific << setprecision(7) << output.estimQvalues[k]<< "\t"
        << scientific << setprecision(7) << output.empirQvalues[k]<< "\t"
        << scientific << setprecision(7) << output.protein_ids[k] << endl;
  }
  of.close();
}

/**
 * Helper function for ProteinProbEstimator::writeXML; looks for PSMs associated
 * with a protein. The protein graph is divided into subgraphs: go through each
 * of them sequentially.
 *
 * @param protein_id the protein to be located
 * @param os stream the associated peptides will be written to
 * @param proteinsToPeptides hash table containing the associations between
 * proteins and peptides as calculated by Percolator
 */
void ProteinHelper::writeXML_writeAssociatedPeptides(string& protein_id,
    ofstream& os, map<string, vector<ScoreHolder*> >& proteinsToPeptides){
  vector<ScoreHolder*>* peptides = &proteinsToPeptides.find(protein_id)->second;
  vector<ScoreHolder*>::iterator peptIt = peptides->begin();
  for(; peptIt<peptides->end(); peptIt++){
    string pept = (*peptIt)->pPSM->getPeptideNoResidues();
    os << "      <peptide_seq seq=\"" << pept << "\"/>"<<endl;
    // check that protein_id is there
    assert((*peptIt)->pPSM->proteinIds.find(protein_id)
        != (*peptIt)->pPSM->proteinIds.end());
  }
}

/**
 * populates a hash table that associates the name of a protein with the list
 * of unique peptides that were associated to it by percolator. Part of this
 * same information (less conveniently indexed) is stored and manipulated by
 * fido in the proteinGraph field
 *
 * @param fullset set of unique peptides with scores computed by Percolator
 * @param proteinsToPeptides hash table to be populated
 */
void ProteinHelper::populateProteinsToPeptidesTable(Scores* fullset,
    ProteinProbEstimator* thisEstimator){
  vector<ScoreHolder>::iterator peptIt = fullset->scores.begin();
  // for each peptide
  for(; peptIt < fullset->scores.end();peptIt++){
    set<string>::iterator protIt = peptIt->pPSM->proteinIds.begin();
    // for each protein
    for(; protIt != peptIt->pPSM->proteinIds.end(); protIt++){
      // look for it in the hash table
      map<string, vector<ScoreHolder*> >::iterator found =
          thisEstimator->proteinsToPeptides.find(*protIt);
      if(found == thisEstimator->proteinsToPeptides.end()){
        // if not found insert a new protein-peptide pair...
        vector<ScoreHolder*> peptides;
        peptides.push_back(&*peptIt);
        pair<string,vector<ScoreHolder*> > protPeptPair(*protIt, peptides);
        thisEstimator->proteinsToPeptides.insert(protPeptPair);
      } else {
        // ... otherwise update
        found->second.push_back(&*peptIt);
      }
    }

  }
}

double ProteinHelper::isDecoyProbability(string protein_id,
    ProteinProbEstimator* estimator){
  vector<ScoreHolder*>::iterator peptIt =
      estimator->proteinsToPeptides.find(protein_id)->second.begin();
  vector<ScoreHolder*>::iterator peptItEnd =
      estimator->proteinsToPeptides.find(protein_id)->second.end();
  unsigned int decoys = 0;
  unsigned int targets = 0;
  for(; peptIt < peptItEnd; peptIt++){
    // check whether the peptide is target or a decoy
    if((*peptIt)->label != 1) decoys++;
    else targets++;
  }
  return decoys/(decoys+targets);
}

/**
 * prints the total number of proteins (targets and decoys) that were outputed
 * by fido, the number of targets at certain q-value thresholds and similar
 * statistics
 */
void ProteinHelper::printStatistics(const fidoOutput& output){
  assert(output.peps.size() != 0);
  if(VERB<=1) return;
  cerr << "There were " << output.totTargets << " targets and "
      << output.totDecoys << " decoy proteins (" << fixed << setprecision(2) <<
      (double)output.totDecoys/(output.totTargets+output.totDecoys)*100 <<
      "% decoys)\n";
  cerr << output.targetsAtThr1 << "\t: proteins found at a q-value of "
      << output.threshold1 <<"\n";
  cerr << output.targetsAtThr2 << "\t: proteins found at a q-value of "
      << output.threshold2 <<"\n";
  cerr << "Used pi_0 = " << output.pi_0 << "\n";
}

///////////////////////////////////////////////////////////////////////////////
// GRID SEARCH
///////////////////////////////////////////////////////////////////////////////

/**
 * point in the space generated by the possible values of the parameters alpha
 * and beta. Each point has two coordinates, two lists of protein names (the
 * true positives and false negatives calculated by fido for the current choice
 * of alpha and beta) and the value of the objective function calculated in
 * that point
 */
struct GridPoint {
    GridPoint(){
      objectiveFnValue = -1;
    }
    GridPoint(double alpha_par, double beta_par){
      objectiveFnValue = -1;
      alpha = alpha_par;
      beta = beta_par;
      targets = Array<string>();
      decoys = Array<string>();
    }
    ~GridPoint(){
      ;
    }
    void calculateObjectiveFn(double lambda, ProteinProbEstimator*
        toBeTested, ostringstream& debug);
    bool operator <(const GridPoint& rhs) const {
      assert(objectiveFnValue!=-1);
      assert(rhs.objectiveFnValue!=-1);
      if(objectiveFnValue < rhs.objectiveFnValue)
        return true;
      else return false;
    }
    Array<string> targets;
    Array<string> decoys;
    double objectiveFnValue;
    double alpha;
    double beta;
};

/**
 * 2D grid to be searched
 */
struct Grid{
    /**
     * constructs a 2D grid given ranges for its two variables
     */
    Grid(const double& l_a, const double& u_a, const double& l_b, const
        double& u_b, const double& i_a, const double& i_b):
          lower_a(l_a), upper_a(u_a), lower_b(l_b), upper_b(u_b),
        current(NULL), debugInfo(ostringstream::in | ostringstream::out),
        incrementAlpha(i_a), incrementBeta(i_b) {
      // before the search starts, the best point seen so far is artificially
      // set to max infinity (since we are minimizing!)
      bestSoFar = new GridPoint();
      bestSoFar->objectiveFnValue = numeric_limits<double>::max();
    }
    /**
     * constructor that builds a Grid in the default range
     */
    Grid(): current(NULL), debugInfo(ostringstream::in | ostringstream::out),
        incrementAlpha(0.4), incrementBeta(0.4){
      //incrementAlpha(0.3), incrementBeta(0.3){
      //lower_a = 0.05, upper_a = 0.3;
      //lower_b = 0.005, upper_b = 0.03;
      lower_a = 0.03, upper_a = 0.2;
      lower_b = 0.003, upper_b = 0.02;
      bestSoFar = new GridPoint();
      bestSoFar->objectiveFnValue = numeric_limits<double>::max();
    }
    Grid(const Grid&);
    ~Grid(){
      // only bestSoFar needs to be deleted: current either points to the same
      // location (if the last point examined was also the best found so far)
      // or has been deleted (if the last point examined was less good than the
      // best so far)
      delete bestSoFar;
    }
    void limitSearch(const int& dimension, const double& value);
    void toCurrentPoint();
    void calculateObjectiveFn(ProteinProbEstimator* toBeTested);
    void updateBest();
    bool wasSuccessful();
    void setToBest(ProteinProbEstimator* toBeTested);
    double getLower_a() const;
    double getUpper_a() const;
    double getLower_b() const;
    double getUpper_b() const;
    double updateCurrent_a();
    double updateCurrent_b();
    void compareAgainstDefault(ProteinProbEstimator* toBeTested);
    double current_a;
    double current_b;
    static int alpha;
    static int beta;
    ostringstream debugInfo;
  private:
    const static double lambda = 0.15;
    const double incrementAlpha;
    const double incrementBeta;
    double lower_a;
    double upper_a;
    double lower_b;
    double upper_b;
    GridPoint* bestSoFar;
    GridPoint* current;
};

int Grid::alpha=0;
int Grid::beta=1;

/**
 * @param dimension variable to be set: values are Grid::alpha and Grid::beta
 * @param value value the dimension should be set to
 */
void Grid::limitSearch(const int& dimension, const double& value){
  if(dimension==alpha) lower_a = upper_a = value;
  else if(dimension==beta) lower_b = upper_b = value;
}

double Grid::getLower_a() const {
  if(ProteinProbEstimator::logScaleSearch){
    // starting loop: print b label(s)
    if(VERB > 1) {
      cerr << endl << "\t\t";
      for(double b=log(lower_b); b<=log(upper_b); b+=incrementBeta)
        cerr << "beta=" << scientific << setprecision(1) << exp(b) << "\t";
    }
    // return value
    return log(lower_a);
  } else{
    // starting loop: print b label(s)
    if(VERB > 1) {
      cerr << endl << "\t\t";
      for(double b=lower_b; b<=upper_b; b+=incrementBeta)
        cerr << "beta=" << scientific << setprecision(1) << b << "\t";
    }
    // return value
    return lower_a;
  }
}

double Grid::getUpper_a() const {
  if(ProteinProbEstimator::logScaleSearch)
    return log(upper_a);
  else
    return upper_a;
}

double Grid::getLower_b() const {
  if(ProteinProbEstimator::logScaleSearch) {
    // starting loop: print a label
    if(VERB > 1) {
      cerr << "\nalpha=" << scientific << std::setprecision(1) << exp(current_a);
    }
    // return value
    return log(lower_b);
  } else {
    // starting loop: print a label
    if(VERB > 1) {
      cerr << "\nalpha=" << scientific << std::setprecision(1) << current_a;
    }
    // return value
    return lower_b;
  }
}

double Grid::getUpper_b() const {
  if(ProteinProbEstimator::logScaleSearch)
    return log(upper_b);
  else
    return upper_b;
}

double Grid::updateCurrent_a(){
  current_a+=incrementAlpha;
}
double Grid::updateCurrent_b(){
  current_b+=incrementBeta;
}

/**
 * updates to position of the current point to the coordinates given by
 * current_a and current_b
 */
void Grid::toCurrentPoint(){
  if(ProteinProbEstimator::logScaleSearch)
    current = new GridPoint(exp(current_a),exp(current_b));
  else
    current = new GridPoint(current_a, current_b);
  if (VERB > 2){
    debugInfo << scientific << setprecision(3) <<
        current->alpha << "\t" << current->beta << "\t";
  }
}

/**
 * calculates the objective function value in the current point.
 */
void Grid::calculateObjectiveFn(ProteinProbEstimator* toBeTested){
  current->calculateObjectiveFn(lambda,toBeTested,debugInfo);
  if(VERB > 1) {
    if(isinf(current->objectiveFnValue))
      cerr << "\t+infinity";
    else if(current->objectiveFnValue > 0) {
      cerr << "\t" << scientific << setprecision(3) <<
          "+" << current->objectiveFnValue;
    }
    else {
      cerr << "\t" << scientific << setprecision(3) <<
          current->objectiveFnValue;
    }
  }
}

/**
 * evaluate the objective function in the default location to compare
 * performances
 */
void Grid::compareAgainstDefault(ProteinProbEstimator* toBeTested){
  if(ProteinProbEstimator::logScaleSearch){
    current_a = log(0.1);
    current_b = log(0.01);
  } else {
    current_a = 0.1;
    current_b = 0.01;
  }
  toCurrentPoint();
  current->calculateObjectiveFn(lambda, toBeTested,debugInfo);
  updateBest();
}

/**
 * if the current point is lower than the bestSoFar, keep current (i.e. delete
 * bestSoFar and move pointer to current) else keep best so far (i.e. delete
 * current and leave bestSoFar pointing to the minimum)
 */
void Grid::updateBest(){
  if(*current < *bestSoFar){
    delete bestSoFar;
    bestSoFar = current;
  } else {
    delete current;
  }
}

/**
 * @return true if at least one meaningful pair of parameters has been found
 */
bool Grid::wasSuccessful(){
  if(bestSoFar->objectiveFnValue == numeric_limits<double>::max())
    return false;
  else return true;
}

/**
 * @param toBeTested protein estimator whose parameters were being grid
 * searched: they are set to the best found.
 */
void Grid::setToBest(ProteinProbEstimator* toBeTested){
  toBeTested->alpha = bestSoFar->alpha;
  toBeTested->beta = bestSoFar->beta;
}

/**
 * Helper function to GridPoint::calculateObjectiveFn: it populates the list of
 * true positive and false negative proteins by looking at the label (decoy or
 * target) of the peptides associated with the protein in the proteinsToPeptides
 * hash table
 */
void populateTargetDecoyLists(GridPoint* point, const fidoOutput& output,
    ProteinProbEstimator* toBeTested){
  assert(point->targets.size()==0);
  assert(point->decoys.size()==0);
  // for each protein in fido's output
  for(int i1=0; i1<output.peps.size(); i1++){
    for(int i2=0; i2<output.protein_ids[i1].size(); i2++){
      // store the protein id
      string protein_id = output.protein_ids[i1][i2];
      double probOfDecoy =
          ProteinHelper::isDecoyProbability(protein_id, toBeTested);
      if(probOfDecoy == 0) point->targets.add(protein_id);
      else if(probOfDecoy == 1) point->decoys.add(protein_id);
      else {
        point->targets.add(protein_id);
        point->decoys.add(protein_id);
      }

    }
  }
}

// forward declarations needed by gridPoint::calculateObjectiveFn

void calculateFDRs(
    const fidoOutput output,
    const Array<string>& targets, const Array<string>& decoys,
    Array<double>& estimatedFdrs, Array<double>& empiricalFdrs);

double calculateMSE_FDR(double threshold,
    const Array<double>& estimatedFdr, const Array<double>& empiricalFdr);

void calculateRoc(const fidoOutput output,
    const Array<string>& targets, const Array<string>& decoys,
    Array<int>& fps, Array<int>& tps);

double calculateROCX(int N, const Array<int>& fps, const Array<int>& tps);


#include "ProteinProbEstimatorDebugger.h"
/**
 * for a given choice of alpha and beta, calculates (1 − λ) MSE_FDR − λ ROCX
 * and stores the result in objectiveFnValue
 *
 * @return expression representing the objective function that was evaluated
 */
void GridPoint::calculateObjectiveFn(double lambda,
    ProteinProbEstimator* toBeTested, ostringstream& debug){
  assert(alpha!=-1);
  assert(beta!=-1);
  toBeTested->alpha = alpha;
  toBeTested->beta = beta;
  // we are in the middle of a grid search: don't want to recursively start one!
  bool startGridSearch = false;
  fidoOutput output = toBeTested->run(startGridSearch);
  if(output.wellFormed == false) {
    // if output in broken there is no need to continue: set the result to inf
    // this might happen for extremely low values of alpha and beta
    objectiveFnValue = numeric_limits<double>::infinity();
    debug << 1-lambda << "* ?" << " - " << lambda << "* ?" << endl;
    return;
  }
  populateTargetDecoyLists(this, output, toBeTested);
  if(ProteinProbEstimator::debugginMode) {
    // output the results of the probability calculation to file
    ProteinHelper::writeOutputToFile(output, string(WRITABLE_DIR) +
        "2percolator_fido_output.txt");
    string s = string(WRITABLE_DIR) + "3percolator_TPFP_lists.txt";
    ofstream o(s.c_str());
    o << targets << decoys << endl;
    o.close();
  }

  // calculate MSE_FDR
  Array<double> estimatedFdrs = Array<double>();
  Array<double> empiricalFdrs = Array<double>();
  calculateFDRs(output, targets, decoys,
      estimatedFdrs, empiricalFdrs);
  if(ProteinProbEstimator::debugginMode) {
    // output the results of the MSE_FDR to file
    string s = string(WRITABLE_DIR) + "4percolator_FDR_lists.txt";
    ofstream o1(s.c_str());
    o1 << "estimatedFdrs\n" <<estimatedFdrs << endl
        << "empiricalFdrs\n" << empiricalFdrs << endl;
    o1.close();
  }
  double threshold = 0.1;
  double mse_fdr = calculateMSE_FDR(threshold, estimatedFdrs, empiricalFdrs);
  if(isinf(mse_fdr)) {
    // if MSE is infinity abort: set the result to inf
    objectiveFnValue = numeric_limits<double>::infinity();
    debug << scientific << setprecision(3) <<
        "infinity" << "\t" << "   ?   " << "\t" <<
        objectiveFnValue << "\t" << "\n";
    return;
  }

  // calculate ROCX
  Array<int> fps = Array<int>();
  Array<int> tps = Array<int>();
  calculateRoc(output, targets, decoys, fps, tps);
  if(ProteinProbEstimator::debugginMode) {
    // output the results of the ROCX calculation to file
    string s = string(WRITABLE_DIR) + "5percolator_ROCX_lists.txt";
    ofstream o2(s.c_str());
    o2 << fps << endl << tps << endl;
    o2.close();
  }
  int N = output.decoysAtThr2;
  if(N==0) {
    // if the number of decoys at 5% is zero abort: set the result to inf
    objectiveFnValue = numeric_limits<double>::infinity();
    debug << scientific << setprecision(3) <<
        mse_fdr << "\t" << "-infinity" << "\t" <<
        "+infinity" << "\t" << "\n";
    return;
  }
  double rocX = calculateROCX(N, fps, tps);

  // combine results to calculate objective function value
  //objectiveFnValue = (1-lambda)*mse_fdr - lambda*rocX;
  objectiveFnValue = (1-lambda)*mse_fdr - lambda*rocX;
  debug << scientific << setprecision(3) <<
      mse_fdr << "\t" << rocX << "\t" <<
      objectiveFnValue << "\t" <<
      fixed << output.targetsAtThr1 << "\t\t" << output.targetsAtThr2 <<
      "\n";

  // debugging plots
  if(ProteinProbEstimator::gridSearchDebugPlotting){
    ProteinDebugger::plotQValues(output,toBeTested);
    ProteinDebugger::plotRoc(output,toBeTested,N);
  }
}

/**
 * quantifies the overlap between positiveNames and atThreshold (used to count
 * target and decoy proteins)
 */
int matchCount(const boost::unordered_set<string>& positiveNames,
    const Array<string> & atThreshold) {
  int count = 0;
  // for each protein
  for (int k=0; k<atThreshold.size(); k++) {
    // if target...
    if (positiveNames.find(atThreshold[k]) != positiveNames.end()){
      // ... done! if ties are being counted as one protein
      if(ProteinProbEstimator::tiesAsOneProtein){
        return 1;
      }
      // ... keep counting otherwise.
      else count++;
    }
  }
  return count;
}

Array<string> matches(const boost::unordered_set<string>& positiveNames,
    const Array<string> & atThreshold) {
  Array<string> result;
  for (int k=0; k<atThreshold.size(); k++) {
    if (positiveNames.find(atThreshold[k]) != positiveNames.end())
      result.add( atThreshold[k] );
  }
  return result;
}

/**
 * calculates empirical and estimated FDRs and stores the results in the
 * estimatedFdr and empiricalFdr Arrays for use in calculateMSE_FDR()
 */
void calculateFDRs(
    const fidoOutput output,
    const Array<string>& targets, const Array<string>& decoys,
    Array<double>& estimatedFdrs, Array<double>& empiricalFdrs) {

  estimatedFdrs.clear();
  empiricalFdrs.clear();
  boost::unordered_set<string> targetsSet(targets.size()),
      decoysSet(decoys.size());
  int k;
  for (k=0; k<targets.size(); k++)
    targetsSet.insert(targets[k]);
  for (k=0; k<decoys.size(); k++)
    decoysSet.insert(decoys[k]);
  Array<string> protsAtThreshold;
  string line;
  double prob, lastProb=-1;
  int decoysCount = 0, targetsCount = 0;
  int numScored = 0;
  Array<string> observedProteins;
  double estFDR = 0.0;
  double empiricalFDR = 0.0;
  double totalFDR = 0.0;
  bool scheduledUpdate = false;

  for(k=0; k<output.peps.size(); k++){
    prob = 1 - output.peps[k];
    protsAtThreshold = output.protein_ids[k];
    numScored += protsAtThreshold.size();
    observedProteins.append(protsAtThreshold);
    int decoysChange = matchCount(decoysSet, protsAtThreshold);
    int targetsChange = matchCount(targetsSet, protsAtThreshold);

    if ( prob != lastProb && lastProb != -1 ){
      scheduledUpdate = true;
    }
    if ( scheduledUpdate ) {
      if ( decoysChange > 0 || targetsChange > 0) {
        estimatedFdrs.add(estFDR);
        empiricalFdrs.add(empiricalFDR);
        scheduledUpdate = false;
      }
    }

    decoysCount += decoysChange;
    targetsCount += targetsChange;
    estFDR = output.estimQvalues[k];
    // calculating FDRs using a pi_0 approximation
    if(ProteinProbEstimator::usePi0){
      empiricalFDR = output.pi_0 * decoysCount/targetsCount;
    } else {
      empiricalFDR = (double)decoysCount/targetsCount;
    }
    double stored = output.empirQvalues[k];
    //double inspectMe = empiricalFDR;
    assert(abs(empiricalFDR-stored)<1e-10);
    /* the same done without pi_0
    totalFDR += (1-prob) * (fpChange + tpChange);
    estFDR = totalFDR / (fpCount + tpCount);
    empiricalFDR = double(fpCount) / (fpCount + tpCount);
     */
    lastProb = prob;
  }
  lastProb = prob;
  {
    estimatedFdrs.add(estFDR);
    empiricalFdrs.add(empiricalFDR);
  }
}

double squareAntiderivativeAt(double m, double b, double xVal) {
  double u = m*m;
  double v = 2*m*b;
  double t = b*b;
  return u*xVal*xVal*xVal/3.0 + v*xVal*xVal/2.0 + t*xVal;
}

double antiderivativeAt(double m, double b, double xVal) {
  return m*xVal*xVal/2.0 + b*xVal;
}

double areaSq(double x1, double y1, double x2, double y2, double threshold) {
  double m = (y2-y1)/(x2-x1);
  double b = y1-m*x1;
  double area = squareAntiderivativeAt(m, b, min(threshold, x2) ) -
      squareAntiderivativeAt(m, b, x1);
  if(isnan(area)) return 0.0;
  else return area;
}

/**
 * calculates the FDR Mean Square Error
 *
 * @return FDR Mean Square Error
 */
double calculateMSE_FDR(double threshold,
    const Array<double>& estimatedFdr, const Array<double>& empiricalFdr) {
  assert(estimatedFdr.size() == empiricalFdr.size());
  Vector diff = Vector(estimatedFdr) - Vector(empiricalFdr);
  double tot = 0.0;
  int k;
  for (k=0; k<diff.size()-1; k++) {
    // stop if no part of the estFDR is < threshold
    // values are monotonically increasing
    if (estimatedFdr[k] >= threshold) {
      if (k == 0)
        tot = 1.0 / 0.0;
      break;
    }
    tot += areaSq(estimatedFdr[k],diff[k],estimatedFdr[k+1],diff[k+1],threshold);
  }
  double xRange = min(threshold, estimatedFdr[k]) - estimatedFdr[0];

  if (isinf(tot))
    return tot;
  else return (tot/xRange);
}

/**
 * calculates the roc curve and stores the results in the fps and tps Arrays
 * for use in calculateROCX
 *
 */
void calculateRoc(const fidoOutput output,
    const Array<string>& targets, const Array<string>& decoys,
    Array<int>& fps, Array<int>& tps) {

  boost::unordered_set<string> targetsSet(targets.size()),
      decoysSet(decoys.size());
  int k;
  for (k=0; k<targets.size(); k++) {
    targetsSet.insert( targets[k] );
  }
  for (k=0; k<decoys.size(); k++) {
    decoysSet.insert( decoys[k] );
  }
  Array<string> protsAtThreshold;
  string line;
  double prob, lastProb=-1;
  int decoysCount = 0, targetsCount = 0;
  int numScored = 0;
  Array<string> observedProteins;
  fps.add(0);
  tps.add(0);
  bool scheduledUpdate = false;
  double totalFDR = 0.0, estFDR = 0.0;

  for(k=0; k<output.peps.size(); k++) {
    prob = 1 - output.peps[k];
    protsAtThreshold = output.protein_ids[k];
    numScored += protsAtThreshold.size();
    observedProteins.append( protsAtThreshold );
    int decoysChange = matchCount(decoysSet, protsAtThreshold);
    int targetsChange = matchCount(targetsSet, protsAtThreshold);
    if ( prob != lastProb && lastProb != -1 ) {
      scheduledUpdate = true;
    }
    if ( scheduledUpdate ) {
      fps.add(decoysCount);
      tps.add(targetsCount);
      scheduledUpdate = false;
      // calculating FDRs using a pi_0 approximation
      estFDR = output.estimQvalues[k];
      /* the same done without pi_0
      totalFDR += (1-prob) * (fpChange + tpChange);
      estFDR = totalFDR / (fpCount + tpCount);
       */
    }
    decoysCount += decoysChange;
    targetsCount += targetsChange;
    lastProb = prob;
  }
  lastProb = prob;
  fps.add( decoysCount );
  tps.add( targetsCount );
  fps.add( decoysSet.size() );
  tps.add( targetsSet.size() );
}

double area(double x1, double y1, double x2, double y2, int N)
{
  double m = (y2-y1)/(x2-x1);
  double b = y1-m*x1;
  double area = antiderivativeAt(m, b, min(double(N), x2) ) -
      antiderivativeAt(m, b, x1);
  if(isnan(area)) return 0.0;
  else return area;
}

/**
 * calculates the area under the roc curve up to N false positives
 *
 * @return rocX
 */
double calculateROCX(int N, const Array<int>& fps, const Array<int>& tps){
  double rocN = 0.0;
  if ( fps.back() < N ) {
    cerr << "There are not enough false positives; needed " << N
        << " and was only given " << fps.back() << endl << endl;
    exit(1);
  }
  for (int k=0; k<fps.size()-1; k++) {
    // find segments where the fp value changes
    if ( fps[k] >= N )
      break;
    if ( fps[k] != fps[k+1] ) {
      // this line segment is a function
      double currentArea = area(fps[k], tps[k], fps[k+1], tps[k+1], N);
      rocN += currentArea;
    }
  }
  double rocX = rocN / (N * tps.back());
  return rocX;
}

#endif /* PROTEINPROBESTIMATORHELPER_H_ */
