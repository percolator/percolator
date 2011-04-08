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
#include "Vector.h"
#include "Globals.h"
using namespace std;
#include <ext/hash_set>


///////////////////////////////////////////////////////////////////////////////
// HELPER FUNCTIONS FOR ProteinProbEstimator METHODS
///////////////////////////////////////////////////////////////////////////////

double isDecoyProbability(string protein_id, ProteinProbEstimator* estimator);

unsigned int countTargets(const Array<string>& protein_ids,
    ProteinProbEstimator* estimator){
  unsigned int count(0);
  for(int p=0; p<protein_ids.size(); p++)
    if(isDecoyProbability(protein_ids[p], estimator)<0.5) count++;
  return count;
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
fidoOutput buildOutput(GroupPowerBigraph* proteinGraph,
    ProteinProbEstimator* estimator){
  // array containing the PPs (Posterior ProbabilkitieS)
  Array<double> pps = proteinGraph->probabilityR;
  assert(pps.size()!=0);
  // array containing the PEPs (Posterior Error ProbabilkitieS)
  Array<double> peps = Array<double>(pps.size());
  // arrays that (will) contain protein ids and corresponding indexes
  Array< Array<string> > protein_ids =  Array< Array<string> >(pps.size());
  Array<int> indices = pps.sort();
  // array that (will) contain q-values
  Array<double> qvalues = Array<double>(pps.size());
  double sumTargetPepSoFar(0);
  unsigned int targetProtSoFar(0), proteinsAtThr1(0), proteinsAtThr2(0);

  // calculating protein_ids, peps, qvalues
  for (int k=0; k<pps.size(); k++) {
    protein_ids[k] = proteinGraph->groupProtNames[indices[k]];
    peps[k] = 1.0 - pps[k];
    // q-value: average pep of the proteins scoring better than the current one
    // (only targets, current one included)
    unsigned int targetsAtCurrentPep = countTargets(protein_ids[k], estimator);
    sumTargetPepSoFar += (peps[k]*targetsAtCurrentPep);
    targetProtSoFar += targetsAtCurrentPep;
    double qValue = sumTargetPepSoFar/targetProtSoFar;
    qvalues[k] = qValue;
    // update protein count under q-value thresholds
    if(qValue<fidoOutput::threshold2){ proteinsAtThr2+=targetsAtCurrentPep;
    if(qValue<fidoOutput::threshold1) proteinsAtThr1+=targetsAtCurrentPep;
    }
  }
  // appending proteins with peps=0 (pps=1)
  if(proteinGraph->severedProteins.size()!=0){
    peps.add(1);
    protein_ids.add(proteinGraph->severedProteins);
    unsigned int targetsAtCurrentPep =
        countTargets(proteinGraph->severedProteins, estimator);
    sumTargetPepSoFar += (1*targetsAtCurrentPep);
    targetProtSoFar += targetsAtCurrentPep;
    qvalues.add(sumTargetPepSoFar/targetProtSoFar);
  }

  double pi_0 = qvalues[pps.size()-1];
  return fidoOutput(peps, protein_ids, qvalues, proteinsAtThr1, proteinsAtThr2,
      targetProtSoFar, pi_0);
}

/**
 * writes protein weights to file.
 *
 * @param proteinGraph proteins and associated probabilities to be outputted
 * @param fileName file that will store the outputted information
 */
void writeOutputToFile(fidoOutput output, string fileName) {
  ofstream of(fileName.c_str());
  int size = output.peps.size();
  of << "PEP\t\t" << "q-value\t\t" << "proteins\n";
  for (int k=0; k<size; k++) {
    of << scientific << setprecision(7) << output.peps[k] << "\t"
        << scientific << setprecision(7) << output.qvalues[k]<< "\t"
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
void writeXML_writeAssociatedPeptides(string& protein_id,
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
void populateProteinsToPeptidesTable(Scores* fullset,
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

double isDecoyProbability(string protein_id, ProteinProbEstimator* estimator){
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

///////////////////////////////////////////////////////////////////////////////
// GRID SEARCH
///////////////////////////////////////////////////////////////////////////////

/**
 * point in the space generated by the possible values of the parameters alpha
 * and beta. Each point has two coordinates, two lists of protein names (the
 * true positives and false negatives calculated by fido with the particular
 * choice of alpha and beta) and the value of the objective function calculated
 * in that point
 */
struct GridPoint {
    GridPoint(){
      objectiveFnValue = -1;
    }
    GridPoint(double alpha_par, double beta_par){
      objectiveFnValue = -1;
      alpha = alpha_par;
      beta = beta_par;
      truePositives = Array<string>();
      falsePositives = Array<string>();
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
    Array<string> truePositives;
    Array<string> falsePositives;
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
        double& u_b): lower_a(l_a), upper_a(u_a), lower_b(l_b), upper_b(u_b),
        current(NULL), debugInfo(ostringstream::in | ostringstream::out) {
      // before the search starts, the best point seen so far is artificially
      // set to max infinity (since we are minimizing!)
      bestSoFar = new GridPoint();
      bestSoFar->objectiveFnValue = numeric_limits<double>::max();
    }
    /**
     * constructor that builds a Grid in the default range
     */
    Grid(): current(NULL), debugInfo(ostringstream::in | ostringstream::out) {
      lower_a = 0.07, upper_a = 0.3;
      lower_b = 0.007, upper_b = 0.1;
      bestSoFar = new GridPoint();
      bestSoFar->objectiveFnValue = numeric_limits<double>::max();
    }
    Grid(const Grid&);
    ~Grid(){
      delete bestSoFar;
      if(current) delete current;
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
    static void testGridRanges();
    void compareAgainstDefault(ProteinProbEstimator* toBeTested);
    double current_a;
    double current_b;
    static int alpha;
    static int beta;
    ostringstream debugInfo;
  private:
    const static double lambda = 0.25;
    const static double incrementAlpha=0.3;
    const static double incrementBeta=0.3;
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
  // starting loop: print b label(s)
  if(VERB > 1) {
    cerr << endl << "\t\t";
    for(double b=log(lower_b); b<=log(upper_b); b+=incrementBeta)
      cerr << "beta=" << scientific << setprecision(1) << exp(b) << "\t";
  }
  // return value
  return log(lower_a);
}
double Grid::getUpper_a() const {
  return log(upper_a);
}

double Grid::getLower_b() const {
  // starting loop: print a label
  if(VERB > 1) {
    cerr << "\nalpha=" << scientific << std::setprecision(1) << exp(current_a);
  }
  // return value
  return log(lower_b);
}
double Grid::getUpper_b() const {
  return log(upper_b);
}

double Grid::updateCurrent_a(){
  current_a+=incrementAlpha;
}
double Grid::updateCurrent_b(){
  current_b+=incrementBeta;
}

/** prints the grid
 */
void Grid::testGridRanges(){
  Grid grid = Grid();
  grid.current_a = grid.getLower_a();
  for(; grid.current_a<=grid.getUpper_a(); grid.updateCurrent_a()){
    grid.current_b = grid.getLower_b();
    for(; grid.current_b<=grid.getUpper_b(); grid.updateCurrent_b()){
    }
  }
}

/**
 * updates to position of the current point to the coordinates given by
 * current_a and current_b
 */
void Grid::toCurrentPoint(){
  current = new GridPoint(exp(current_a),exp(current_b));
  if (VERB > 2){
    debugInfo << "alpha=" << scientific << setprecision(1) << current->alpha
        << ", " << "beta=" <<scientific << setprecision(1) << current->beta
        << ": ";
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
      cerr << "\t" << fixed << std::setprecision(5)
      << "+" << current->objectiveFnValue;
    }
    else {
      cerr << "\t" << fixed << std::setprecision(5)
      << current->objectiveFnValue;
    }
  }
}

/**
 * evaluate the objective function in the default location to compare
 * performances
 */
void Grid::compareAgainstDefault(ProteinProbEstimator* toBeTested){
  current_a = log(0.1);
  current_b = log(0.01);
  toCurrentPoint();
  current->calculateObjectiveFn(lambda, toBeTested,debugInfo);
  if(VERB > 2) {
    cerr << "Objective fn estimation for default values "
        <<"(alpha=0.1, beta=0.01): "
        << current->objectiveFnValue;
  }
  updateBest();
}

/**
 * if the current point is lower than the bestSoFar, update. (minimizing!)
 */
void Grid::updateBest(){
  if(*current < *bestSoFar) bestSoFar = current;
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
void populateTPandFNLists(GridPoint* point, const fidoOutput& output,
    ProteinProbEstimator* toBeTested){
  assert(point->truePositives.size()==0);
  assert(point->falsePositives.size()==0);
  // for each protein in fido's output
  for(int i1=0; i1<output.peps.size(); i1++){
    for(int i2=0; i2<output.protein_ids[i1].size(); i2++){
      // store the protein id
      string protein_id = output.protein_ids[i1][i2];
      double probOfDecoy = isDecoyProbability(protein_id, toBeTested);
      if(probOfDecoy == 0) point->truePositives.add(protein_id);
      else if(probOfDecoy == 1) point->falsePositives.add(protein_id);
      else {
        point->truePositives.add(protein_id);
        point->falsePositives.add(protein_id);
      }

    }
  }
}

// forward declarations needed by gridPoint::calculateObjectiveFn

void calculateFDRs(
    const fidoOutput output,
    const Array<string>& truePositives, const Array<string>& falsePositives,
    Array<double>& estimatedFdrs, Array<double>& empiricalFdrs);

double calculateMSE_FDR(double threshold,
    const Array<double>& estimatedFdr, const Array<double>& empiricalFdr);

void calculateRoc(const fidoOutput output,
    const Array<string>& truePositives, const Array<string>& falsePositives,
    Array<int>& fps, Array<int>& tps);

double calculateROC50(int N, const Array<int>& fps, const Array<int>& tps);

/**
 * for a given choice of alpha and beta, calculates (1 − λ) MSE_FDR − λ ROC50
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
  fidoOutput output = toBeTested->calculateProteinProb(false);
  populateTPandFNLists(this, output, toBeTested);
  if(VERB > 4) {
    // output the results of the probability calculation to file
    writeOutputToFile(output, "/tmp/2percolator_fido_output.txt");
    ofstream o("/tmp/3percolator_TPFP_lists.txt");
    o << truePositives << falsePositives << endl;
    o.close();
  }

  // calculate MSE_FDR
  Array<double> estimatedFdrs = Array<double>();
  Array<double> empiricalFdrs = Array<double>();
  calculateFDRs(output, truePositives, falsePositives,
      estimatedFdrs, empiricalFdrs);
  if(VERB > 4) {
    // output the results of the MSE_FDR to file
    ofstream o1 ("/tmp/4percolator_FDR_lists.txt");
    o1 << "estimatedFdrs\n" <<estimatedFdrs << endl
        << "empiricalFdrs\n" << empiricalFdrs << endl;
    o1.close();
  }
  double threshold = 0.1;
  double mse_fdr = calculateMSE_FDR(threshold, estimatedFdrs, empiricalFdrs);

  if(isinf(mse_fdr)) {
    // if MSE is infinity there is no need to continue: set the result to inf
    objectiveFnValue = mse_fdr;
    debug << 1-lambda << "* infinity" << " - " << lambda << "* ?" << endl;
  }

  // calculate ROC50
  Array<int> fps = Array<int>();
  Array<int> tps = Array<int>();
  calculateRoc(output, truePositives, falsePositives, fps, tps);
  if(VERB > 4) {
    // output the results of the ROC50 calculation to file
    ofstream o2("/tmp/5percolator_ROC50_lists.txt");
    o2 << fps << endl << tps << endl;
    o2.close();
  }
  int N = 50;
  double roc50 = calculateROC50(N, fps, tps);

  objectiveFnValue = (1-lambda)*mse_fdr - lambda*roc50;
  debug << 1-lambda << "*" << scientific << setprecision(1) << mse_fdr << " - "
      << lambda << "*" << scientific << setprecision(1) << roc50 << " = "
      << scientific << setprecision(1) << objectiveFnValue <<"\n";
}

namespace __gnu_cxx {
template<> struct hash< std::string > {
    size_t operator()( const std::string & x ) const {
      return hash< const char* >()( x.c_str() );
    }
};
}

int matchCount(const __gnu_cxx::hash_set<string> & positiveNames,
    const Array<string> & atThreshold) {
  int count = 0;
  for (int k=0; k<atThreshold.size(); k++) {
    if ( positiveNames.count( atThreshold[k] ) > 0 )
      count++;
  }
  return count;
}

Array<string> matches(const __gnu_cxx::hash_set<string> & positiveNames,
    const Array<string> & atThreshold) {
  Array<string> result;
  for (int k=0; k<atThreshold.size(); k++) {
    if ( positiveNames.count( atThreshold[k] ) > 0 )
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
    const Array<string>& truePositives, const Array<string>& falsePositives,
    Array<double>& estimatedFdrs, Array<double>& empiricalFdrs) {

  estimatedFdrs.clear();
  empiricalFdrs.clear();
  __gnu_cxx::hash_set<string> truePosSet(truePositives.size()),
      falsePosSet(falsePositives.size());
  int k;
  for (k=0; k<truePositives.size(); k++)
    truePosSet.insert(truePositives[k]);
  for (k=0; k<falsePositives.size(); k++)
    falsePosSet.insert(falsePositives[k]);
  Array<string> protsAtThreshold;
  string line;
  double prob, lastProb=-1;
  int fpCount = 0, tpCount = 0;
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
    int fpChange = matchCount(falsePosSet, protsAtThreshold);
    int tpChange = matchCount(truePosSet, protsAtThreshold);

    if ( prob != lastProb && lastProb != -1 ){
      scheduledUpdate = true;
    }
    if ( scheduledUpdate ) {
      if ( fpChange > 0 || tpChange > 0) {
        estimatedFdrs.add(estFDR);
        empiricalFdrs.add(empiricalFDR);
        scheduledUpdate = false;
      }
    }

    fpCount += fpChange;
    tpCount += tpChange;
    // calculating FDRs using a pi_0 approximation
    estFDR = output.qvalues[k];
    empiricalFDR = output.pi_0 * fpCount/tpCount;
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
  double area = squareAntiderivativeAt(m, b, min(threshold, x2) )
                          - squareAntiderivativeAt(m, b, x1);
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
 * for use in calculateROC50
 *
 */
void calculateRoc(const fidoOutput output,
    const Array<string>& truePositives, const Array<string>& falsePositives,
    Array<int>& fps, Array<int>& tps) {

  __gnu_cxx::hash_set<string> truePosSet(truePositives.size()),
      falsePosSet(falsePositives.size());
  int k;
  for (k=0; k<truePositives.size(); k++) {
    truePosSet.insert( truePositives[k] );
  }
  for (k=0; k<falsePositives.size(); k++) {
    falsePosSet.insert( falsePositives[k] );
  }
  Array<string> protsAtThreshold;
  string line;
  double prob, lastProb=-1;
  int fpCount = 0, tpCount = 0;
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
    int fpChange = matchCount(falsePosSet, protsAtThreshold);
    int tpChange = matchCount(truePosSet, protsAtThreshold);
    if ( prob != lastProb && lastProb != -1 ) {
      scheduledUpdate = true;
    }
    if ( scheduledUpdate ) {
      fps.add( fpCount );
      tps.add( tpCount );
      scheduledUpdate = false;
      // calculating FDRs using a pi_0 approximation
      estFDR = output.qvalues[k];
      /* the same done without pi_0
      totalFDR += (1-prob) * (fpChange + tpChange);
      estFDR = totalFDR / (fpCount + tpCount);
       */
    }
    fpCount += fpChange;
    tpCount += tpChange;
    lastProb = prob;
  }
  lastProb = prob;
  fps.add( fpCount );
  tps.add( tpCount );
  fps.add( falsePosSet.size() );
  tps.add( truePosSet.size() );
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
 * calculates the area under the roc curve up to 50 false positives
 *
 * @return roc50
 */
double calculateROC50(int N, const Array<int>& fps, const Array<int>& tps){
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
  double roc50 = rocN / (N * tps.back());
  return roc50;
}

#endif /* PROTEINPROBESTIMATORHELPER_H_ */
