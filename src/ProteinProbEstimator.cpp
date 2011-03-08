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

#include "ProteinProbEstimator.h"
#include <iostream>
#include <fstream>

ProteinProbEstimator* ProteinProbEstimator::instance=0;


ProteinProbEstimator* ProteinProbEstimator::getInstance(){
  if(instance==0)
    instance = new ProteinProbEstimator();
  return instance;
}

/**
 * writes protein weights to file.
 *
 * @param proteinGraph proteins and associated probabilities to be outputted
 * @param fileName file that will store the outputted information
 */
void writeProteinProbToFile(GroupPowerBigraph* proteinGraph, string fileName) {
  ofstream proteinProbToFile;
  proteinProbToFile.open(fileName.c_str());
  Array<double> sorted = proteinGraph->probabilityR;
  Array<int> indices = sorted.sort();
  cerr << "\nWriting Protein level probabilities to file:";
  cerr << fileName << endl;
  for (int k=0; k<sorted.size(); k++) {
    proteinProbToFile << sorted[k] << " "
        << proteinGraph->groupProtNames[ indices[k] ] << endl;
  }
  proteinProbToFile.close();
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
void populateProteinsToPeptides(Scores* fullset,
    map<string, vector<ScoreHolder*> >& proteinsToPeptides){
  vector<ScoreHolder>::iterator peptIt = fullset->scores.begin();
  // for each peptide
  for(; peptIt < fullset->scores.end();peptIt++){
    set<string>::iterator protIt = peptIt->pPSM->proteinIds.begin();
    // for each protein
    for(; protIt != peptIt->pPSM->proteinIds.end(); protIt++){
      // look for it in the hash table
      map<string, vector<ScoreHolder*> >::iterator found =
          proteinsToPeptides.find(*protIt);
      if(found == proteinsToPeptides.end()){
        // if not found insert a new protein-peptide pair...
        vector<ScoreHolder*> peptides;
        peptides.push_back(&*peptIt);
        pair<string,vector<ScoreHolder*> > protPeptPair(*protIt, peptides);
        proteinsToPeptides.insert(protPeptPair);
      } else {
        // ... otherwise update
        found->second.push_back(&*peptIt);
      }
    }

  }
}

/**
 * choose the set of parameters that jointly maximizes the ROC50 score (the
 * average sensitivity when allowing between zero and 50 false positives) and
 * minimizes the mean squared error (MSE) from an ideally calibrated
 * probability.
 */
void ProteinProbEstimator::gridSearchAlphaBeta(){
  //TODO
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
int ProteinProbEstimator::calculateProteinProb(Scores* fullset,
    bool gridSearch){
  populateProteinsToPeptides(fullset, proteinsToPeptides);
  srand(time(NULL));
  cout.precision(8);
  cerr.precision(8);

  // run grid search for <alpha,beta> values by default
  if(gridSearch) gridSearchAlphaBeta();

  // if not already set, set parameters values to default
  if(alpha == -1) alpha = 0.1;
  if(beta == -1) beta = 0.01;

  //GroupPowerBigraph::LOG_MAX_ALLOWED_CONFIGURATIONS = ;
  proteinGraph = new GroupPowerBigraph ( RealRange(alpha, 1, alpha),
      RealRange(beta, 1, beta), gamma );
  proteinGraph->read(fullset);
  proteinGraph->getProteinProbs();

  // print the protein weights to file
  writeProteinProbToFile(proteinGraph, "/tmp/fido/fidoOut.txt");

  return 0;
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
    if ((*peptIt)->pPSM->proteinIds.find(protein_id) == (*peptIt)->pPSM->proteinIds.end()){
      cout << "MERDA" << endl;
      exit(-1);
    }
  }
  /*
  // the following code was trying to achieve the same as the above exclusively
  // with information coming from fido (without making use of external
  // information coming from Percolator
  bool found = false;
  int peptidesIndex = -1;
  int proteinIndex = 0;
  // look for protein_id in the subgraphs
  while(peptidesIndex == -1 && proteinIndex<proteinGraph->subgraphs.size()){
    StringTable proteins;
    proteins = StringTable::AddElements(
        proteinGraph->subgraphs[proteinIndex].proteinsToPSMs.names);
    peptidesIndex = proteins.lookup(protein_id);
    if(peptidesIndex == -1) proteinIndex++;
  }
  // if it was found
  if(peptidesIndex!=-1){
    Set peptides = proteinGraph->subgraphs[proteinIndex].
        proteinsToPSMs.associations[peptidesIndex];
    // for each PSM associated with the protein, print the peptides
    for (int k=0; k<peptides.size(); k++) {
      string pept = proteinGraph->subgraphs[proteinIndex].
          PSMsToProteins.names[peptides[k]];
      os << "      <peptide_seq seq=\"" << pept << "\"/>"<<endl;
    }
  } else {
    // it seems to me that every protein should have associated peptides
    // as they were given in input to calculateProteinProb. At present this
    // is not so. Consider throwing an error here.
  }
  return;
  */
}

/**
 * output protein level probabilites results in xml format
 *
 * @param os stream to which the xml is directed
 */
void ProteinProbEstimator::writeXML(ofstream& os){
  assert(proteinGraph!=0);
  // print protein weights to cerr
  // proteinGraph->printProteinWeights();
  os << "  <proteins>" << endl;
  Array<double> probabilities = proteinGraph->probabilityR;
  Array<int> indices = probabilities.sort();

  // for each probability
  for (int k=0; k<probabilities.size(); k++) {
    // for each protein with a certain probability
    for(int k2=0; k2<proteinGraph->groupProtNames[indices[k]].size(); k2++) {
      string protein_id = proteinGraph->groupProtNames[indices[k]][k2];
      os << "    <protein p:protein_id=\"" << protein_id << "\">" << endl;
      double pep = probabilities[k];
      os << "      <pep>" << pep << "</pep>" << endl;
      os << "      <q_value>" << 0.0 << "</q_value>" << endl;
      writeXML_writeAssociatedPeptides(protein_id, os, proteinsToPeptides);
      os << "    </protein>" << endl;
    }
  }
  for(int k=0; k<proteinGraph->severedProteins.size(); k++) {
    string protein_id = proteinGraph->severedProteins[k];
    os << "    <protein p:protein_id=\"" << protein_id << "\">" << endl;
    os << "      <pep>" << 0.0 << "</pep>" << endl;
    os << "      <q_value>" << 0.0 << "</q_value>" << endl;
    writeXML_writeAssociatedPeptides(protein_id, os, proteinsToPeptides);
    os << "    </protein>" << endl;
  }
  os << "  </proteins>" << endl;
}
