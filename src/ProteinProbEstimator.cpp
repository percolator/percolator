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

ProteinProbEstimator* ProteinProbEstimator::instance=0;


ProteinProbEstimator* ProteinProbEstimator::getInstance(){
  if(instance==0)
    instance = new ProteinProbEstimator();
  return instance;
}

/**
 * calculate protein level probabilities
 *
 * @param fullset set of unique peptides with scores computed by Percolator
 */
int ProteinProbEstimator::calculateProteinProb(Scores& fullset){
  srand(time(NULL));
  cout.precision(8);
  cerr.precision(8);

  double gamma = 0.5; //gridSearchGamma();
  double alpha = 0.1; //gridSearchAlpha();
  double beta = 0.01; //gridSearchBeta();

  //GroupPowerBigraph::LOG_MAX_ALLOWED_CONFIGURATIONS = ;
  proteinGraph = new GroupPowerBigraph ( RealRange(alpha, 1, alpha),
      RealRange(beta, 1, beta), gamma );
  proteinGraph->read(fullset);
  proteinGraph->getProteinProbs();
  return 0;
}

/**
 * Helper function for ProteinProbEstimator::writeXML; looks for PSMs associated
 * with a protein. The protein graph is devided into subgraphs: go through each
 * of them sequentially.
 *
 * @param protein the protein to be located
 */
void writeXML_writeAssociatedPeptides(GroupPowerBigraph* proteinGraph,
    string& protein_id, ofstream& os){
  bool found = false;
  int peptidesIndex = -1;
  int proteinIndex = 0 ;
  // look for protein_id in the subgraphs
  while(peptidesIndex == -1 && proteinIndex<proteinGraph->subgraphs.size()){
    StringTable proteins;
    proteins = StringTable::AddElements(
        proteinGraph->subgraphs[proteinIndex].proteinsToPSMs.names);
    peptidesIndex = proteins.lookup(protein_id);
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
      writeXML_writeAssociatedPeptides(proteinGraph, protein_id, os);
      os << "    </protein>" << endl;
    }
  }
  for(int k=0; k<proteinGraph->severedProteins.size(); k++) {
    string protein_id = proteinGraph->severedProteins[k];
    os << "    <protein p:protein_id=\"" << protein_id << "\">" << endl;
    os << "      <pep>" << 0.0 << "</pep>" << endl;
    os << "      <q_value>" << 0.0 << "</q_value>" << endl;
    writeXML_writeAssociatedPeptides(proteinGraph, protein_id, os);
    os << "    </protein>" << endl;
  }
  os << "  </proteins>" << endl;
}
