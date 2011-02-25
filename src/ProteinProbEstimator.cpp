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
      // look for PSMs associated with the protein. The protein graph is devided into subgraphs:
      // go through each of them sequentially
      bool foundIt = false;
      for(int k3=0; k3<proteinGraph->subgraphs.size() && !foundIt; k3++){
        StringTable proteins;
        proteins = StringTable::AddElements(proteinGraph->subgraphs[k3].proteinsToPSMs.names);
        int peptidesIndex = proteins.lookup(protein_id);
        // if the protein has been found in the present subgraph
        if(peptidesIndex!=-1){
          foundIt = true;
          Set peptides = proteinGraph->subgraphs[k3].proteinsToPSMs.associations[peptidesIndex];
          int kkk = peptides.size();
          // for each PSM associated with the protein, print the peptides
          for (int k4=0; k4<peptides.size(); k4++) {
            string prot = protein_id;
            string pept =  proteinGraph->subgraphs[k3].PSMsToProteins.names[peptides[k4]];
            os << "      <peptide_seq seq=\"" << pept << "\"/>"<<endl;
          }
        }
      }
      os << "    </protein>" << endl;
    }
  }
  for(int k=0; k<proteinGraph->severedProteins.size(); k++) {
    string protein_id = proteinGraph->severedProteins[k];
    os << "    <protein p:protein_id=\"" << protein_id << "\">" << endl;
    os << "      <pep>" << 0.0 << "</pep>" << endl;
    os << "      <q_value>" << 0.0 << "</q_value>" << endl;
    //os << "      <peptide_seq
    os << "    </protein>" << endl;
  }
  os << "  </proteins>" << endl;
}
