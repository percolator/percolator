/*******************************************************************************
 Copyright 2006-2012 Lukas Käll <lukas.kall@scilifelab.se>

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

#ifndef FIDO_MAIN_H_
#define FIDO_MAIN_H_

#include "FidoInterface.h"

using namespace std;



class Fido {
  
      /** FIDO PARAMETERS **/
    
    /** when grouping proteins discard all possible combinations for each group*/
    const static bool trivialGrouping = false;
    /** compute peptide level prior probability instead of using default = 0.1 **/
    const static bool computePriors = false;
    /** threshold used for fido to remove poor PSMs **/
    const static double psmThreshold = 0.0;
    const static double reduced_psmThreshold = 0.2;
    /** threshold used for fido to classify a peptide as very low confidence **/
    const static double peptideThreshold = 0.001;
    const static double reduced_peptideThreshold = 0.2;
    /** threshold used for fido to classify a protein as very low confidence and prune it **/
    const static double proteinThreshold = 0.01; 
    const static double reduced_proteinThreshold = 0.2;
    /** default value for peptide prior probability used in fido to compute the peptide likehood **/
    const static double peptidePrior = 0.1; 
    /** number of maximum of tree configurations allowed in fido **/
    const static double max_allow_configurations = 18;
    /** allow the presence of peptides with the same sequence but different label (target/decoy) **/
    const static bool allow_multiple_labeled_peptides = false;
    
    /** GRID SEARCH PARAMETERS **/
    
    /** value that balances the objective function equation (lambda * rocR) - (1-lambda) * (fdr_mse) **/
    const static double lambda = 0.15;
    
    /** number of false positives allowed while estiaming the ROC curve score **/
    /** if updateRocN is true the N value will be estimated automatically according to the number of FP found at a certain threshold **/
    const static unsigned default_rocN = 50;
    const static bool updateRocN = true;
    
    /** activate the optimization of the parameters to see the best boundaries**/
    const static bool optimize = false;

 public:
   
	Fido();
	virtual ~Fido();
	string greeter();
	bool parseOptions(int argc, char** argv);
	static std::string Usage()
	{
	  ostringstream endnote;
	  endnote << "Usage: Fido [options] <graph file>" << endl;
	  return endnote.str();
	}
	int run();
	
 private:
  
	Fido *fido;
	std::string fname;
	std::string targetOut;
	std::string decoyOut;
	/*fido parameters*/
	double fido_alpha;
	double fido_beta;
	double fido_gamma;
	bool fido_grouProteins; 
	bool fido_noprune;
	bool fido_noseparate;
	bool fido_reduceTree;
	bool fido_truncate;
	unsigned fido_depth;
	double fido_mse_threshold;
	/* general protein probabilities options */
	bool tiesAsOneProtein;
	bool usePi0;
	std::string decoy_prefix;
	FidoInterface *protEstimator;
};

int main(int argc, char **argv);

#endif /* FIDO_MAIN_H_ */