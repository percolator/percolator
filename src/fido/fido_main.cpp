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
#include "Option.h"
#include "Globals.h"
#include "fido_main.h"


Fido::Fido() {
 
}

Fido::~Fido() {
  if(fido)
    delete fido;
  fido = 0;
  
  if(protEstimator)
    delete protEstimator;
  protEstimator = 0;
}

string Fido::greeter()
{
  ostringstream oss;
  oss << "Copyright (c) 2008-9 University of Washington. All rights reserved.\n"
      << "Written by Oliver R. Serang (orserang@u.washington.edu) in the\n"
      << "Department of Genome Sciences at the University of Washington.\n" << std::endl;
  return oss.str();
}

bool Fido::parseOptions(int argc, char** argv) {
  // init
  ostringstream intro;
  intro << greeter() << endl;
  intro << "Usage:" << endl;
  intro << "   Fido [options] <graph file>" << endl;
  CommandLineParser cmd(intro.str());
  // finally parse and handle return codes (display help etc...)
  
  cmd.defineOption("v",
      "verbose",
      "Set verbosity of output: 0=no processing info, 5=all, default is 2",
      "level");
  cmd.defineOption("r",
      "results",
      "Output tab delimited results to a file instead of stdout",
      "filename");
  cmd.defineOption("B",
      "decoy-results",
      "Output tab delimited results for decoys into a file",
      "filename");
  cmd.defineOption("a",
      "alpha",
      "Probability with which a present protein emits an associated peptide (to be used jointly with the -A option) \
       Set by grid search if not specified.",
      "value");
  cmd.defineOption("b",
      "beta",
      "Probability of the creation of a peptide from noise (to be used jointly with the -A option). Set by grid search if not specified",
      "value");
  cmd.defineOption("G",
      "gamma",
      "Prior probability of that a protein is present in the sample ( to be used with the -A option). Set by grid search if not specified",
      "value");
  cmd.defineOption("g",
      "allow-protein-group",
      "treat ties as if it were one protein (Only valid if option -A is active).",
      "",
      TRUE_IF_SET);
  cmd.defineOption("I",
      "protein-level-pi0",
      "use pi_0 value when calculating empirical q-values (no effect if option Q is activated) (Only valid if option -A is active).",
      "", 
      TRUE_IF_SET);
  cmd.defineOption("N",
      "group-proteins", 		   
      "activates the grouping of proteins with similar connectivity, \
       for example if proteins P1 and P2 have the same peptides matching both of them, P1 and P2 can be grouped as one protein \
       (Only valid if option -A is active).",
      "",
      TRUE_IF_SET);
  cmd.defineOption("E",
      "no-separate-proteins", 		   
      "Proteins graph will not be separated in sub-graphs (Only valid if option -A is active).",
      "",
      TRUE_IF_SET); 
  cmd.defineOption("C",
      "no-prune-proteins", 		   
      "it does not prune peptides with a very low score (~0.0) which means that if a peptide with a very low score is matching two proteins,\
       when we prune the peptide,it will be duplicated to generate two new protein groups (Only valid if option -A is active).",
      "",
      TRUE_IF_SET);
  cmd.defineOption("d",
      "gridsearch-depth",
      "Setting depth 0 or 1 or 2 or 3 from high depth to low depth(less computational time) \
       of the grid search for the estimation Alpha,Beta and Gamma parameters for fido(Only valid if option -A is active). Default value is 3",
      "value");
  cmd.defineOption("P",
      "pattern",
      "Define the text pattern to identify the decoy proteins and/or PSMs, set this up if the label that idenfifies the decoys in the database \
       is not the default (by default : ramdom) (Only valid if option -A  is active).",
      "value");
  cmd.defineOption("T",
      "fido-reduce-tree-in-gridsearch",
      "Reduce the tree of proteins (removing low scored proteins) in order to estimate alpha,beta and gamma faster.(Only valid if option -A is active).",
      "",
      TRUE_IF_SET);
  cmd.defineOption("H",
      "grid-search-mse-threshold",
      "Q-value threshold that will be used in the computation of the MSE and ROC AUC score in the grid search (recommended 0.05 for normal size datasets and 0.1 for big size datasets).(Only valid if option -A is active).",
      "",
      "value");
  cmd.defineOption("W",
      "no-truncation",
      "Proteins with a very low score (< 0.001) will not be truncated (assigned 0.0 probability).(Only valid if option -A is active).",
      "",
      FALSE_IF_SET);
  
  
  // finally parse and handle return codes (display help etc...)
  cmd.parseArgs(argc, argv);
  // now query the parsing results
  if (cmd.optionSet("B")) {
    decoyOut = cmd.options["B"];
  }
  
  if (cmd.optionSet("r")) {
    targetOut = cmd.options["r"];
  }
  
  if (cmd.optionSet("v")) {
    Globals::getInstance()->setVerbose(cmd.getInt("v", 0, 10));
  }
  
  tiesAsOneProtein = cmd.optionSet("g");
  usePi0 = cmd.optionSet("I");
  fido_grouProteins = cmd.optionSet("N"); 
  fido_noprune = cmd.optionSet("C");
  fido_noseparate = cmd.optionSet("E");
  fido_reduceTree = cmd.optionSet("T");
  fido_truncate = cmd.optionSet("W");
  if (cmd.optionSet("P"))  decoy_prefix = cmd.options["P"];
  if (cmd.optionSet("d"))  fido_depth = cmd.getInt("d", 0, 3);
  if (cmd.optionSet("a"))  fido_alpha = cmd.getDouble("a", 0.00, 1.0);
  if (cmd.optionSet("b"))  fido_beta = cmd.getDouble("b", 0.00, 1.0);
  if (cmd.optionSet("G"))  fido_gamma = cmd.getDouble("G", 0.00, 1.0);
  if (cmd.optionSet("H"))  fido_mse_threshold = cmd.getDouble("H",0.001,1.0);

  fname = cmd.arguments[0];

  return true;
}


int Fido::run() {
  
  time_t startTime;
  clock_t startClock;
  time(&startTime);
  startClock = clock();  


  protEstimator = new FidoInterface(fido_alpha,fido_beta,fido_gamma,fido_grouProteins,fido_noseparate,
				      fido_noprune,fido_depth,fido_reduceTree,fido_truncate,fido_mse_threshold,
				      tiesAsOneProtein,usePi0,false/*outputEmpirQVal*/,decoy_prefix);
  
  if (VERB > 0)
  {
    cerr << "\nCalculating protein level probabilities with Fido\n";
    cerr << protEstimator->printCopyright();
  }
  
  protEstimator->run();
  protEstimator->computeProbabilities(std::string(fname));
  protEstimator->computeStatistics();
  
  time_t procStart;
  clock_t procStartClock = clock();
  time(&procStart);
  double diff_time = difftime(procStart, startTime);
  
  if (VERB > 1) 
  {  
    cerr << "Estimating Protein Probabilities took : "
    << ((double)(procStartClock - startClock)) / (double)CLOCKS_PER_SEC
    << " cpu seconds or " << diff_time << " seconds wall time" << endl;
  }
  
  protEstimator->printOut(targetOut,decoyOut);

  return true;
}


int main(int argc, char** argv) 
{
    Fido* pCaller = new Fido();
    int retVal = -1;
    try
    {
     if (pCaller->parseOptions(argc, argv)) {
       retVal=pCaller->run();
     }
    } 
    catch (const std::exception& e) 
    {
      std::cerr << e.what() << endl;
      retVal = -1;
    }
    catch(...)
    {
      std::cerr << "Unknown exception, contact the developer.." << std::endl;
      retVal = -1;
    }
    delete pCaller;
    Globals::clean();
    return retVal;
}
