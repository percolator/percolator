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


int Fido::run() {
  
  time_t startTime;
  clock_t startClock;
  time(&startTime);
  startClock = clock();  


  protEstimator = new FidoInterface(fido_alpha,fido_beta,fido_gamma,fido_grouProteins,fido_noseparate,
				      fido_noprune,fido_depth,fido_reduceTree,fido_truncate,fido_mse_threshold,
				      tiesAsOneProtein,usePi0,outputEmpirQVal,decoy_prefix);
  
  if (VERB > 0)
  {
    cerr << "\nCalculating protein level probabilities with Fido\n";
    cerr << protEstimator->printCopyright();
  }
  
  ifstream fin(fname);
  
  protEstimator->run();
  protEstimator->computeProbabilitiesFromFile(fin);
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


int main(int argc, char** argv) {
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