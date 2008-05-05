#include "Option.h"
#include "PosteriorEstimator.h"


int main(int argc, char **argv){
  PosteriorEstimator *pCaller = new PosteriorEstimator();
  int retVal = -1;
  if(pCaller->parseOptions(argc,argv))
  {
    pCaller->run();
  }
  delete pCaller;
  return retVal;
}   
