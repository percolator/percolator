#include <vector>
#include <set>
#include <map>
#include <string>
using namespace std;
#include "DataSet.h"
#include "Scores.h"
#include "SetHandler.h"
#include "Caller.h"
#include "Globals.h"

int main(int argc, char **argv){
  Caller *pCaller = new Caller();
  int retVal = -1;
  if(pCaller->parseOptions(argc,argv))
  {
    retVal = pCaller->run();
  }
  delete pCaller;
  Globals::clean();
  return retVal;
}   
