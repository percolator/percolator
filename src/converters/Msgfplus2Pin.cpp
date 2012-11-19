#include "Msgfplus2Pin.h"


Msgfplus2pin::Msgfplus2pin() {
 
}

Msgfplus2pin::~Msgfplus2pin() {
  if(reader)
    delete reader;
  reader = 0;
}



int Msgfplus2pin::run() {
  
  // Content of sqt files is merged: preparing to write it to xml file
  
  ofstream xmlOutputStream;
  xmlOutputStream.open(xmlOutputFN.c_str());
  if(!xmlOutputStream && xmlOutputFN != ""){
    cerr << "Error: invalid path to output file: " << xmlOutputFN << endl;
    cerr << "Please invoke sqt2pin with a valid -o option" << endl;
    return 0;
  }
  
  //initialize reader
  parseOptions.targetFN = targetFN;
  parseOptions.decoyFN = decoyFN;
  parseOptions.call = call;
  parseOptions.spectrumFN = spectrumFile;
  parseOptions.xmlOutputFN = xmlOutputFN;
  reader = new MsfgplusReader(&parseOptions);
  
  reader->init();
  reader->print(xmlOutputStream);
  
  if (VERB>2)
    cerr << "\nAll the input files have been successfully processed"<< endl;

  return true;
}
int main(int argc, char** argv) {
  
  Msgfplus2pin* pmsgfplus2pin = new Msgfplus2pin();
  int retVal = -1;
  
  try
  {
    if (pmsgfplus2pin->parseOpt(argc, argv, Msgfplus2pin::Usage())) {
      retVal = pmsgfplus2pin->run();
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
  }  
    
  delete pmsgfplus2pin;
  Globals::clean();
  return retVal;
}
