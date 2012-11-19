#include "Sequest2Pin.h"


Sequest2Pin::Sequest2Pin() {
 
}

Sequest2Pin::~Sequest2Pin() {
  if(reader)
    delete reader;
  reader = 0;
}


int Sequest2Pin::run() {
  
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
  reader = new SequestReader(&parseOptions);
  
  reader->init();
  reader->print(xmlOutputStream);
  
  if (VERB>2)
    cerr << "\nAll the input files have been successfully processed"<< endl;

  return true;
}
int main(int argc, char** argv) {
  
  Sequest2Pin* pSequest2Pin = new Sequest2Pin();
  int retVal = -1;
  
  try
  {
    if (pSequest2Pin->parseOpt(argc, argv, Sequest2Pin::Usage())) {
      retVal = pSequest2Pin->run();
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
    
  delete pSequest2Pin;
  Globals::clean();
  return retVal;
}
