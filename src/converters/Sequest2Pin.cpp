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
  
  ofstream outputStream;
  outputStream.open(outputFN.c_str());
  if(!outputStream && outputFN != ""){
    cerr << "Error: invalid path to output file: " << outputFN << endl;
    cerr << "Please invoke sqt2pin with a valid -o option" << endl;
    return 0;
  }
  
  //initialize reader
  parseOptions.targetFN = targetFN;
  parseOptions.decoyFN = decoyFN;
  parseOptions.call = call;
  parseOptions.spectrumFN = spectrumFile;
  parseOptions.xmlOutputFN = outputFN;
  reader = new SequestReader(&parseOptions);
  
  reader->init();
  reader->print((outputFN == "") ? std::cout : outputStream, xmlOutput);
  outputStream.close();
  
  if (VERB>2)
    cerr << "\nAll the input files have been successfully processed"<< endl;

  return true;
}
int main(int argc, char** argv) {
  
  Sequest2Pin* pSequest2Pin = new Sequest2Pin();
  int retVal = EXIT_FAILURE;
  
  try
  {
    if (pSequest2Pin->parseOpt(argc, argv, Sequest2Pin::Usage())) {
      if(pSequest2Pin->run()) retVal = EXIT_SUCCESS;
    }
  }
  catch (const std::exception& e) 
  {
    std::cerr << e.what() << endl;
    retVal = EXIT_FAILURE;
  }
  catch(...)
  {
    std::cerr << "Unknown exception, contact the developer.." << std::endl;
    retVal = EXIT_FAILURE;
  }  
    
  delete pSequest2Pin;
  Globals::clean();
  return retVal;
}
