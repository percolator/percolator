#include "Sqt2Pin.h"

Sqt2Pin::Sqt2Pin() {
  reader = NULL;
}

Sqt2Pin::~Sqt2Pin() {
  if(reader)
    delete reader;
  reader = 0;
}


int Sqt2Pin::run() {
  
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
  reader = new SqtReader(&parseOptions);
  
  reader->init();
  reader->print((outputFN == "") ? std::cout : outputStream, xmlOutput);
  outputStream.close();
  
  if (VERB>2)
    cerr << "\nAll the input files have been successfully processed"<< endl;

  return true;
}
int main(int argc, char** argv) {
  
  Sqt2Pin* pSqt2Pin = new Sqt2Pin();
  int retVal = EXIT_FAILURE;
  
  try
  {
    if (pSqt2Pin->parseOpt(argc, argv, Sqt2Pin::Usage())) {
      if(pSqt2Pin->run()) retVal = EXIT_SUCCESS;
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
  delete pSqt2Pin;

  Globals::clean();
  return retVal;
}
