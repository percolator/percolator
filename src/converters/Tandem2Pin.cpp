#include "Tandem2Pin.h"

Tandem2Pin::Tandem2Pin() {
  
}

Tandem2Pin::~Tandem2Pin() {
  if(reader)
    delete reader;
  reader = 0;
}

int Tandem2Pin::run() {
  // Content of tandem files is merged: preparing to write it to xml file
  ofstream outputStream;
  outputStream.open(outputFN.c_str());
  if(!outputStream && outputFN != ""){
    cerr << "Error: invalid path to output file: " << outputFN << endl;
    cerr << "Please invoke tandem2pin with a valid -o option" << endl;
    return 0;
  }
  
  //initialize reader
  parseOptions.targetFN = targetFN;
  parseOptions.decoyFN = decoyFN;
  parseOptions.call = call;
  parseOptions.spectrumFN = spectrumFile;
  parseOptions.xmlOutputFN = outputFN;
  reader = new TandemReader(&parseOptions);
  
  reader->init();
  reader->print((outputFN == "") ? std::cout : outputStream, xmlOutput);
  outputStream.close();
  
  if (VERB>2)
    cerr << "\nAll the input files have been successfully processed"<< endl;

  return true;
}

int main(int argc, char** argv) {
  Tandem2Pin* ptandem2Pin = new Tandem2Pin();
  int retVal = -1;
  
  try
  {
    if (ptandem2Pin->parseOpt(argc, argv, Tandem2Pin::Usage())) {
      retVal = ptandem2Pin->run();
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
  
  delete ptandem2Pin;
  Globals::clean();
  return retVal;
}
