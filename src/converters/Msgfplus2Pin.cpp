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

  ofstream outputStream;
  outputStream.open(outputFN.c_str());
  if(!outputStream && outputFN != ""){
    cerr << "Error: invalid path to output file: " << outputFN << endl;
    cerr << "Please invoke msgf2pin with a valid -o option" << endl;
    return 0;
  }

  //initialize reader
  parseOptions.targetFN = targetFN;
  parseOptions.decoyFN = decoyFN;
  parseOptions.call = call;
  parseOptions.spectrumFN = spectrumFile;
  parseOptions.xmlOutputFN = outputFN;
  reader = new MsgfplusReader(&parseOptions);

  reader->init();
  reader->print((outputFN == "") ? std::cout : outputStream, xmlOutput);
  outputStream.close();

  return true;
}
int main(int argc, char** argv) {

  Msgfplus2pin* pmsgfplus2pin = new Msgfplus2pin();
  int retVal = EXIT_FAILURE;

  try
  {
    if (pmsgfplus2pin->parseOpt(argc, argv, Msgfplus2pin::Usage())) {
      if(pmsgfplus2pin->run()) retVal = EXIT_SUCCESS;
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

  delete pmsgfplus2pin;
  Globals::clean();
  return retVal;
}
