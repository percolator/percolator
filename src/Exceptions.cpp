#include "Exceptions.h"

const char* MS2Exception::what() const throw(){
  return "test";
  //return "The ms2 in input does not appear to contain retention time "
  //    + "information. Please run without -2 option.";
}
