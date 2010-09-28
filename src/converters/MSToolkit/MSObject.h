#ifndef _MSOBJECT_H
#define _MSOBJECT_H

#include "MSToolkitTypes.h"
#include "Spectrum.h"

#include <string>

using namespace std;

class MSObject {
 public:
  //Constructors & Destructors
  MSObject();
  MSObject(const MSObject&);
  ~MSObject(); 

  //Operator Functions
  MSObject& operator=(const MSObject&);
  
  //Functions
  void add(Spectrum&);
  bool addToHeader(char*);
  bool addToHeader(string);
  Spectrum& at(unsigned int);
  Peak_T& at(unsigned int, unsigned int);
  void clear();
  void erase(unsigned int);
  void erase(unsigned int, unsigned int);
  MSHeader& getHeader();
  void setHeader(const MSHeader& h);
  int size();
  
 protected:
 private:
  vector<Spectrum> *vSpectrum;
  string fileName;
	MSHeader header;
  MSSpectrumType fileType;
  
};

#endif

