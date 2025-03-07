#include "MSObject.h"
#include <iostream>

using namespace std;

MSObject::MSObject(){
  vSpectrum = new vector<Spectrum>;
  fileName="";
  fileType=Unspecified;
	for(int i=0;i<16;i++)	strncpy(header.header[i],"\0", sizeof(header.header[i]));
};

MSObject::~MSObject(){
  delete vSpectrum;
};

MSObject::MSObject(const MSObject& m){
  std::size_t i;
  vSpectrum = new vector<Spectrum>;

  for(i=0;i<m.vSpectrum->size();i++){
    vSpectrum->push_back(m.vSpectrum->at(i));
  };

  fileType = m.fileType;
  fileName = m.fileName;
	for(i=0;i<16;i++)	strncpy(header.header[i],m.header.header[i], sizeof(header.header[i]));

};

MSObject& MSObject::operator=(const MSObject& m){
  std::size_t i;
  if (this!=&m){
    delete vSpectrum;
    vSpectrum = new vector<Spectrum>;
    for(i=0;i<m.vSpectrum->size();i++){
      vSpectrum->push_back(m.vSpectrum->at(i));
    };
    fileType = m.fileType;
    fileName = m.fileName;
		for(i=0;i<16;i++)	strncpy(header.header[i],m.header.header[i], sizeof(header.header[i]));
  };
  return *this;
};

  
void MSObject::add(Spectrum& s){
  vSpectrum->push_back(s);
};

bool MSObject::addToHeader(char* c){
	if(strlen(c)>127) return false;
	for(int i=0;i<16;i++){
		if(header.header[i][0]=='\0'){
			strncpy(header.header[i], c, sizeof(header.header[i]));
			return true;
		};
	};
	return false;
};

bool MSObject::addToHeader(string s){
	if(s.size()>127) return false;
	for(int i=0;i<16;i++){
		if(header.header[i][0]=='\0'){
			strncpy(header.header[i],&s[0],sizeof(header.header[i]));
			return true;
		};
	};
	return false;
};

Spectrum& MSObject::at(unsigned int i){
  return vSpectrum->at(i);
};

Peak_T& MSObject::at(unsigned int i, unsigned int j){
  return vSpectrum->at(i).at(j);
};

void MSObject::clear(){
  delete vSpectrum;
  vSpectrum = new vector<Spectrum>;
	for(int i=0;i<16;i++) strncpy(header.header[i],"\0", 2);
};

MSHeader& MSObject::getHeader(){
	return header;
};

void MSObject::setHeader(const MSHeader& h){
	for(int i=0;i<16;i++)	strncpy(header.header[i],h.header[i], sizeof(header.header[i]));
};

int MSObject::size(){
  return static_cast<int>(vSpectrum->size());
};

