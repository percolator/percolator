#include "Spectrum.h"
#include <cstdlib>
#include <iostream>

Spectrum::Spectrum(){
  rTime=0;
  charge=0;
  scanNumber=0;
  scanNumber2=0;
	mz=0;
	fileType=Unspecified;
  vPeaks = new vector<Peak_T>;
	vZ = new vector<ZState>;
  actMethod=na;
};

Spectrum::~Spectrum(){
  delete vPeaks;
	delete vZ;
};

Spectrum::Spectrum(const Spectrum& s){
  size_t i;

  rTime = s.rTime;
  charge = s.charge;
  scanNumber = s.scanNumber;
  scanNumber2 = s.scanNumber2;
	mz = s.mz;
	fileType = s.fileType;
  vPeaks = new vector<Peak_T>;
  for(i=0;i<s.vPeaks->size();i++){
    vPeaks->push_back(s.vPeaks->at(i));
  };
	vZ = new vector<ZState>;
  for(i=0;i<s.vZ->size();i++){
    vZ->push_back(s.vZ->at(i));
  };
};

Spectrum& Spectrum::operator=(const Spectrum& s){
  size_t i;
  if (this != &s) {
    delete vPeaks;
		delete vZ;
    vPeaks = new vector<Peak_T>;
    for(i=0;i<s.vPeaks->size();i++){
      vPeaks->push_back(s.vPeaks->at(i));
    };
		vZ = new vector<ZState>;
		for(i=0;i<s.vZ->size();i++){
			vZ->push_back(s.vZ->at(i));
		};
    rTime = s.rTime;
    charge = s.charge;
    scanNumber = s.scanNumber;
    scanNumber2 = s.scanNumber2;
		mz = s.mz;
		fileType = s.fileType;
  };
  return *this;
};

Peak_T& Spectrum::operator[](const int& i) {
	return vPeaks->operator[](i);
};


/* ----- Begin Functions ----- */


/* Adds Result struct to end of spectrum. */
void Spectrum::add(Peak_T& p){
  vPeaks->push_back(p);
};

void Spectrum::add(double d1, float d2){
  Peak_T p;
  p.mz=d1;
  p.intensity=d2;
  vPeaks->push_back(p);
};

void Spectrum::addZState(ZState& z){
	vZ->push_back(z);
};

void Spectrum::addZState(int i, double d){
	ZState z;
	z.z=i;
	z.mz=d;
	vZ->push_back(z);
};

/* Returns Result struct of single element in the spectrum. */
Peak_T& Spectrum::at(const int& i){
  return vPeaks->operator [](i);
};

Peak_T& Spectrum::at(const unsigned int& i){
  return vPeaks->operator [](i);
};

ZState& Spectrum::atZ(const int& i){
	return vZ->operator [](i);
};

ZState& Spectrum::atZ(const unsigned int& i){
	return vZ->operator [](i);
};

/* Clears the spectrum */
void Spectrum::clear(){
	delete vPeaks;
	vPeaks = new vector<Peak_T>;
	delete vZ;
	vZ = new vector<ZState>;
	scanNumber = 0;
  scanNumber2 = 0;
	rTime = 0;
	charge = 0;
	mz = 0;
	fileType = Unspecified;
  actMethod=na;
};

/* Erases element i in the spectrum. */
void Spectrum::erase(unsigned int i){
  vector<Peak_T>::iterator vi;
  vi=vPeaks->begin()+i;
  vPeaks->erase(vi);
};

/* Erases element i to element j, inclusive, in the spectrum. */
void Spectrum::erase(unsigned int i, unsigned int j){
  vector<Peak_T>::iterator vi1;
  vector<Peak_T>::iterator vi2;
  vi1=vPeaks->begin()+i;
  vi2=vPeaks->begin()+j+1;
  vPeaks->erase(vi1,vi2);
};

void Spectrum::eraseZ(unsigned int i){
  vector<ZState>::iterator vi;
  vi=vZ->begin()+i;
  vZ->erase(vi);
};

/* Erases element i to element j, inclusive, in the spectrum. */
void Spectrum::eraseZ(unsigned int i, unsigned int j){
  vector<ZState>::iterator vi1;
  vector<ZState>::iterator vi2;
  vi1=vZ->begin()+i;
  vi2=vZ->begin()+j+1;
  vZ->erase(vi1,vi2);
};

int Spectrum::getCharge(){
  return charge;
};

MSSpectrumType Spectrum::getFileType(){
	return fileType;
};

double Spectrum::getMZ(){
  return mz;
};

float Spectrum::getRTime(){
  return rTime;
};

int Spectrum::getScanNumber(bool second){
  if(second) return scanNumber2;
  else return scanNumber;
};

void Spectrum::setCharge(int i){
  charge=i;
};

void Spectrum::setFileType(MSSpectrumType f){
	fileType=f;
};

void Spectrum::setMZ(double d){
  mz=d;
};

void Spectrum::setRTime(float d){
  rTime=d;
};

void Spectrum::setScanNumber(int i, bool second){
  if(second)scanNumber2=i;
  else scanNumber=i;
};

/* Returns the number of elements in the spectrum. */
int Spectrum::size(){
  return vPeaks->size();
};

int Spectrum::sizeZ(){
	return vZ->size();
};

/* Sorts the spectrum by Data. */
void Spectrum::sortIntensity(){
  qsort(&vPeaks->at(0),vPeaks->size(),sizeof(Peak_T),compareIntensity);
};

/* Sorts the spectrum by Mass. */
void Spectrum::sortMZ(){
  qsort(&vPeaks->at(0),vPeaks->size(),sizeof(Peak_T),compareMZ);
};

/* Sorts the spectrum by Data. */
void Spectrum::sortIntensityRev(){
  qsort(&vPeaks->at(0),vPeaks->size(),sizeof(Peak_T),compareIntensityRev);
};

/* Sorts the spectrum by Mass. */
void Spectrum::sortMZRev(){
  qsort(&vPeaks->at(0),vPeaks->size(),sizeof(Peak_T),compareMZRev);
};

vector<Peak_T>* Spectrum::getPeaks(){
	return vPeaks;
};

MSActivation Spectrum::getActivationMethod(){
  return actMethod;
};

void Spectrum::setActivationMethod(MSActivation m){
  actMethod=m;
};


//Private Functions

/* For the qsort */
int Spectrum::compareIntensity(const void *p1, const void *p2){
  const Peak_T d1 = *(Peak_T *)p1;
  const Peak_T d2 = *(Peak_T *)p2;
  if(d1.intensity<d2.intensity) return -1;
  else if(d1.intensity>d2.intensity) return 1;
  else return 0;
};

/* For the qsort */
int Spectrum::compareMZ(const void *p1, const void *p2){
  const Peak_T d1 = *(Peak_T *)p1;
  const Peak_T d2 = *(Peak_T *)p2;
  if(d1.mz<d2.mz) return -1;
  else if(d1.mz>d2.mz) return 1;
  else return 0;
};

/* For the qsort */
int Spectrum::compareIntensityRev(const void *p1, const void *p2){
  const Peak_T d1 = *(Peak_T *)p1;
  const Peak_T d2 = *(Peak_T *)p2;
  if(d1.intensity>d2.intensity) return -1;
  else if(d1.intensity<d2.intensity) return 1;
  else return 0;
};

/* For the qsort */
int Spectrum::compareMZRev(const void *p1, const void *p2){
  const Peak_T d1 = *(Peak_T *)p1;
  const Peak_T d2 = *(Peak_T *)p2;
  if(d1.mz>d2.mz) return -1;
  else if(d1.mz<d2.mz) return 1;
  else return 0;
};
