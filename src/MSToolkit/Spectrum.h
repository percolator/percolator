#ifndef _SPECTRUM_H
#define _SPECTRUM_H

#include "MSToolkitTypes.h"
#include <vector>
#include <cstring>
#include <iomanip>

using namespace std;

class Spectrum {
 public:
  //Constructors & Destructors
  Spectrum();
  Spectrum(char*);
  Spectrum(char, unsigned int);
  Spectrum(const Spectrum&);
  ~Spectrum();

  //Operator Functions
  Spectrum& operator=(const Spectrum&);
	Peak_T& operator[](const int&);

  //Functions
  void	    			add(Peak_T&);
  void			    	add(double,float);
  void    				addZState(int,double);
  void		    		addZState(ZState&);
  Peak_T&			    at(const int&);
  Peak_T&	    		at(const unsigned int&);
  ZState&			    atZ(const int&);
  ZState&	    		atZ(const unsigned int&);
  void			    	clear();
  void				    erase(unsigned int);
  void				    erase(unsigned int, unsigned int);
  void				    eraseZ(unsigned int);
  void				    eraseZ(unsigned int, unsigned int);
  MSActivation    getActivationMethod();
  float           getBPI();
  double          getBPM();
  int				      getCharge();
  double          getConversionA();
  double          getConversionB();
  MSSpectrumType  getFileType();
  float           getIonInjectionTime();
  double    			getMZ();
  float		    		getRTime();
  int	      			getScanNumber(bool second=false);
  double          getTIC();
  int             getMsLevel();
  void            setActivationMethod(MSActivation);
  void            setBPI(float);
  void            setBPM(double);
  void			    	setCharge(int);
  void            setConversionA(double);
  void            setConversionB(double);
  void    				setFileType(MSSpectrumType);
  void            setIonInjectionTime(float);
  void		    		setMZ(double);
  void				    setRTime(float);
  void    				setScanNumber(int, bool second=false);
  void            setTIC(double);
  void            setMsLevel(int level);
  int			      	size();
  int     				sizeZ();
  void		    		sortIntensity();
  void				    sortIntensityRev();
  void    				sortMZ();
  void            setPeaks( std::vector<Peak_T> peaks);
  void		    		sortMZRev();

  //for sqlite format
  void setScanID(int scanID);
  int getScanID();

  vector<Peak_T>* getPeaks();
  //void setPeaks(vector<Peak_T> peaks);
  float getTotalIntensity();
 
  //for debugging
  void printMe();

 protected:

 //Data Members
  vector<Peak_T> *vPeaks;
  vector<ZState> *vZ;
  int		          charge;
  float		        rTime;
  int		          scanNumber;
  int             scanNumber2;
  int             msLevel;
  double	        mz;
  MSSpectrumType  fileType;
  MSActivation    actMethod;
  int             scanID;       //index for sqlite
  float           IIT;
  float           BPI;          //Base Peak Intensity
  double          convA;
  double          convB;
  double          TIC;
  double          BPM;          //Base Peak Mass

  //private:
  //Functions
  static int compareIntensity(const void *p1,const void *p2);
  static int compareMZ(const void *p1,const void *p2);
  static int compareIntensityRev(const void *p1,const void *p2);
  static int compareMZRev(const void *p1,const void *p2);

};


#endif

