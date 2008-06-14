#ifndef _SPECTRUM_H
#define _SPECTRUM_H

#include "MSToolkitTypes.h"
#include <vector>
#include <cstring>

using namespace std;

class Spectrum {
 public:
  //Constructors & Destructors
  Spectrum();
  Spectrum(char *);
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
  int				      getCharge();
  MSSpectrumType  getFileType();
  double    			getMZ();
  float		    		getRTime();
  int	      			getScanNumber(bool second=false);
  void            setActivationMethod(MSActivation);
  void			    	setCharge(int);
  void    				setFileType(MSSpectrumType);
  void		    		setMZ(double);
  void				    setRTime(float);
  void    				setScanNumber(int, bool second=false);
  int			      	size();
  int     				sizeZ();
  void		    		sortIntensity();
  void				    sortIntensityRev();
  void    				sortMZ();
  void		    		sortMZRev();

  vector<Peak_T>* getPeaks();

 protected:

 private:
  //Data Members
  vector<Peak_T> *vPeaks;
  vector<ZState> *vZ;
  int		          charge;
  float		        rTime;
  int		          scanNumber;
  int             scanNumber2;
  double	        mz;
  MSSpectrumType  fileType;
  MSActivation    actMethod;

  //Functions
  static int compareIntensity(const void *p1,const void *p2);
  static int compareMZ(const void *p1,const void *p2);
  static int compareIntensityRev(const void *p1,const void *p2);
  static int compareMZRev(const void *p1,const void *p2);

};

#endif

