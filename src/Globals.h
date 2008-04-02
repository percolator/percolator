/*******************************************************************************
 * Percolator unofficial version
 * Copyright (c) 2006-7 University of Washington. All rights reserved.
 * Written by Lukas Käll (lukall@u.washington.edu) in the 
 * Department of Genome Science at the University of Washington. 
 *
 * $Id: Globals.h,v 1.8 2008/04/02 06:36:41 cegrant Exp $
 *******************************************************************************/
#ifndef GLOBALS_H_
#define GLOBALS_H_

#ifdef WIN32
  #define C_DARRAY(name,nelem) double *name = (double *) _malloca((nelem) * sizeof(double));
  #define D_DARRAY(name) _freea(name);
#else
  #define C_DARRAY(name,nelem) double name[nelem];
  #define D_DARRAY(name)
#endif

#define VERB (Globals::getInstance()->getVerbose())

class Globals
{
public:
	virtual ~Globals();
    static Globals * getInstance();
    static void clean();
    int getVerbose() {return verbose;}
    void setVerbose(int verb) {verbose=verb;}
    void decVerbose() {verbose--;}
    void incVerbose() {verbose++;}
private:
	Globals();
    int verbose;
    static Globals * glob; 
};

#endif /*GLOBALS_H_*/
