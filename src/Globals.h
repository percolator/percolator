/*******************************************************************************
 Copyright (c) 2008 Lukas Käll

 Permission is hereby granted, free of charge, to any person
 obtaining a copy of this software and associated documentation
 files (the "Software"), to deal in the Software without
 restriction, including without limitation the rights to use,
 copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the
 Software is furnished to do so, subject to the following
 conditions:

 The above copyright notice and this permission notice shall be
 included in all copies or substantial portions of the Software. 
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 OTHER DEALINGS IN THE SOFTWARE.
 
 $Id: Globals.h,v 1.9 2008/05/07 21:25:08 lukall Exp $
 
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
