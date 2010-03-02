/*******************************************************************************
 Copyright 2006-2009 Lukas Käll <lukas.kall@cbr.su.se>

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

 *******************************************************************************/

#ifndef GLOBALS_H_
#define GLOBALS_H_

#ifdef WIN32
#define C_DARRAY(name,nelem) double *name = (double *) _malloca((nelem) * sizeof(double));
#define D_DARRAY(name) _freea(name);
#include <float.h>
#define isfinite _finite
#else
#define C_DARRAY(name,nelem) double name[nelem];
#define D_DARRAY(name)
#endif

#define VERB (Globals::getInstance()->getVerbose())

#define PERCOLATOR_VERSION 1.00

class Globals {
  public:
    virtual ~Globals();
    static Globals * getInstance();
    static void clean();
    int getVerbose() {
      return verbose;
    }
    void setVerbose(int verb) {
      verbose = verb;
    }
    void decVerbose() {
      verbose--;
    }
    void incVerbose() {
      verbose++;
    }
  private:
    Globals();
    int verbose;
    static Globals * glob;
};

#endif /*GLOBALS_H_*/
