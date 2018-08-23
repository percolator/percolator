/*******************************************************************************
Copyright 2006-2012 Lukas Käll <lukas.kall@scilifelab.se>

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
#ifdef _OPENMP
#include <omp.h>
#endif

// #define NOMAP

#ifndef GLOBALS_H_
#define GLOBALS_H_

#ifndef PIN_VERSION_MAJOR
  #define PIN_VERSION_MAJOR "@PIN_VERSION_MAJOR@"
#endif

#ifndef POUT_VERSION_MAJOR
  #define POUT_VERSION_MAJOR "@POUT_VERSION_MAJOR@"
#endif

#ifndef PIN_VERSION_MINOR
  #define PIN_VERSION_MINOR "@PIN_VERSION_MINOR@"
#endif

#ifndef POUT_VERSION_MINOR
  #define POUT_VERSION_MINOR "@POUT_VERSION_MINOR@"
#endif

#ifndef PERCOLATOR_IN_NAMESPACE
  #define PERCOLATOR_IN_NAMESPACE "@percolator-in-namespace@"
#endif

#ifndef PERCOLATOR_OUT_NAMESPACE
  #define PERCOLATOR_OUT_NAMESPACE "@percolator-out-namespace@"
#endif

#ifndef WRITABLE_DIR
  #define WRITABLE_DIR "@WRITABLE_DIR@"
#endif

#ifndef PIN_SCHEMA_LOCATION
  #define PIN_SCHEMA_LOCATION "@PIN_SCHEMA_LOCATION@"
#endif

#ifndef POUT_SCHEMA_LOCATION
  #define POUT_SCHEMA_LOCATION "@POUT_SCHEMA_LOCATION@"
#endif

#ifndef ELUDE_MODELS_PATH
  #define ELUDE_MODELS_PATH "@ELUDE_MODELS_PATH@"
#endif

#ifndef M_PI // for some reason Windows removed the Pi constant from math.h
  #define M_PI 3.14159265358979323846
#endif

#ifdef WIN32
  #define C_DARRAY(name,nelem) double *name = (double *) _malloca((nelem) * sizeof(double));
  #define D_DARRAY(name) _freea(name);
  #include <float.h>
#else
  #define C_DARRAY(name,nelem) double name[nelem];
  #define D_DARRAY(name)
#endif

#ifdef _WIN32
  #ifndef __stdcall
    //#define __stdcall
  #endif
#endif
#ifdef _WIN32
  #include <fcntl.h>
#endif

#ifdef _MSC_VER
  #ifndef isfinite
    #define isfinite _finite
  #endif
#endif



#ifdef __APPLE__
  #ifdef __INTEL_COMPILER
    #define isnan(x) _isnan(x)
    #define isinf(x) _isinf(x)
//  #else
//    #define isinf(x) std::isinf(x)
//    #define isnan(x) x!=x
  #endif
#endif


#define VERB (Globals::getInstance()->getVerbose())
#define NO_TERMINATE (Globals::getInstance()->getNoTerminate())
#define STD_CERR (*(Globals::getInstance()->getLogger())) 

#include <time.h>
#include <string>
#include "Logger.h"
#include "MyException.h"




class Globals {

  
  public:
    
    virtual ~Globals();
    static Globals* getInstance();
    static void clean();
    Logger* getLogger();
    
    const std::string& getLogFile(){
      return fileLog;
    }
    int getVerbose() {
      return verbose;
    }
    void setVerbose(int verb) {
      verbose = verb;
    }
    int getNumThreads() {
      return nr_thread;
    }
    void setNumThreads(int threads) {
#ifdef _OPENMP
      // divide supplied threads amongst number of cross-validation folds
      if(nr_thread > omp_get_max_threads()){
        cout << "Num threads " << nr_thread << " greater than " <<
        omp_get_max_threads() << " max system threads, defaulting to system max.\n";
        nr_thread = omp_get_max_threads() / 3;
      } else {
        nr_thread = threads / 3;
      }

      if(nr_thread > 1){ // turn on nested parallelism
        omp_set_nested(1); 
      	omp_set_dynamic(0);
      }
#else
      nr_thread = threads;
#endif
    }

    void decVerbose() {
      verbose--;
    }
    void incVerbose() {
      verbose++;
    }
    void setNoTerminate(bool noTerminate) {
      noTerminate_ = noTerminate;
    }
    bool getNoTerminate() {
      return noTerminate_;
    }
    void setLogFile(const std::string& filename);
    void initLogger();
    int redirectBuffer();
    void unredirectBuffer();
    const std::string getXMLDir(bool isConverter = false);    
  private:
    Globals();
    int verbose;
    int nr_thread;
    bool noTerminate_;
    static Globals* glob;
    Logger *log;
    std::string fileLog;
    bool buffer_redirected;
    streambuf* save_sbuf_cerr;
    ofstream ferr;
};

#endif /*GLOBALS_H_*/

