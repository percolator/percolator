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

#ifndef VERSION_MAJOR
#define VERSION_MAJOR "@CPACK_PACKAGE_VERSION_MAJOR@"
#endif
#ifndef VERSION_MINOR
#define VERSION_MINOR "@CPACK_PACKAGE_VERSION_MINOR@"
#endif
#ifndef VERSION
#define VERSION "@CPACK_PACKAGE_VERSION_MAJOR@.@CPACK_PACKAGE_VERSION_MINOR@"
#endif

#ifndef PIN_VERSION_MAJOR
#define PIN_VERSION_MAJOR @PIN_VERSION_MAJOR@
#endif

#ifndef POUT_VERSION_MAJOR
#define POUT_VERSION_MAJOR @POUT_VERSION_MAJOR@
#endif

#ifndef PIN_VERSION_MINOR
#define PIN_VERSION_MINOR @PIN_VERSION_MINOR@
#endif

#ifndef POUT_VERSION_MINOR
#define POUT_VERSION_MINOR @POUT_VERSION_MINOR@
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

#ifndef TEMP_DIR
#define TEMP_DIR "@TEMP_DIR@"
#endif

#ifdef WIN32
#ifndef isfinite
#define isfinite _finite
#endif
#define C_DARRAY(name,nelem) double *name = (double *) _malloca((nelem) * sizeof(double));
#define D_DARRAY(name) _freea(name);
#include <float.h>
#else
#define C_DARRAY(name,nelem) double name[nelem];
#define D_DARRAY(name)
#endif

#define VERB (Globals::getInstance()->getVerbose())
#define STD_CERR (*(Globals::getInstance()->getLogger())) 

#include <time.h>
#include <string>
#include "Logger.h"

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
    void decVerbose() {
      verbose--;
    }
    void incVerbose() {
      verbose++;
    }
    clock_t checkTimeClock;
    bool timeCheckPoint;
    void checkTime(const std::string& message);
    void setLogFile(const std::string& filename);
    void initLogger();
    int redirectBuffer();
    void unredirectBuffer();
    
  private:
    Globals();
    int verbose;
    static Globals* glob;
    Logger *log;
    std::string fileLog;
    bool buffer_redirected;
    streambuf* save_sbuf_cerr;
    ofstream ferr;
};

#endif /*GLOBALS_H_*/
