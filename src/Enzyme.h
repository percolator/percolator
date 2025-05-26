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

#ifndef ENZYME_H_
#define ENZYME_H_

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <cstring> // needed for stricmp()
#include <string>
#include <assert.h>

#include "MyException.h"

class Enzyme {
  public:
    enum EnzymeType {
      NO_ENZYME, TRYPSIN, TRYPSINP, CHYMOTRYPSIN, THERMOLYSIN, PROTEINASEK, PEPSIN, ELASTASE, 
      LYSN, LYSC, ARGC, ASPN, GLUC
    };
    Enzyme() {}
    virtual ~Enzyme() {}
    
    static Enzyme* createEnzyme(EnzymeType enz);
    static Enzyme* createEnzyme(std::string enzyme);
    
    EnzymeType getEnzymeType() const {
      return getET();
    }
    size_t countEnzymatic(std::string& peptide) const;
    bool isEnzymatic(const char& n, const char& c) const {
      return isEnz(n, c);
    }
    bool isEnzymatic(std::string peptide) const {
      return (isEnz(peptide[0], peptide[2])
          && isEnz(peptide[peptide.length() - 3],
                                peptide[peptide.length() - 1]));
    }
    
    std::string getStringEnzyme() const {
      return toString();
    };
    
  protected:
    virtual bool isEnz(const char& n, const char& c) const = 0;
    virtual std::string toString() const = 0;
    virtual EnzymeType getET() const = 0;
    
  private:
    
    inline int stricmp (const std::string &s1,const std::string &s2)
    {
      return stricmp (s1.c_str(), s2.c_str()); // C's stricmp
    }
};

class NoEnzyme : public Enzyme {
  public:
    virtual ~NoEnzyme() {
      ;
    }
    NoEnzyme() {
      ;
    }
    static std::string getString() {
      return "no_enzyme";
    }
  protected:
    virtual std::string toString() const {
      return getString();
    }
    virtual bool isEnz(const char& n, const char& c) const {
      (void)n;
      (void)c;
      return true;
    }
    virtual EnzymeType getET() const {
      return NO_ENZYME;
    }
};

class Trypsin : public Enzyme {
  public:
    virtual ~Trypsin() {
      ;
    }
    Trypsin() {
      ;
    }
    static std::string getString() {
      return "trypsin";
    }
  protected:
    virtual std::string toString() const {
      return getString();
    }
    virtual bool isEnz(const char& n, const char& c) const {
      return (((n == 'K' || n == 'R') && c != 'P') || n == '-' || c == '-');
    }
    virtual EnzymeType getET() const {
      return TRYPSIN;
    }
};

class TrypsinP : public Enzyme {
  public:
    virtual ~TrypsinP() {
      ;
    }
    TrypsinP() {
      ;
    }
    static std::string getString() {
      return "trypsinp";
    }
  protected:
    virtual std::string toString() const {
      return getString();
    }
    virtual bool isEnz(const char& n, const char& c) const {
      return ((n == 'K' || n == 'R') || n == '-' || c == '-');
    }
    virtual EnzymeType getET() const {
      return TRYPSINP;
    }
};

class Chymotrypsin : public Enzyme {
  public:
    virtual ~Chymotrypsin() {
      ;
    }
    Chymotrypsin() {
      ;
    }
    static std::string getString() {
      return "chymotrypsin";
    }
  protected:
    virtual std::string toString() const {
      return getString();
    }
    virtual bool isEnz(const char& n, const char& c) const {
      return (((n == 'F' || n == 'W' || n == 'Y' || n == 'L' )
        && c != 'P') || n == '-' || c == '-');
    }
    virtual EnzymeType getET() const {
      return CHYMOTRYPSIN;
    }
};

class Thermolysin : public Enzyme {
  public:
    virtual ~Thermolysin() {
      ;
    }
    Thermolysin() {
      ;
    }
    static std::string getString() {
      return "thermolysin";
    }
  protected:
    virtual std::string toString() const {
      return getString();
    }
    virtual bool isEnz(const char& n, const char& c) const {
      return (((c == 'A' || c == 'F' || c == 'I' || c == 'L' || c == 'M'
          || c == 'V' || (n == 'R' && c == 'G')) && n != 'D' && n != 'E') || n == '-' || c == '-');
    }
    virtual EnzymeType getET() const {
      return THERMOLYSIN;
    }
};

class Proteinasek : public Enzyme {
  public:
    virtual ~Proteinasek() {
      ;
    }
    Proteinasek() {
      ;
    }
    static std::string getString() {
      return "proteinasek";
    }
  protected:
    virtual std::string toString() const {
      return getString();
    }
    virtual bool isEnz(const char& n, const char& c) const {
      return ((n == 'A' || n == 'E' || n == 'F' || n == 'I' || n == 'L'
          || n == 'T' || n == 'V' || n == 'W' || n == 'Y' ) || n == '-' || c == '-');
    }
    virtual EnzymeType getET() const {
      return PROTEINASEK;
    }
};

class Pepsin : public Enzyme {
  public:
    virtual ~Pepsin() {
      ;
    }
    Pepsin() {
      ;
    }
    static std::string getString() {
      return "pepsin";
    }
  protected:
    virtual std::string toString() const {
      return getString();
    }
    virtual bool isEnz(const char& n, const char& c) const {
      return (((c == 'F' || c == 'L' || c == 'W' || c == 'Y' || n == 'F'
          || n == 'L' || n == 'W' || n == 'Y') && n != 'R') || n == '-' || c == '-');
    }
    virtual EnzymeType getET() const {
      return PEPSIN;
    }
};

class Elastase : public Enzyme {
  public:
    virtual ~Elastase() {
      ;
    }
    Elastase() {
      ;
    }
    static std::string getString() {
      return "elastase";
    }
  protected:
    virtual std::string toString() const {
      return getString();
    }
    virtual bool isEnz(const char& n, const char& c) const {
      return (((n == 'L' || n == 'V' || n == 'A' || n == 'G') && c != 'P')
          || n == '-' || c == '-');
    }
    virtual EnzymeType getET() const {
      return ELASTASE;
    }
};

class LysN : public Enzyme {
  public:
    virtual ~LysN() {
      ;
    }
    LysN() {
      ;
    }
    static std::string getString() {
      return "lys-n";
    }
  protected:
    virtual std::string toString() const {
      return getString();
    }
    virtual bool isEnz(const char& n, const char& c) const {
      return ((c == 'K')
          || n == '-' || c == '-');
    }
    virtual EnzymeType getET() const {
      return LYSN;
    }
};

class LysC : public Enzyme {
  public:
    virtual ~LysC() {
      ;
    }
    LysC() {
      ;
    }
    static std::string getString() {
      return "lys-c";
    }
  protected:
    virtual std::string toString() const {
      return getString();
    }
    virtual bool isEnz(const char& n, const char& c) const {
      return (((n == 'K') && c != 'P')
          || n == '-' || c == '-');
    }
    virtual EnzymeType getET() const {
      return LYSC;
    }
};

class ArgC : public Enzyme {
  public:
    virtual ~ArgC() {
      ;
    }
    ArgC() {
      ;
    }
    static std::string getString() {
      return "arg-c";
    }
  protected:
    virtual std::string toString() const {
      return getString();
    }
    virtual bool isEnz(const char& n, const char& c) const {
      return (((n == 'R') && c != 'P')
          || n == '-' || c == '-');
    }
    virtual EnzymeType getET() const {
      return ARGC;
    }
};

class AspN : public Enzyme {
  public:
    virtual ~AspN() {
      ;
    }
    AspN() {
      ;
    }
    static std::string getString() {
      return "asp-n";
    }
  protected:
    virtual std::string toString() const {
      return getString();
    }
    virtual bool isEnz(const char& n, const char& c) const {
      return ((c == 'D')
          || n == '-' || c == '-');
    }
    virtual EnzymeType getET() const {
      return ASPN;
    }
};

class GluC : public Enzyme {
  public:
    virtual ~GluC() {
      ;
    }
    GluC() {
      ;
    }
    static std::string getString() {
      return "glu-c";
    }
  protected:
    virtual std::string toString() const {
      return getString();
    }
    virtual bool isEnz(const char& n, const char& c) const {
      return (((n == 'E') && (c != 'P'))
          || n == '-' || c == '-');
    }
    virtual EnzymeType getET() const {
      return GLUC;
    }
};

#endif /* ENZYME_H_ */
