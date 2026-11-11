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
#include <memory>
#include <string>
#include <assert.h>

#include "MyException.h"

/* A protease: its name as --protein-enzyme takes it, and its cleavage rule.
 * The enzymes and their rules are the rows of one table in Enzyme.cpp. */
class Enzyme {
  public:
    enum EnzymeType {
      NO_ENZYME, TRYPSIN, TRYPSINP, CHYMOTRYPSIN, THERMOLYSIN, PROTEINASEK, PEPSIN, ELASTASE, 
      LYSN, LYSC, ARGC, ASPN, GLUC
    };

    // a value that is not an EnzymeType gives trypsin
    static std::unique_ptr<Enzyme> createEnzyme(EnzymeType enz);
    // by name, in any case; throws MyException for an unknown name
    static std::unique_ptr<Enzyme> createEnzyme(std::string enzyme);
    // the name of an enzyme, e.g. "trypsinp" or "lys-c"
    static std::string getString(EnzymeType enz);

    EnzymeType getEnzymeType() const;
    std::string getStringEnzyme() const;

    size_t countEnzymatic(std::string& peptide) const;
    // whether the enzyme cleaves between the residues n and c; '-' on either
    // side is the end of the protein, where every enzyme "cleaves"
    bool isEnzymatic(const char& n, const char& c) const;
    // for a peptide with flanks, e.g. "K.AEPTIDER.A": both ends are cleaved
    bool isEnzymatic(std::string peptide) const {
      return (isEnzymatic(peptide[0], peptide[2])
          && isEnzymatic(peptide[peptide.length() - 3],
                         peptide[peptide.length() - 1]));
    }

    struct Rule;  // a row of the table
  private:
    explicit Enzyme(const Rule* rule) : rule_(rule) {}
    const Rule* rule_;
};

#endif /* ENZYME_H_ */
