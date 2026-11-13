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

#include "Enzyme.h"

#include <algorithm>
#include <cctype>
#include <sstream>

#include "MyException.h"

struct Enzyme::Rule {
  EnzymeType type;
  const char* name;
  // whether the enzyme cleaves between the residues n and c inside a protein
  bool (*cleaves)(char n, char c);
};

namespace {

inline bool isOneOf(char residue, const char* residues) {
  for (; *residues != '\0'; ++residues) {
    if (residue == *residues) return true;
  }
  return false;
}

// One row per enzyme, in the order of EnzymeType.
const Enzyme::Rule kRules[] = {
  {Enzyme::NO_ENZYME, "no_enzyme",
   [](char, char) { return true; }},
  {Enzyme::TRYPSIN, "trypsin",
   [](char n, char c) { return isOneOf(n, "KR") && c != 'P'; }},
  {Enzyme::TRYPSINP, "trypsinp",
   [](char n, char) { return isOneOf(n, "KR"); }},
  {Enzyme::CHYMOTRYPSIN, "chymotrypsin",
   [](char n, char c) { return isOneOf(n, "FWYL") && c != 'P'; }},
  {Enzyme::THERMOLYSIN, "thermolysin",
   [](char n, char c) {
     return (isOneOf(c, "AFILMV") || (n == 'R' && c == 'G')) &&
            n != 'D' && n != 'E';
   }},
  {Enzyme::PROTEINASEK, "proteinasek",
   [](char n, char) { return isOneOf(n, "AEFILTVWY"); }},
  {Enzyme::PEPSIN, "pepsin",
   [](char n, char c) {
     return (isOneOf(c, "FLWY") || isOneOf(n, "FLWY")) && n != 'R';
   }},
  {Enzyme::ELASTASE, "elastase",
   [](char n, char c) { return isOneOf(n, "ALIV") && c != 'P'; }},
  {Enzyme::LYSN, "lys-n",
   [](char, char c) { return c == 'K'; }},
  {Enzyme::LYSC, "lys-c",
   [](char n, char c) { return n == 'K' && c != 'P'; }},
  {Enzyme::ARGC, "arg-c",
   [](char n, char c) { return n == 'R' && c != 'P'; }},
  {Enzyme::ASPN, "asp-n",
   [](char, char c) { return c == 'D'; }},
  {Enzyme::GLUC, "glu-c",
   [](char n, char c) { return isOneOf(n, "DE") && c != 'P'; }},
};

const Enzyme::Rule* findRule(Enzyme::EnzymeType enz) {
  for (const Enzyme::Rule& rule : kRules) {
    if (rule.type == enz) return &rule;
  }
  return &kRules[Enzyme::TRYPSIN];
}

}  // namespace

std::unique_ptr<Enzyme> Enzyme::createEnzyme(EnzymeType enz) {
  return std::unique_ptr<Enzyme>(new Enzyme(findRule(enz)));
}

std::unique_ptr<Enzyme> Enzyme::createEnzyme(std::string enzyme) {
  std::transform(enzyme.begin(), enzyme.end(), enzyme.begin(), ::tolower);
  for (const Rule& rule : kRules) {
    if (enzyme == rule.name) {
      return std::unique_ptr<Enzyme>(new Enzyme(&rule));
    }
  }
  std::ostringstream temp;
  temp << "The selected enzyme \"" << enzyme << "\" has no corresponding class" << std::endl;
  throw MyException(temp.str());
}

std::string Enzyme::getString(EnzymeType enz) {
  return findRule(enz)->name;
}

Enzyme::EnzymeType Enzyme::getEnzymeType() const {
  return rule_->type;
}

std::string Enzyme::getStringEnzyme() const {
  return rule_->name;
}

bool Enzyme::isEnzymatic(const char& n, const char& c) const {
  return rule_->cleaves(n, c) || n == '-' || c == '-';
}

size_t Enzyme::countEnzymatic(std::string& peptide) const {
  size_t count = 0;
  for (size_t ix = 1; ix < peptide.size(); ++ix) {
    if (isEnzymatic(peptide[ix - 1], peptide[ix])) {
      ++count;
    }
  }
  return count;
}

