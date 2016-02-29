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
#include <cmath>
#include <assert.h>

#include "Globals.h"
#include "PSMDescription.h"
#include "DescriptionOfCorrect.h"

PSMDescription::PSMDescription() :
    features(NULL), expMass(0.), calcMass(0.), scan(0),
    id_(""), peptide("") {
}

PSMDescription::PSMDescription(const std::string& pep) :
    features(NULL), expMass(0.), calcMass(0.), scan(0),
    id_(""), peptide(pep) {
}

PSMDescription::~PSMDescription() {}

void PSMDescription::deletePtr(PSMDescription* psm) {
  psm->deleteRetentionFeatures();
  delete psm;
  psm = NULL;
}

bool PSMDescription::isNotEnzymatic() {
  std::string peptideSeq = removePTMs(peptide);
  std::string peptideSeqNoFlanks = removeFlanks(peptide);
  return !(Enzyme::isEnzymatic(peptideSeq[0], peptideSeq[2])
      && Enzyme::isEnzymatic(peptideSeq[peptideSeq.size() - 3],
                             peptideSeq[peptideSeq.size() - 1])
      && Enzyme::countEnzymatic(peptideSeqNoFlanks) == 0);
}

std::string PSMDescription::removePTMs(const string& peptide) {
  std::string peptideSequence = peptide;
  peptideSequence = peptide.substr(2, peptide.size()- 4);
  for (unsigned int ix = 0; ix < peptideSequence.size(); ++ix) {
    if (peptideSequence[ix] == '[') {
      size_t posEnd = peptideSequence.substr(ix).find_first_of(']');
      if (posEnd == string::npos) {
        ostringstream temp;
        temp << "Error : Peptide sequence " << peptide << " contains an invalid modification" << endl;
        throw MyException(temp.str());
      } else {
        peptideSequence.erase(ix--, posEnd + 1);
      }
    }
  }
  return peptide.substr(0,1) + std::string(".") + peptideSequence + std::string(".") + peptide.substr(peptide.size() - 1,1);
}

void PSMDescription::printProteins(std::ostream& out) {
  std::vector<std::string>::const_iterator it = proteinIds.begin();
  for ( ; it != proteinIds.end(); ++it) {
    out << '\t' << *it;
  }
}
