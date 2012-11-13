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

#include "Enzyme.h"

Enzyme* Enzyme::theEnzyme = NULL;

Enzyme* Enzyme::getEnzyme() {
  if (theEnzyme == NULL) {
    theEnzyme = new Trypsin(); // Use trypsin per default
  }
  return theEnzyme;
}

void Enzyme::setEnzyme(EnzymeType enz) {
  if (theEnzyme) {
    delete theEnzyme;
  }
  theEnzyme = NULL;
  switch (enz) {
    case CHYMOTRYPSIN:
      theEnzyme = new Chymotrypsin();
      return;
    case THERMOLYSIN:
      theEnzyme = new Thermolysin();
      return;
    case PROTEINASEK:
      theEnzyme = new Proteinasek();
      return;
    case PEPSIN:
      theEnzyme = new Pepsin();
      return;
    case ELASTASE:
      theEnzyme = new Elastase();
      return;
    case LYSN:
      theEnzyme = new LysN();
      return;
    case LYSC:
      theEnzyme = new LysC();
      return;
    case ARGC:
      theEnzyme = new ArgC();
      return;
    case ASPN:
      theEnzyme = new AspN();
      return;
    case GLUC:
      theEnzyme = new GluC();
      return;
    case NO_ENZYME:
      theEnzyme = new Enzyme();
      return;
    case TRYPSIN:
    default:
      theEnzyme = new Trypsin;
      return;
  }
}

void Enzyme::setEnzyme(std::string enzyme) {
  if (theEnzyme) {
    delete theEnzyme;
  }
  theEnzyme = NULL;
  
  if (boost::iequals(enzyme,Chymotrypsin::getString())) {
      theEnzyme = new Chymotrypsin();
  } else if (boost::iequals(enzyme,Thermolysin::getString())) {
      theEnzyme = new Thermolysin();
  } else if (boost::iequals(enzyme,Proteinasek::getString())) {
      theEnzyme = new Proteinasek();
  } else if (boost::iequals(enzyme,Pepsin::getString())) {
      theEnzyme = new Pepsin();
  } else if (boost::iequals(enzyme,Elastase::getString())) {
      theEnzyme = new Elastase();
  } else if (boost::iequals(enzyme,LysN::getString())) {
      theEnzyme = new LysN();
  } else if (boost::iequals(enzyme,LysC::getString())) {
      theEnzyme = new LysC();
  } else if (boost::iequals(enzyme,ArgC::getString())) {
      theEnzyme = new ArgC();
  } else if (boost::iequals(enzyme,AspN::getString())) {
      theEnzyme = new AspN();
  } else if (boost::iequals(enzyme,GluC::getString())) {
      theEnzyme = new GluC();
  } else if (boost::iequals(enzyme,Enzyme::getString())) {
      theEnzyme = new Enzyme();
  } else if (boost::iequals(enzyme,Trypsin::getString())) {
      theEnzyme = new Trypsin();
  }
  else {
    ostringstream temp;
    temp << "The selected enzyme have no corresponding class" << std::endl;
    throw MyException(temp.str());
  }
}

size_t Enzyme::countEnzymatic(std::string& peptide) {
  size_t count = 0;
  for (size_t ix = 1; ix < peptide.size(); ++ix) {
    if (isEnzymatic(peptide[ix - 1], peptide[ix])) {
      ++count;
    }
  }
  return count;
}

