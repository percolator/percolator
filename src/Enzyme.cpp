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

Enzyme* Enzyme::theEnzyme = NULL;

Enzyme* Enzyme::getEnzyme() {
  if (theEnzyme == NULL) {
    theEnzyme = new Trypsin(); // Use trypsin per default
  }
  return theEnzyme;
}

void Enzyme::destroy() {
  if (theEnzyme) {
    delete theEnzyme;
  }
  theEnzyme = NULL;
}

void Enzyme::setEnzyme(EnzymeType enz) {
  destroy();
  
  switch (enz) {
    case TRYPSINP:
      theEnzyme = new TrypsinP();
      return;
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
      theEnzyme = new Trypsin();
      return;
  }
}

void Enzyme::setEnzyme(std::string enzyme) {
  destroy();
  
  std::transform(enzyme.begin(), enzyme.end(), enzyme.begin(), ::tolower);
  if (enzyme == Chymotrypsin::getString()) {
    theEnzyme = new Chymotrypsin();
  } else if (enzyme == Thermolysin::getString()) {
    theEnzyme = new Thermolysin();
  } else if (enzyme == Proteinasek::getString()) {
    theEnzyme = new Proteinasek();
  } else if (enzyme == Pepsin::getString()) {
    theEnzyme = new Pepsin();
  } else if (enzyme == Elastase::getString()) {
    theEnzyme = new Elastase();
  } else if (enzyme == LysN::getString()) {
    theEnzyme = new LysN();
  } else if (enzyme == LysC::getString()) {
    theEnzyme = new LysC();
  } else if (enzyme == ArgC::getString()) {
    theEnzyme = new ArgC();
  } else if (enzyme == AspN::getString()) {
    theEnzyme = new AspN();
  } else if (enzyme == GluC::getString()) {
    theEnzyme = new GluC();
  } else if (enzyme == Enzyme::getString()) {
    theEnzyme = new Enzyme();
  } else if (enzyme == TrypsinP::getString()) {
    theEnzyme = new TrypsinP();
  } else if (enzyme == Trypsin::getString()) {
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

