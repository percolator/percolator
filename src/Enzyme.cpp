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

Enzyme* Enzyme::createEnzyme(EnzymeType enz) {
  Enzyme* theEnzyme;
  switch (enz) {
    case TRYPSINP:
      theEnzyme = new TrypsinP();
      break;
    case CHYMOTRYPSIN:
      theEnzyme = new Chymotrypsin();
      break;
    case THERMOLYSIN:
      theEnzyme = new Thermolysin();
      break;
    case PROTEINASEK:
      theEnzyme = new Proteinasek();
      break;
    case PEPSIN:
      theEnzyme = new Pepsin();
      break;
    case ELASTASE:
      theEnzyme = new Elastase();
      break;
    case LYSN:
      theEnzyme = new LysN();
      break;
    case LYSC:
      theEnzyme = new LysC();
      break;
    case ARGC:
      theEnzyme = new ArgC();
      break;
    case ASPN:
      theEnzyme = new AspN();
      break;
    case GLUC:
      theEnzyme = new GluC();
      break;
    case NO_ENZYME:
      theEnzyme = new NoEnzyme();
      break;
    case TRYPSIN:
    default:
      theEnzyme = new Trypsin();
      break;
  }
  return theEnzyme;
}

Enzyme* Enzyme::createEnzyme(std::string enzyme) {
  Enzyme* theEnzyme;
  
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
  } else if (enzyme == NoEnzyme::getString()) {
    theEnzyme = new NoEnzyme();
  } else if (enzyme == TrypsinP::getString()) {
    theEnzyme = new TrypsinP();
  } else if (enzyme == Trypsin::getString()) {
    theEnzyme = new Trypsin();
  }
  else {
    ostringstream temp;
    temp << "The selected enzyme \"" << enzyme << "\" has no corresponding class" << std::endl;
    throw MyException(temp.str());
  }
  
  return theEnzyme;
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

