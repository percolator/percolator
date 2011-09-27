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
    case ELASTASE:
      theEnzyme = new Elastase();
      return;
    case LYSN:
      theEnzyme = new LysN();
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
  // boost::iequals: case insensitive string comparison
  //TODO NOT VERY PORTABLE, USING stricmp INSTEAD
//   if (boost::iequals(enzyme,Chymotrypsin::getString())) {
//       theEnzyme = new Chymotrypsin();
//   } else if (boost::iequals(enzyme,Elastase::getString())) {
//       theEnzyme = new Elastase();
//   } else if (boost::iequals(enzyme,Enzyme::getString())) {
//       theEnzyme = new Enzyme();
//   } else if (boost::iequals(enzyme,Trypsin::getString())) {
//       theEnzyme = new Trypsin();
  if (stricmp(enzyme,Chymotrypsin::getString())) {
      theEnzyme = new Chymotrypsin();
  } else if (stricmp(enzyme,Elastase::getString())) {
      theEnzyme = new Elastase();
  } else if (stricmp(enzyme,Enzyme::getString())) {
      theEnzyme = new Enzyme();
  } else if (stricmp(enzyme,Trypsin::getString())) {
      theEnzyme = new Trypsin();
  } else {
    std::cerr << "The selected enzyme have no corresponding class" << std::endl;
    std::exit(-1);
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

bool CompareCaseInsensitive(std::string strFirst, std::string strSecond)
{
  // Convert both strings to upper case by transfrom() before compare.
  std::transform(strFirst.begin(), strFirst.end(), strFirst.begin(), toupper);
  std::transform(strSecond.begin(), strSecond.end(), strSecond.begin(), toupper);
  if(strFirst == strSecond) return true; else return false;
}
