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

int mystricmp(const char* str1, const char* str2)
{
    if (str1 == str2)
      return 0;
    else if (str1 == NULL)
      return -1;
    else if (str2 == NULL)
      return 1;
    else {
      while (tolower(*str1) == tolower(*str2) && *str1 != 0 && *str2 != 0)
    {
	++str1;
	++str2;
    }
    if (*str1 < *str2)
	return -1;
    else if (*str1 > *str2)
	return 1;
    else
	return 0;
    }
}

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

  //TODO NOT VERY PORTABLE, USING mystricmp INSTEAD
  
//   if (boost::iequals(enzyme,Chymotrypsin::getString())) {
//       theEnzyme = new Chymotrypsin();
//   } else if (boost::iequals(enzyme,Elastase::getString())) {
//       theEnzyme = new Elastase();
//   } else if (boost::iequals(enzyme,Enzyme::getString())) {
//       theEnzyme = new Enzyme();
//   } else if (boost::iequals(enzyme,Trypsin::getString())) {
//       theEnzyme = new Trypsin();
  
  if (mystricmp(enzyme.c_str(),Chymotrypsin::getString().c_str())) {
      theEnzyme = new Chymotrypsin();
  } else if (mystricmp(enzyme.c_str(),Elastase::getString().c_str())) {
      theEnzyme = new Elastase();
  } else if (mystricmp(enzyme.c_str(),Enzyme::getString().c_str())) {
      theEnzyme = new Enzyme();
  } else if (mystricmp(enzyme.c_str(),Trypsin::getString().c_str())) {
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

