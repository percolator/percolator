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

#include "Protein.h"

Protein::Protein(string namenew, double qnew, double qempnew, double pepnew, 
    double pnew, bool isdecoy_new, Protein::Peptide* __peptide)
	: name(namenew), q(qnew), qemp(qempnew), pep(pepnew), p(pnew),
	  isDecoy(isdecoy_new) {
  if (__peptide) peptides.push_back(__peptide);
}

Protein::~Protein() {
  for (unsigned i = 0; i < peptides.size(); i++)
    if (peptides[i]) delete peptides[i];
}

