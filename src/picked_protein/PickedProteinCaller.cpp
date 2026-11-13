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

#include "PickedProteinCaller.h"

#include <algorithm>
#include <fstream>
#include <iostream>

#include "Globals.h"

using namespace std;
using PickedProtein::DigestConstraint;
using PickedProtein::DigestType;
using PickedProtein::FastaProtein;
using PickedProtein::PeptideIndex;

PickedProteinCaller::PickedProteinCaller()
    : enzyme_(Enzyme::createEnzyme(Enzyme::TRYPSIN)),
      digestion_(DigestType::Full), min_peptide_length_(6),
      max_peptide_length_(50), max_miscleavages_(0),
      decoyPattern_("decoy_"), fasta_has_decoys_(false),
      database_read_(false) {}

PickedProteinCaller::~PickedProteinCaller() {}

void PickedProteinCaller::initConstraints(Enzyme::EnzymeType enzyme,
    DigestType digestion, int min_peptide_length, int max_peptide_length,
    int max_miscleavages) {
  enzyme_ = Enzyme::createEnzyme(enzyme);
  digestion_ = digestion;
  min_peptide_length_ = min_peptide_length;
  max_peptide_length_ = max_peptide_length;
  max_miscleavages_ = max_miscleavages;
}

//!
//! fills \p fragment_protein_map which contains groups of proteins with same or
//! subset peptides
//!
//! @param[in] index peptide index of the digest
//! @param[in] protein_idx index of the protein in the database
//! @param[in] num_peptides_per_protein helps to find out which proteins have
//!   identical sets and which ones are proper subsets.
//! @param[out] fragment_protein_map groups of proteins with same or subset peptides
//!
void PickedProteinCaller::findFragmentProteins(const PeptideIndex& index,
    size_t protein_idx, std::map<size_t, size_t>& num_peptides_per_protein,
    std::map<size_t, std::vector<size_t> >& fragment_protein_map) {
  // find all proteins of which the current protein's peptides are a subset (possibly identical)
  std::vector<size_t> protein_idx_intersection;
  index.proteinsContainingAllPeptides(protein_idx, protein_idx_intersection);

  // if there are still proteins left in the intersection, it means that the
  // current protein is a subset of at least one another protein
  if (protein_idx_intersection.size() > 1) {
    addToFragmentProteinMap(protein_idx, protein_idx_intersection,
        num_peptides_per_protein, fragment_protein_map);
  }
}

void PickedProteinCaller::addToFragmentProteinMap(
    const size_t protein_idx, std::vector<size_t>& protein_idx_intersection,
    std::map<size_t, size_t>& num_peptides_per_protein,
    std::map<size_t, std::vector<size_t> >& fragment_protein_map) {
  // sort in descending order such that the first element will be the last one
  // to be visited
  std::sort(protein_idx_intersection.rbegin(), protein_idx_intersection.rend());
  fragment_protein_map[protein_idx_intersection[0]].push_back(protein_idx);

  if (protein_idx_intersection[0] == protein_idx) {
    size_t largestProteinIdx = 0;
    size_t largestProteinNumPeptides = 0;
    for (std::vector<size_t>::iterator it = protein_idx_intersection.begin(); it != protein_idx_intersection.end(); ++it) {
      if (num_peptides_per_protein[*it] > largestProteinNumPeptides) {
        largestProteinIdx = *it;
        largestProteinNumPeptides = num_peptides_per_protein[*it];
      }
    }

    if (largestProteinIdx != protein_idx) {
      // if there are other proteins that are subsets of the current protein,
      // move them together with the current protein
      fragment_protein_map[largestProteinIdx].insert(
          fragment_protein_map[largestProteinIdx].end(),
          fragment_protein_map[protein_idx].begin(),
          fragment_protein_map[protein_idx].end());
      fragment_protein_map.erase(protein_idx);
    }
  }
}

void PickedProteinCaller::findFragmentsAndDuplicates(const Proteins& proteins,
    std::map<size_t, std::vector<size_t> >& fragment_protein_map,
    std::map<size_t, size_t>& num_peptides_per_protein,
    std::map<std::string, std::string>& fragment_map,
    std::map<std::string, std::string>& duplicate_map) {
  std::map<size_t, std::vector<size_t> >::iterator it;
  for (it = fragment_protein_map.begin(); it != fragment_protein_map.end(); ++it) {
    size_t i = it->first;
    const std::string& this_protein_id = proteins[i].id;

    for (std::vector<size_t>::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
      size_t j = *it2;

      if (i != j) {
        const std::string& that_protein_id = proteins[j].id;

        if (num_peptides_per_protein[i] == num_peptides_per_protein[j]) {
          duplicate_map[that_protein_id] = this_protein_id;
        } else {
          fragment_map[that_protein_id] = this_protein_id;
        }
      }
    }
  }
}

//! creates a local peptide->protein map for each candidate protein group that
//! which was based on a full digest with max_len <= 50 and max_misclv <= 2
void PickedProteinCaller::findFragmentsAndDuplicatesExtraDigest(
    const Proteins& proteins, const DigestConstraint& constraint,
    std::map<size_t, std::vector<size_t> >& fragment_protein_map,
    std::map<std::string, std::string>& fragment_map,
    std::map<std::string, std::string>& duplicate_map) {
  std::map<size_t, std::vector<size_t> >::iterator it;
  for (it = fragment_protein_map.begin(); it != fragment_protein_map.end(); ++it) {
    size_t i = it->first;
    if (std::find(it->second.begin(), it->second.end(), i) == it->second.end()) {
      it->second.push_back(i);
    }
    std::sort(it->second.begin(), it->second.end());

    PeptideIndex index(proteins, it->second, *enzyme_, constraint);
    std::map<size_t, size_t> num_peptides_per_protein_local;
    for (std::vector<size_t>::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
      num_peptides_per_protein_local[*it2] = index.numPeptides(*it2);
    }

    std::map<size_t, std::vector<size_t> > fragment_protein_map_local;
    for (std::vector<size_t>::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
      findFragmentProteins(index, *it2, num_peptides_per_protein_local,
          fragment_protein_map_local);
    }

    findFragmentsAndDuplicates(proteins, fragment_protein_map_local,
        num_peptides_per_protein_local, fragment_map, duplicate_map);
  }
}

void PickedProteinCaller::findFragmentsAndDuplicatesNonSpecificDigest(
    const Proteins& proteins,
    std::map<size_t, std::vector<size_t> >& fragment_protein_map,
    std::map<std::string, std::string>& fragment_map,
    std::map<std::string, std::string>& duplicate_map) {
  std::map<size_t, std::vector<size_t> >::iterator it;
  for (it = fragment_protein_map.begin(); it != fragment_protein_map.end(); ++it) {
    size_t i = it->first;
    if (std::find(it->second.begin(), it->second.end(), i) == it->second.end()) {
      it->second.push_back(i);
    }
    std::sort(it->second.begin(), it->second.end());

    std::map<size_t, size_t> num_peptides_per_protein_local;
    for (std::vector<size_t>::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
      num_peptides_per_protein_local[*it2] = proteins[*it2].sequence.size();
    }

    // In a non-specific digest the only possibility for one protein to be subset
    // of another is if the entire string is contained (except for some very
    // unlikely cases where a region longer than max_len is repeated more than twice)
    std::map<size_t, std::vector<size_t> > fragment_protein_map_local;
    for (std::vector<size_t>::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
      size_t protein_idx = *it2;
      const std::string& sequence = proteins[protein_idx].sequence;

      std::vector<size_t> protein_idx_intersection;
      for (std::vector<size_t>::iterator it3 = it->second.begin(); it3 != it->second.end(); ++it3) {
        if (proteins[*it3].sequence.find(sequence) != std::string::npos) {
          protein_idx_intersection.push_back(*it3);
        }
      }

      if (protein_idx_intersection.size() > 1) {
        addToFragmentProteinMap(protein_idx, protein_idx_intersection,
            num_peptides_per_protein_local, fragment_protein_map_local);
      }
    }

    findFragmentsAndDuplicates(proteins, fragment_protein_map_local,
        num_peptides_per_protein_local, fragment_map, duplicate_map);
  }
}

bool PickedProteinCaller::getProteinFragmentsAndDuplicates(
    std::map<std::string, std::string>& fragment_map,
    std::map<std::string, std::string>& duplicate_map,
    bool reverseProteinSeqs) {
  time_t startTime;
  time(&startTime);
  clock_t startClock = clock();

  if (!database_read_) {
    proteins_.clear();
    std::ifstream fasta(protein_db_file_.c_str(), std::ios::binary);
    if (!fasta.is_open() || !PickedProtein::readFasta(fasta, proteins_)) {
      std::cerr << "Failed to parse database, cannot create index for "
                << protein_db_file_ << std::endl;
      proteins_.clear();
      return true;
    }
    database_read_ = true;
  }

  if (reverseProteinSeqs) {
    // the decoy database, in place: reversed sequences, identifiers with the
    // decoy prefix; a later call reads the file again
    for (FastaProtein& protein : proteins_) {
      std::reverse(protein.sequence.begin(), protein.sequence.end());
      // MT: the crux interface will change the protein identifier. If we are
      // not inside the crux environment we do this separately here.
      if (protein.id.compare(0, decoyPattern_.size(), decoyPattern_) != 0) {
        protein.id = decoyPattern_ + protein.id;
      }
    }
    database_read_ = false;
  } else if (!fasta_has_decoys_) {
    for (const FastaProtein& protein : proteins_) {
      if (protein.id.compare(0, decoyPattern_.size(), decoyPattern_) == 0) {
        fasta_has_decoys_ = true;
        break;
      }
    }
  }
  const Proteins& proteins = proteins_;

  if (VERB > 3) reportProgress("Creating database", startTime, startClock);

  // First do a "basic" digest to get candidates for protein grouping*
  //
  // * there are some rare cases in which fragment proteins have
  // distinct subsets in a full digest, but form a proper subset in a
  // partial or non-enzymatic digest (e.g. if one would add a single
  // amino acid to the beginning of a protein sequence). It's too hard
  // to check for these at the moment...
  //
  std::map<size_t, std::vector<size_t> > fragment_protein_map;
  std::map<size_t, size_t> num_peptides_per_protein;
  {
    const DigestConstraint basic_constraint{DigestType::Full,
        min_peptide_length_, (std::min)(50, max_peptide_length_),
        (std::min)(2, max_miscleavages_)};
    const PeptideIndex index(proteins, std::vector<size_t>(), *enzyme_,
                             basic_constraint);
    for (size_t protein_idx = 0; protein_idx < proteins.size(); ++protein_idx) {
      num_peptides_per_protein.emplace_hint(num_peptides_per_protein.end(),
          protein_idx, index.numPeptides(protein_idx));
    }

    if (VERB > 3) {
      reportProgress("Creating protein peptide map", startTime, startClock);
    }

    // Find all proteins whose peptides form a subset (possibly identical)
    // of another protein
    for (size_t protein_idx = 0; protein_idx < proteins.size(); ++protein_idx) {
      findFragmentProteins(index, protein_idx, num_peptides_per_protein,
          fragment_protein_map);
    }
  }

  if (VERB > 3) {
    reportProgress("Creating fragment protein map", startTime, startClock);
  }

  // If the actual digest was more permissive than the basic digest (see above),
  // check validity of each candidate protein group
  if (digestion_ == DigestType::Full && max_miscleavages_ <= 2
          && max_peptide_length_ <= 50) {
    findFragmentsAndDuplicates(proteins, fragment_protein_map,
                  num_peptides_per_protein, fragment_map, duplicate_map);
  } else if (digestion_ != DigestType::NonSpecific) {
    const DigestConstraint extra_constraint{digestion_, min_peptide_length_,
        max_peptide_length_, max_miscleavages_};
    findFragmentsAndDuplicatesExtraDigest(proteins, extra_constraint,
                  fragment_protein_map, fragment_map, duplicate_map);
  } else {
    findFragmentsAndDuplicatesNonSpecificDigest(proteins,
                  fragment_protein_map, fragment_map, duplicate_map);
  }

  if (VERB > 2) {
    reportProgress("Resolving protein fragments and duplicates",
                      startTime, startClock);
  }

  return false;
}

void PickedProteinCaller::reportProgress(const std::string& msg,
    const time_t& startTime, const clock_t& startClock) {
  time_t procStart;
  clock_t procStartClock = clock();
  time(&procStart);
  double diff = difftime(procStart, startTime);
  if (VERB > 3) cerr << msg << " took "
      << ((double)(procStartClock - startClock)) / (double)CLOCKS_PER_SEC
      << " cpu seconds or " << diff << " seconds wall time" << endl;
}
