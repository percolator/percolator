/*******************************************************************************
 Copyright 2006-2010 Lukas Kall <lukas.kall@scilifelab.se>

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
/*
 * @ Created by Matthew The
 * Sep, 2015
 */
/*
 * PickedProteinCaller finds the fragment and duplicate proteins of a FASTA
 * database: proteins whose digest is a subset of, or the same as, the digest
 * of another protein.
 */

#ifndef PICKED_PROTEIN_PICKED_PROTEIN_CALLER_H_
#define PICKED_PROTEIN_PICKED_PROTEIN_CALLER_H_

#include <ctime>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "Enzyme.h"
#include "ProteinDigest.h"

class PickedProteinCaller {
 public:
  PickedProteinCaller();
  ~PickedProteinCaller();

  void initConstraints(Enzyme::EnzymeType enzyme,
                       PickedProtein::DigestType digestion,
                       int min_peptide_length, int max_peptide_length,
                       int max_miscleavages);

  bool fastaHasDecoys() const { return fasta_has_decoys_; }

  void setFastaDatabase(const std::string& protein_db_file,
                        const std::string& decoyPattern) {
    protein_db_file_ = protein_db_file;
    decoyPattern_ = decoyPattern;
    proteins_.clear();
    database_read_ = false;
  }

  /**
   * Adds the fragment and duplicate proteins of the database to the maps,
   * each mapped to the protein that represents it. With reverseProteinSeqs,
   * the protein sequences are reversed and the identifiers get the decoy
   * prefix first. The FASTA file is read on the first call and kept for the
   * next one, except after a call with reverseProteinSeqs, which reverses
   * the kept proteins in place. Returns true if the file could not be read.
   */
  bool getProteinFragmentsAndDuplicates(
      std::map<std::string, std::string>& fragment_map,
      std::map<std::string, std::string>& duplicate_map,
      bool reverseProteinSeqs);

 private:
  typedef std::vector<PickedProtein::FastaProtein> Proteins;

  std::unique_ptr<Enzyme> enzyme_;
  PickedProtein::DigestType digestion_;
  int min_peptide_length_, max_peptide_length_, max_miscleavages_;
  std::string decoyPattern_;
  bool fasta_has_decoys_;

  std::string protein_db_file_;
  Proteins proteins_;
  bool database_read_;

  void findFragmentProteins(
      const PickedProtein::PeptideIndex& index, size_t protein_idx,
      std::map<size_t, size_t>& num_peptides_per_protein,
      std::map<size_t, std::vector<size_t> >& fragment_protein_map);
  void addToFragmentProteinMap(
      const size_t protein_idx, std::vector<size_t>& protein_idx_intersection,
      std::map<size_t, size_t>& num_peptides_per_protein,
      std::map<size_t, std::vector<size_t> >& fragment_protein_map);

  void findFragmentsAndDuplicates(
      const Proteins& proteins,
      std::map<size_t, std::vector<size_t> >& fragment_protein_map,
      std::map<size_t, size_t>& num_peptides_per_protein,
      std::map<std::string, std::string>& fragment_map,
      std::map<std::string, std::string>& duplicate_map);
  void findFragmentsAndDuplicatesExtraDigest(
      const Proteins& proteins,
      const PickedProtein::DigestConstraint& constraint,
      std::map<size_t, std::vector<size_t> >& fragment_protein_map,
      std::map<std::string, std::string>& fragment_map,
      std::map<std::string, std::string>& duplicate_map);
  void findFragmentsAndDuplicatesNonSpecificDigest(
      const Proteins& proteins,
      std::map<size_t, std::vector<size_t> >& fragment_protein_map,
      std::map<std::string, std::string>& fragment_map,
      std::map<std::string, std::string>& duplicate_map);

  void reportProgress(const std::string& msg,
    const time_t& startTime, const clock_t& startClock);
};

#endif /* PICKED_PROTEIN_PICKED_PROTEIN_CALLER_H_ */
