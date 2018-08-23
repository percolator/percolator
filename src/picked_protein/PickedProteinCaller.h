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
 * This file stores the class PickedProteinCaller which is managing the interface with the user
 */

#ifndef PICKED_PROTEIN_PICKED_PROTEIN_CALLER_H_
#define PICKED_PROTEIN_PICKED_PROTEIN_CALLER_H_

#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <iostream>
#include <algorithm>

#include "Database.h"
#include "PeptideConstraint.h"
#include "ProteinPeptideIterator.h"
#include "Protein.h"

class PickedProteinCaller{
 public:
  PickedProteinCaller();
  ~PickedProteinCaller();
  
  void initConstraints(PercolatorCrux::ENZYME_T enzyme, 
      PercolatorCrux::DIGEST_T digestion, int min_peptide_length, 
      int max_peptide_length, int max_miscleavages);
  
  std::string greeter() const;
  
  bool parseOptions(int argc, char** argv);
  
  bool fastaHasDecoys() const { return fasta_has_decoys_; }
  
  void setFastaDatabase(const std::string& protein_db_file, 
                        const std::string& decoyPattern) {
    protein_db_file_ = protein_db_file;
    decoyPattern_ = decoyPattern;
  }
  
  bool getProteinFragmentsAndDuplicates(
      std::map<std::string, std::string>& fragment_map,
      std::map<std::string, std::string>& duplicate_map,
      bool reverseProteinSeqs);
  
 private:
  PercolatorCrux::ENZYME_T enzyme_;
  PercolatorCrux::DIGEST_T digestion_;
  int min_peptide_length_, max_peptide_length_, max_miscleavages_;
  std::string decoyPattern_;
  bool fasta_has_decoys_;
  
  std::string protein_db_file_, peptide_input_file_, protein_output_file_;
  
  void addToPeptideProteinMap(PercolatorCrux::Database& db, 
    size_t protein_idx, PercolatorCrux::PeptideConstraint& peptide_constraint,
    std::map<std::string, std::vector<size_t> >& peptide_protein_map,
    std::map<size_t, size_t>& num_peptides_per_protein,
    bool reverseProteinSeqs);
  
  void findFragmentProteins(PercolatorCrux::Database& db, 
    size_t protein_idx, PercolatorCrux::PeptideConstraint& peptide_constraint,
    std::map<std::string, std::vector<size_t> >& peptide_protein_map,
    std::map<size_t, size_t>& num_peptides_per_protein,
    std::map<size_t, std::vector<size_t> >& fragment_protein_map);
  void addToFragmentProteinMap(
    const size_t protein_idx, std::vector<size_t>& protein_idx_intersection,
    std::map<size_t, size_t>& num_peptides_per_protein,
    std::map<size_t, std::vector<size_t> >& fragment_protein_map);
  
  void findFragmentsAndDuplicates(
    PercolatorCrux::Database& db,
    std::map<size_t, std::vector<size_t> >& fragment_protein_map,
    std::map<size_t, size_t>& num_peptides_per_protein,
    std::map<std::string, std::string>& fragment_map, 
    std::map<std::string, std::string>& duplicate_map);
  void findFragmentsAndDuplicatesExtraDigest(
    PercolatorCrux::Database& db, 
    PercolatorCrux::PeptideConstraint& peptide_constraint,
    std::map<size_t, std::vector<size_t> >& fragment_protein_map,
    std::map<std::string, std::string>& fragment_map, 
    std::map<std::string, std::string>& duplicate_map);
  void findFragmentsAndDuplicatesNonSpecificDigest(
    PercolatorCrux::Database& db, 
    std::map<size_t, std::vector<size_t> >& fragment_protein_map,
    std::map<std::string, std::string>& fragment_map, 
    std::map<std::string, std::string>& duplicate_map);
  
  void reportProgress(const std::string& msg,
    const time_t& startTime, const clock_t& startClock);
};

#endif /* PICKED_PROTEIN_PICKED_PROTEIN_CALLER_H_ */
