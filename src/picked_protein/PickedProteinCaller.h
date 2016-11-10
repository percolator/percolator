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
  
  void initConstraints(PercolatorCrux::ENZYME_T enzyme, PercolatorCrux::DIGEST_T digestion, 
    int min_peptide_length, int max_peptide_length, int max_miscleavages);
  /* introductory message */
  std::string greeter() const;
  /* parse the command line arguments */
  bool parseOptions(int argc, char** argv);
  
  bool getProteinFragmentsAndDuplicates(
    std::map<std::string, std::string>& fragment_map,
    std::map<std::string, std::string>& duplicate_map,
    bool generateDecoys);
  
  void setFastaDatabase(const std::string& protein_db_file) {
    protein_db_file_ = protein_db_file;
  }
 private:
  PercolatorCrux::ENZYME_T enzyme_;
  PercolatorCrux::DIGEST_T digestion_;
  int min_peptide_length_, max_peptide_length_, max_miscleavages_;
  std::string decoyPattern_;
  
  std::string protein_db_file_, peptide_input_file_, protein_output_file_;
  
  bool getPeptideProteinMap(PercolatorCrux::Database& db, 
    PercolatorCrux::PeptideConstraint& peptide_constraint,
    std::map<std::string, std::vector<size_t> >& peptide_protein_map,
    bool generateDecoys);
  void addProteinToPeptideProteinMap(PercolatorCrux::Database& db, 
    size_t protein_idx, PercolatorCrux::PeptideConstraint& peptide_constraint,
    std::map<std::string, std::vector<size_t> >& peptide_protein_map,
    bool generateDecoys);
  
  bool getFragmentProteinMap(PercolatorCrux::Database& db, 
    PercolatorCrux::PeptideConstraint& peptide_constraint,
    std::map<std::string, std::vector<size_t> >& peptide_protein_map,
    std::map<size_t, std::vector<size_t> >& fragment_protein_map,
    std::map<size_t, size_t>& num_peptides_per_protein);
  void addProteinToFragmentProteinMap(PercolatorCrux::Database& db, 
    size_t protein_idx, PercolatorCrux::PeptideConstraint& peptide_constraint,
    std::map<std::string, std::vector<size_t> >& peptide_protein_map,
    std::map<size_t, std::vector<size_t> >& fragment_protein_map,
    std::map<size_t, size_t>& num_peptides_per_protein);
  
  bool getProteinFragmentsAndDuplicatesExtraDigest(PercolatorCrux::Database& db, 
    PercolatorCrux::PeptideConstraint& peptide_constraint,
    std::map<size_t, std::vector<size_t> >& fragment_protein_map,
    std::map<size_t, size_t>& num_peptides_per_protein,
    std::map<std::string, std::string>& fragment_map, 
    std::map<std::string, std::string>& duplicate_map);
  bool getProteinFragmentsAndDuplicates(PercolatorCrux::Database& db,
    std::map<size_t, std::vector<size_t> >& fragment_protein_map,
    std::map<size_t, size_t>& num_peptides_per_protein,
    std::map<std::string, std::string>& fragment_map, 
    std::map<std::string, std::string>& duplicate_map);
};

#endif /* PICKED_PROTEIN_PICKED_PROTEIN_CALLER_H_ */
