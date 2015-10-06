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
 * This file stores the class FisherCaller which is managing the interface with the user
 */

#ifndef FISHER_FISHERCALLER_H_
#define FISHER_FISHERCALLER_H_

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <algorithm>

#include "Database.h"
#include "DatabaseProteinIterator.h"
#include "PeptideConstraint.h"
#include "ProteinPeptideIterator.h"
#include "Protein.h"

class FisherCaller{
 public:
  FisherCaller();
  FisherCaller(PeptideConstraint* peptide_constraint);
  ~FisherCaller();
  
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
  Database* database_;
  PeptideConstraint* peptide_constraint_;
  
  std::string protein_db_file_, peptide_input_file_, protein_output_file_;
  
  bool getPeptideProteinMap(
    std::map<std::string, std::vector<size_t> >& peptide_protein_map,
    bool generateDecoys);
  
  bool getFragmentProteinMap(
    std::map<std::string, std::vector<size_t> >& peptide_protein_map,
    std::map<size_t, std::vector<size_t> >& fragment_protein_map,
    std::map<size_t, size_t>& num_peptides_per_protein);
  
  bool getProteinFragmentsAndDuplicates(
    std::map<size_t, std::vector<size_t> >& fragment_protein_map,
    std::map<size_t, size_t>& num_peptides_per_protein,
    std::map<std::string, std::string>& fragment_map, 
    std::map<std::string, std::string>& duplicate_map);
};

#endif /* FISHER_FISHERCALLER_H_ */
