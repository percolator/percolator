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
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <algorithm>

#include "Database.h"
#include "DatabaseProteinIterator.h"
#include "PeptideConstraint.h"
#include "ProteinPeptideIterator.h"
#include "Protein.h"
using namespace std;

int main(int argc, char** argv) {
  
  // now create a database
  Database* database = new Database(argv[1], false);
  
  if(!database->parse()){
    std::cerr << "Failed to parse database, cannot create new index for " << argv[1] << std::endl;
    return EXIT_FAILURE;
  } else {
    std::cerr << "Read " << database->getNumProteins() << " proteins." << std::endl;
    std::cerr << database->getProteinAtIdx(5)->getSequencePointer() << std::endl;
  }
  
  // set a new protein iterator
  DatabaseProteinIterator* database_protein_iterator = new DatabaseProteinIterator(database);
  if(database_protein_iterator == NULL){
    std::cerr << "Could not create protein iterator." << std::endl;
    return EXIT_FAILURE;
  }

  // set peptide constraint
  
  PeptideConstraint* peptide_constraint = new PeptideConstraint(TRYPSIN, FULL_DIGEST, 6, 50, 0);
  
  Crux::Protein* protein = NULL;
  
  size_t protein_idx = 0;
  std::map<std::string, std::vector<size_t> > peptide_protein_map;
  // check if there are any proteins to create peptides from
  while(database_protein_iterator->hasNext()){

    protein = database_protein_iterator->next();
  
    //std::cerr << protein->getSequencePointer() << std::endl;
    
    // set new protein peptide iterator
    ProteinPeptideIterator* cur_protein_peptide_iterator =
      new ProteinPeptideIterator(protein, peptide_constraint);
    
    // if first protein does not contain a match peptide, reinitailize
    while(cur_protein_peptide_iterator->hasNext()) {
      Crux::Peptide* peptide = cur_protein_peptide_iterator->next();
      char* c_seq = peptide->getSequence();
      std::string sequence(c_seq);
      peptide_protein_map[sequence].push_back(protein_idx);
      //std::cerr << sequence << std::endl; 
      //std::cerr << peptide->getSequence() << std::endl;
      std::free(c_seq);
      delete peptide;
    }
    
    ++protein_idx;
    // free old iterator
    delete (cur_protein_peptide_iterator);
  }
  
  delete database_protein_iterator;
  
  std::cerr << peptide_protein_map["VRPLAR"].size() << std::endl;
  std::cerr << peptide_protein_map["EHHEHASAPLLPPPPTSALSSIASTTAASSAHAK"].size() << std::endl;
  std::cerr << database->getProteinAtIdx(peptide_protein_map["EHHEHASAPLLPPPPTSALSSIASTTAASSAHAK"][0])->getSequencePointer() << std::endl;
  
  std::map<size_t, std::vector<size_t> > fragment_protein_map;
  std::map<size_t, size_t> num_peptides_per_protein;
  
  for (size_t i = 0; i < database->getNumProteins(); ++i) {
    protein = database->getProteinAtIdx(i);
    
    // set new protein peptide iterator
    ProteinPeptideIterator* cur_protein_peptide_iterator =
      new ProteinPeptideIterator(protein, peptide_constraint);
    
    bool is_first = true;
    size_t num_sequences = 0;
    std::vector<size_t> protein_idx_intersection;
    // if first protein does not contain a match peptide, reinitailize
    while(cur_protein_peptide_iterator->hasNext()) {
      Crux::Peptide* peptide = cur_protein_peptide_iterator->next();
      char* c_seq = peptide->getSequence();
      std::string sequence(c_seq);
      ++num_sequences;
      
      if (is_first) {
        protein_idx_intersection = peptide_protein_map[sequence];
        is_first = false;
      } else {
        std::vector<size_t> newprotein_idx_intersection(protein_idx_intersection.size());
        std::vector<size_t>::iterator it = std::set_intersection(
            protein_idx_intersection.begin(), protein_idx_intersection.end(), 
            peptide_protein_map[sequence].begin(), peptide_protein_map[sequence].end(), 
            newprotein_idx_intersection.begin());
        newprotein_idx_intersection.resize(it - newprotein_idx_intersection.begin());
        protein_idx_intersection = newprotein_idx_intersection;
      }
      
      //std::cerr << sequence << std::endl; 
      //std::cerr << peptide->getSequence() << std::endl;
      
      std::free(c_seq);
      delete peptide;
      
      if (protein_idx_intersection.size() < 2) break;
    }
    
    if (protein_idx_intersection.size() > 1) {
      num_peptides_per_protein[i] = num_sequences;
      std::sort(protein_idx_intersection.begin(), protein_idx_intersection.end()); // sort in descending order
      
      size_t j = protein_idx_intersection.size() - 1;
      if (protein_idx_intersection[j] == i) --j;
      
      if (fragment_protein_map.find(i) != fragment_protein_map.end()) {
        fragment_protein_map[protein_idx_intersection[j]].insert(
            fragment_protein_map[protein_idx_intersection[j]].end(),
            fragment_protein_map[i].begin(), fragment_protein_map[i].end());
        fragment_protein_map.erase(i);
      }
      fragment_protein_map[protein_idx_intersection[j]].push_back(i);
    }
    
    // free old iterator
    delete (cur_protein_peptide_iterator);
  }
  
  std::map<std::string, std::string> fragment_map, duplicate_map;
  std::map<size_t, std::vector<size_t> >::iterator it;
  for (it = fragment_protein_map.begin(); it != fragment_protein_map.end(); ++it) {
    size_t i = it->first;
    
    char* c_this_id = database->getProteinAtIdx(i)->getId();
    std::string this_protein_id(c_this_id);
    std::free(c_this_id);
    
    for (std::vector<size_t>::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
      size_t j = *it2;
      
      if (i != j) {
        char* c_that_id = database->getProteinAtIdx(j)->getId();
        std::string that_protein_id(c_that_id);
        std::free(c_that_id);
        
        if (num_peptides_per_protein[i] == num_peptides_per_protein[j]) {
          duplicate_map[that_protein_id] = this_protein_id;
          //std::cerr << "Duplicate: " << that_protein_id << " " << this_protein_id << std::endl;
        } else {
          fragment_map[that_protein_id] = this_protein_id;
          //std::cerr << "Fragment: " << that_protein_id << " " << this_protein_id << std::endl;
        }
      }
    }
  }

  delete peptide_constraint;
  delete database;
  
  return EXIT_SUCCESS;
}
