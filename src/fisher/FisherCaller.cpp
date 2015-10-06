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

#include "FisherCaller.h"
#include "Option.h"
#include "Globals.h"
using namespace std;

FisherCaller::FisherCaller(PeptideConstraint* peptide_constraint) : 
    database_(NULL), peptide_constraint_(peptide_constraint) {}

FisherCaller::FisherCaller() : database_(NULL) {
  peptide_constraint_ = new PeptideConstraint(TRYPSIN, FULL_DIGEST, 10, 30, 0);
}

FisherCaller::~FisherCaller() {
  if (database_ != NULL) {
    delete database_;
  }
  if (peptide_constraint_ != NULL) {
    delete peptide_constraint_;
  }
}

/* introductory message */
string FisherCaller::greeter() const {
  ostringstream oss;
  oss << "Fisher version " << VERSION << ", ";
  oss << "Build Date " << __DATE__ << " " << __TIME__ << endl;
  oss << "Distributed under MIT License" << endl;
  oss << "Written by Lukas Kall (lukas.kall@scilifelab.se) "
         "and Matthew The (matthew.the@scilifelab.se)" << endl;
  oss << "Usage:" << endl;
  oss << "   fisher [options]" << endl << endl;
  return oss.str();
}

/* parse the command line arguments */
bool FisherCaller::parseOptions(int argc, char** argv) {
  ostringstream intro;
  intro << greeter() << endl << "Usage:" << endl;
  intro << "   fisher -i \"percolator_peptide_output\" -d \"fasta_protein_db\" " << endl;
  CommandLineParser cmd(intro.str());
  // define available options
  cmd.defineOption("v",
                   "verbose",
                   "Set verbosity of output: 0 = no processing info, 5 = all, default is 2.",
                   "level");
  cmd.defineOption("i",
                   "peptide_in",
                   "Specifies the file with the peptide tab-delimited output file from percolator.",
                   "filename");
  cmd.defineOption("d",
                   "database",
                   "Specifies the file with protein sequences in fasta format.",
                   "filename");
  cmd.defineOption("o",
                   "protein_out",
                   "Specifies the file with the inferred proteins.",
                   "filename");

  cmd.parseArgs(argc, argv);

  // process options
  if (cmd.optionSet("v")) {
    Globals::getInstance()->setVerbose(cmd.getInt("v", 0, 10));
  }
  if (cmd.optionSet("d")) {
    protein_db_file_ = cmd.options["d"];
  }
  if (cmd.optionSet("i")) {
    peptide_input_file_ = cmd.options["i"];
  }
  if (cmd.optionSet("o")) {
    protein_output_file_ = cmd.options["o"];
  }
  
  return true;
}

bool FisherCaller::getPeptideProteinMap(
    std::map<std::string, std::vector<size_t> >& peptide_protein_map,
    bool generateDecoys) {
  // set a new protein iterator
  DatabaseProteinIterator* database_protein_iterator = new DatabaseProteinIterator(database_);
  if(database_protein_iterator == NULL){
    std::cerr << "Could not create protein iterator." << std::endl;
    return false;
  }
  
  Crux::Protein* protein = NULL;
  
  size_t protein_idx = 0;
  // check if there are any proteins to create peptides from
  while(database_protein_iterator->hasNext()){

    protein = database_protein_iterator->next();
    
    if (generateDecoys) {
      protein->shuffle(PROTEIN_REVERSE_DECOYS);
    }
    //std::cerr << protein->getSequencePointer() << std::endl;
    
    // set new protein peptide iterator
    ProteinPeptideIterator* cur_protein_peptide_iterator =
      new ProteinPeptideIterator(protein, peptide_constraint_);
    
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
  return true;
}

bool FisherCaller::getFragmentProteinMap(
    std::map<std::string, std::vector<size_t> >& peptide_protein_map,
    std::map<size_t, std::vector<size_t> >& fragment_protein_map,
    std::map<size_t, size_t>& num_peptides_per_protein) {
  Crux::Protein* protein = NULL;
  
  for (size_t i = 0; i < database_->getNumProteins(); ++i) {
    protein = database_->getProteinAtIdx(i);
    
    // set new protein peptide iterator
    ProteinPeptideIterator* cur_protein_peptide_iterator =
      new ProteinPeptideIterator(protein, peptide_constraint_);
    
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
}

bool FisherCaller::getProteinFragmentsAndDuplicates(
    std::map<size_t, std::vector<size_t> >& fragment_protein_map,
    std::map<size_t, size_t>& num_peptides_per_protein,
    std::map<std::string, std::string>& fragment_map, 
    std::map<std::string, std::string>& duplicate_map) {
  std::map<size_t, std::vector<size_t> >::iterator it;
  for (it = fragment_protein_map.begin(); it != fragment_protein_map.end(); ++it) {
    size_t i = it->first;
    
    char* c_this_id = database_->getProteinAtIdx(i)->getId();
    std::string this_protein_id(c_this_id);
    std::free(c_this_id);
    
    for (std::vector<size_t>::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
      size_t j = *it2;
      
      if (i != j) {
        char* c_that_id = database_->getProteinAtIdx(j)->getId();
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
}

bool FisherCaller::getProteinFragmentsAndDuplicates(
    std::map<std::string, std::string>& fragment_map,
    std::map<std::string, std::string>& duplicate_map,
    bool generateDecoys) {
  // now create a database_
  database_ = new Database(protein_db_file_.c_str(), false);
  
  if(!database_->parse()){
    std::cerr << "Failed to parse database_, cannot create new index for " << protein_db_file_ << std::endl;
    return EXIT_FAILURE;
  } /*else {
    std::cerr << "Read " << database_->getNumProteins() << " proteins." << std::endl;
    std::cerr << database_->getProteinAtIdx(5)->getSequencePointer() << std::endl;
  }*/
  
  std::map<std::string, std::vector<size_t> > peptide_protein_map;
  bool success = getPeptideProteinMap(peptide_protein_map, generateDecoys);
  
  if (!success) {
    std::cerr << "Failed to create peptide protein map." << std::endl;
    return EXIT_FAILURE;
  } /*else {
    std::cerr << peptide_protein_map["VRPLAR"].size() << std::endl;
    std::cerr << peptide_protein_map["EHHEHASAPLLPPPPTSALSSIASTTAASSAHAK"].size() << std::endl;
    std::cerr << database_->getProteinAtIdx(peptide_protein_map["EHHEHASAPLLPPPPTSALSSIASTTAASSAHAK"][0])->getSequencePointer() << std::endl;
  }*/
  
  std::map<size_t, std::vector<size_t> > fragment_protein_map;
  std::map<size_t, size_t> num_peptides_per_protein;
  
  success = getFragmentProteinMap(peptide_protein_map, fragment_protein_map, num_peptides_per_protein);
  
  if (!success) {
    std::cerr << "Failed to create fragment protein map." << std::endl;
    return EXIT_FAILURE;
  }
  
  success = getProteinFragmentsAndDuplicates(fragment_protein_map, num_peptides_per_protein, fragment_map, duplicate_map);
  
  if (!success) {
    std::cerr << "Failed to get protein fragments and duplicates." << std::endl;
    return EXIT_FAILURE;
  }
  
  return EXIT_SUCCESS;
}
