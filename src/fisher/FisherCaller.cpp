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

FisherCaller::FisherCaller() : enzyme_(TRYPSIN), digestion_(FULL_DIGEST),
    min_peptide_length_(6), max_peptide_length_(50), max_miscleavages_(0),
    decoyPattern_("decoy_") {}

FisherCaller::~FisherCaller() {}

void FisherCaller::initConstraints(ENZYME_T enzyme, DIGEST_T digestion, 
    int min_peptide_length, int max_peptide_length, int max_miscleavages) {
  enzyme_ = enzyme;
  digestion_ = digestion;
  min_peptide_length_ = min_peptide_length;
  max_peptide_length_ = max_peptide_length;
  max_miscleavages_ = max_miscleavages;
}

/* introductory message */
string FisherCaller::greeter() const {
  ostringstream oss;
  oss << "Fisher version " << VERSION << ", ";
  oss << "Build Date " << __DATE__ << " " << __TIME__ << endl;
  oss << "Distributed under Apache License" << endl;
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

bool FisherCaller::getPeptideProteinMap(Database& db, 
    PeptideConstraint& peptide_constraint,
    std::map<std::string, std::vector<size_t> >& peptide_protein_map,
    bool generateDecoys) {
  std::vector< std::pair<std::string, size_t> > peptideProteinIdxMap;
  for (size_t protein_idx = 0; protein_idx < db.getNumProteins(); 
       ++protein_idx) {
    addProteinToPeptideProteinMap(db, protein_idx, peptide_constraint, 
        peptide_protein_map, generateDecoys);
    /*Crux::Protein* protein = db.getProteinAtIdx(protein_idx);
    
    if (generateDecoys) {
      protein->shuffle(PROTEIN_REVERSE_DECOYS);
      
      // MT: the crux interface will change the protein identifier. If we are
      // not inside the crux environment we do this separately here.
      std::string currentId(protein->getIdPointer());
      if (currentId.substr(0, decoyPattern_.size()) != decoyPattern_) {
        currentId = decoyPattern_ + currentId;
        protein->setId(currentId.c_str());
      }
    }
    
    ProteinPeptideIterator cur_protein_peptide_iterator(protein, &peptide_constraint);
    while (cur_protein_peptide_iterator.hasNext()) {
      Crux::Peptide* peptide = cur_protein_peptide_iterator.next();
      std::string sequence(peptide->getSequencePointer(), peptide->getLength());
      peptideProteinIdxMap.push_back(std::make_pair(sequence,protein_idx));
      delete peptide;
    }*/
  }
  /*
  std::string prevPeptide = "";
  size_t count = 0;
  std::vector<size_t> proteinIdxs;
  std::sort(peptideProteinIdxMap.begin(), peptideProteinIdxMap.end());
  for (size_t i = 0; i < peptideProteinIdxMap.size(); ++i) {
    if (peptideProteinIdxMap[i].first == prevPeptide) {
      ++count;
      proteinIdxs.push_back(peptideProteinIdxMap[i-1].second);
    } else {
      if (count > 0) {
        proteinIdxs.push_back(peptideProteinIdxMap[i-1].second);
        peptide_protein_map[prevPeptide] = proteinIdxs;
        proteinIdxs.clear();
        count = 0;
      }
      prevPeptide = peptideProteinIdxMap[i].first;
    }
  }
  
  if (count > 0) {
    proteinIdxs.push_back(peptideProteinIdxMap[i-1].second);
    peptide_protein_map[prevPeptide] = proteinIdxs;
  }
  */
  return true;
}

void FisherCaller::addProteinToPeptideProteinMap(Database& db, 
    size_t protein_idx, PeptideConstraint& peptide_constraint,
    std::map<std::string, std::vector<size_t> >& peptide_protein_map,
    bool generateDecoys) {
  Crux::Protein* protein = db.getProteinAtIdx(protein_idx);
    
  if (generateDecoys) {
    protein->shuffle(PROTEIN_REVERSE_DECOYS);
    
    // MT: the crux interface will change the protein identifier. If we are
    // not inside the crux environment we do this separately here.
    std::string currentId(protein->getIdPointer());
    if (currentId.substr(0, decoyPattern_.size()) != decoyPattern_) {
      currentId = decoyPattern_ + currentId;
      protein->setId(currentId.c_str());
    }
  }
  
  ProteinPeptideIterator cur_protein_peptide_iterator(protein, &peptide_constraint);
  while (cur_protein_peptide_iterator.hasNext()) {
    Crux::Peptide* peptide = cur_protein_peptide_iterator.next();
    std::string sequence(peptide->getSequencePointer(), peptide->getLength());
    peptide_protein_map[sequence].push_back(protein_idx);
    delete peptide;
  }
}

bool FisherCaller::getFragmentProteinMap(Database& db, 
    PeptideConstraint& peptide_constraint,
    std::map<std::string, std::vector<size_t> >& peptide_protein_map,
    std::map<size_t, std::vector<size_t> >& fragment_protein_map,
    std::map<size_t, size_t>& num_peptides_per_protein) {
  for (size_t protein_idx = 0; protein_idx < db.getNumProteins(); 
       ++protein_idx) {
    addProteinToFragmentProteinMap(db, protein_idx, peptide_constraint,
        peptide_protein_map, fragment_protein_map, num_peptides_per_protein);
  }
  return true;
}

void FisherCaller::addProteinToFragmentProteinMap(Database& db, 
    size_t protein_idx, PeptideConstraint& peptide_constraint,
    std::map<std::string, std::vector<size_t> >& peptide_protein_map,
    std::map<size_t, std::vector<size_t> >& fragment_protein_map,
    std::map<size_t, size_t>& num_peptides_per_protein) {
  Crux::Protein* protein = db.getProteinAtIdx(protein_idx);
    
  // set new protein peptide iterator
  ProteinPeptideIterator cur_protein_peptide_iterator(protein, &peptide_constraint);
  
  bool is_first = true;
  size_t num_sequences = 0;
  std::vector<size_t> protein_idx_intersection;
  while (cur_protein_peptide_iterator.hasNext()) {
    Crux::Peptide* peptide = cur_protein_peptide_iterator.next();
    std::string sequence(peptide->getSequencePointer(), peptide->getLength());    
    delete peptide;
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
    
    if (protein_idx_intersection.size() < 2) break;
  }
  
  if (protein_idx_intersection.size() > 1) {
    num_peptides_per_protein[protein_idx] = num_sequences;
    std::sort(protein_idx_intersection.begin(), protein_idx_intersection.end()); // sort in descending order
    
    size_t j = protein_idx_intersection.size() - 1;
    if (protein_idx_intersection[j] == protein_idx) --j;
    
    if (fragment_protein_map.find(protein_idx) != fragment_protein_map.end()) {
      fragment_protein_map[protein_idx_intersection[j]].insert(
          fragment_protein_map[protein_idx_intersection[j]].end(),
          fragment_protein_map[protein_idx].begin(), 
          fragment_protein_map[protein_idx].end());
      fragment_protein_map.erase(protein_idx);
    }
    fragment_protein_map[protein_idx_intersection[j]].push_back(protein_idx);
  }
}

bool FisherCaller::getProteinFragmentsAndDuplicates(Database& db,
    std::map<size_t, std::vector<size_t> >& fragment_protein_map,
    std::map<size_t, size_t>& num_peptides_per_protein,
    std::map<std::string, std::string>& fragment_map, 
    std::map<std::string, std::string>& duplicate_map) {
  std::map<size_t, std::vector<size_t> >::iterator it;
  for (it = fragment_protein_map.begin(); it != fragment_protein_map.end(); ++it) {
    size_t i = it->first;
    std::string this_protein_id(db.getProteinAtIdx(i)->getIdPointer());
    
    for (std::vector<size_t>::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
      size_t j = *it2;
      
      if (i != j) {
        std::string that_protein_id(db.getProteinAtIdx(j)->getIdPointer());
        
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
  return true;
}

bool FisherCaller::getProteinFragmentsAndDuplicatesExtraDigest(Database& db, 
    PeptideConstraint& peptide_constraint,
    std::map<size_t, std::vector<size_t> >& fragment_protein_map,
    std::map<size_t, size_t>& num_peptides_per_protein,
    std::map<std::string, std::string>& fragment_map, 
    std::map<std::string, std::string>& duplicate_map) {
  bool generateDecoys = false;
  std::map<size_t, std::vector<size_t> >::iterator it;
  for (it = fragment_protein_map.begin(); it != fragment_protein_map.end(); ++it) {
    size_t i = it->first;
    std::map<std::string, std::vector<size_t> > peptide_protein_map;
    addProteinToPeptideProteinMap(db, i, peptide_constraint, 
        peptide_protein_map, generateDecoys);
    for (std::vector<size_t>::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
      size_t j = *it2;
      if (i != j) {
        addProteinToPeptideProteinMap(db, j, peptide_constraint, 
            peptide_protein_map, generateDecoys);
      }
    }
    
    std::map<size_t, std::vector<size_t> > fragment_protein_map_local;
    std::map<size_t, size_t> num_peptides_per_protein_local;
    addProteinToFragmentProteinMap(db, i, peptide_constraint, peptide_protein_map, 
        fragment_protein_map_local, num_peptides_per_protein_local);
    for (std::vector<size_t>::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
      size_t j = *it2;
      if (i != j) {
        addProteinToFragmentProteinMap(db, j, peptide_constraint, 
            peptide_protein_map, fragment_protein_map_local, 
            num_peptides_per_protein_local);
      }
    }
    
    getProteinFragmentsAndDuplicates(db, fragment_protein_map_local, 
        num_peptides_per_protein_local, fragment_map, duplicate_map);
  }
  return true;
}

bool FisherCaller::getProteinFragmentsAndDuplicates(
    std::map<std::string, std::string>& fragment_map,
    std::map<std::string, std::string>& duplicate_map,
    bool generateDecoys) {
  time_t startTime;
  time(&startTime);
  clock_t startClock = clock();
  
  // now create a database
  Database db(protein_db_file_.c_str(), false);
  
  if(!db.parse()){
    std::cerr << "Failed to parse database, cannot create new index for " << protein_db_file_ << std::endl;
    return EXIT_FAILURE;
  }
  
  time_t procStart;
  clock_t procStartClock = clock();
  time(&procStart);
  double diff = difftime(procStart, startTime);
  if (VERB > 3) cerr << "Creating database took "
      << ((double)(procStartClock - startClock)) / (double)CLOCKS_PER_SEC
      << " cpu seconds or " << diff << " seconds wall time" << endl;
  
  // First do a full digest to get candidates for protein grouping
  PeptideConstraint peptide_constraint(enzyme_, FULL_DIGEST, 
      min_peptide_length_, max_peptide_length_, max_miscleavages_);
  std::map<std::string, std::vector<size_t> > peptide_protein_map;
  bool success = getPeptideProteinMap(db, peptide_constraint, peptide_protein_map, generateDecoys);
  
  if (!success) {
    std::cerr << "Failed to create peptide protein map." << std::endl;
    return EXIT_FAILURE;
  }
  
  procStartClock = clock();
  time(&procStart);
  diff = difftime(procStart, startTime);
  if (VERB > 3) cerr << "Creating protein peptide map took "
      << ((double)(procStartClock - startClock)) / (double)CLOCKS_PER_SEC
      << " cpu seconds or " << diff << " seconds wall time" << endl;
      
  std::map<size_t, std::vector<size_t> > fragment_protein_map;
  std::map<size_t, size_t> num_peptides_per_protein;
  
  success = getFragmentProteinMap(db, peptide_constraint, peptide_protein_map, fragment_protein_map, num_peptides_per_protein);
  
  if (!success) {
    std::cerr << "Failed to create fragment protein map." << std::endl;
    return EXIT_FAILURE;
  }
  
  procStartClock = clock();
  time(&procStart);
  diff = difftime(procStart, startTime);
  if (VERB > 3) cerr << "Creating fragment protein map took "
      << ((double)(procStartClock - startClock)) / (double)CLOCKS_PER_SEC
      << " cpu seconds or " << diff << " seconds wall time" << endl;
  
  // If not a full digest, check validity of each candidate protein group
  if (digestion_ == FULL_DIGEST) {
    success = getProteinFragmentsAndDuplicates(db, fragment_protein_map, num_peptides_per_protein, fragment_map, duplicate_map);
  } else {
    PeptideConstraint peptide_constraint_extra_digest(enzyme_, digestion_, 
        min_peptide_length_, max_peptide_length_, max_miscleavages_);
    success = getProteinFragmentsAndDuplicatesExtraDigest(db, 
        peptide_constraint_extra_digest, fragment_protein_map, 
        num_peptides_per_protein, fragment_map, duplicate_map);
  }
  
  procStartClock = clock();
  time(&procStart);
  diff = difftime(procStart, startTime);
  if (VERB > 2) cerr << "Resolving protein fragments and duplicates map took "
      << ((double)(procStartClock - startClock)) / (double)CLOCKS_PER_SEC
      << " cpu seconds or " << diff << " seconds wall time" << endl;
  
  
  if (!success) {
    std::cerr << "Failed to get protein fragments and duplicates." << std::endl;
    return EXIT_FAILURE;
  }
  
  return EXIT_SUCCESS;
}
