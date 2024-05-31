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
#include "Version.h"
#include "Option.h"
#include "Globals.h"

using namespace std;
using namespace PercolatorCrux;

PickedProteinCaller::PickedProteinCaller() : enzyme_(TRYPSIN), digestion_(FULL_DIGEST),
    min_peptide_length_(6), max_peptide_length_(50), max_miscleavages_(0),
    decoyPattern_("decoy_"), fasta_has_decoys_(false) {}

PickedProteinCaller::~PickedProteinCaller() {}

void PickedProteinCaller::initConstraints(ENZYME_T enzyme, DIGEST_T digestion, 
    int min_peptide_length, int max_peptide_length, int max_miscleavages) {
  enzyme_ = enzyme;
  digestion_ = digestion;
  min_peptide_length_ = min_peptide_length;
  max_peptide_length_ = max_peptide_length;
  max_miscleavages_ = max_miscleavages;
}

/* introductory message */
string PickedProteinCaller::greeter() const {
  ostringstream oss;
  oss << "Picked-Protein version " << VERSION << ", ";
  oss << "Build Date " << __DATE__ << " " << __TIME__ << endl;
  oss << "Distributed under Apache License" << endl;
  oss << "Written by Lukas Kall (lukas.kall@scilifelab.se) "
         "and Matthew The (matthew.the@scilifelab.se)" << endl;
  oss << "Usage:" << endl;
  oss << "   picked-protein [options]" << endl << endl;
  return oss.str();
}

/* parse the command line arguments */
bool PickedProteinCaller::parseOptions(int argc, char** argv) {
  ostringstream intro;
  intro << greeter() << endl << "Usage:" << endl;
  intro << "   picked-protein -i \"percolator_peptide_output\" -d \"fasta_protein_db\" " << endl;
  CommandLineParser cmd(intro.str());
  // define available options
  cmd.defineOption("v",
                   "verbose",
                   "Set verbosity of output: 0 = no processing info, 5 = all, default is 2.",
                   "level");
  cmd.defineOption("i",
                   "peptide-in",
                   "Specifies the file with the peptide tab-delimited output file from percolator.",
                   "filename");
  cmd.defineOption("d",
                   "database",
                   "Specifies the file with protein sequences in fasta format.",
                   "filename");
  cmd.defineOption("o",
                   "protein-out",
                   "Specifies the file with the inferred proteins.",
                   "filename");

  cmd.parseArgs(argc, argv);

  // process options
  if (cmd.isOptionSet("verbose")) {
    Globals::getInstance()->setVerbose(cmd.getInt("verbose", 0, 10));
  }
  if (cmd.isOptionSet("database")) {
    protein_db_file_ = cmd.options["database"];
  }
  if (cmd.isOptionSet("peptide-in")) {
    peptide_input_file_ = cmd.options["peptide-in"];
  }
  if (cmd.isOptionSet("protein-out")) {
    protein_output_file_ = cmd.options["protein-out"];
  }
  
  return true;
}

void PickedProteinCaller::addToPeptideProteinMap(Database& db, 
    size_t protein_idx, PeptideConstraint& peptide_constraint,
    std::map<std::string, std::vector<size_t> >& peptide_protein_map,
    std::map<size_t, size_t>& num_peptides_per_protein,
    bool reverseProteinSeqs) {
  PercolatorCrux::Protein* protein = db.getProteinAtIdx(static_cast<unsigned int>(protein_idx));
  
  if (reverseProteinSeqs) {
    protein->shuffle(PROTEIN_REVERSE_DECOYS);
    
    // MT: the crux interface will change the protein identifier. If we are
    // not inside the crux environment we do this separately here.
    std::string currentId(protein->getIdPointer());
    if (currentId.substr(0, decoyPattern_.size()) != decoyPattern_) {
      currentId = decoyPattern_ + currentId;
      protein->setId(currentId.c_str());
    }
  } else if (!fasta_has_decoys_) {
    std::string currentId(protein->getIdPointer());
    if (currentId.substr(0, decoyPattern_.size()) == decoyPattern_) {
      fasta_has_decoys_ = true;
    }
  }
  
  ProteinPeptideIterator cur_protein_peptide_iterator(protein, &peptide_constraint);
  size_t numPeptides = 0;
  while (cur_protein_peptide_iterator.hasNext()) {
    PercolatorCrux::Peptide* peptide = cur_protein_peptide_iterator.next();
    std::string sequence(peptide->getSequencePointer(), peptide->getLength());
    peptide_protein_map[sequence].push_back(protein_idx);
    if (sequence[0] == 'M' && peptide->getNTermFlankingAA() == '-'
          && peptide->getLength() - 1 >= min_peptide_length_) {
      std::string metCleavedSequence(sequence.substr(1, static_cast<std::size_t>(peptide->getLength() - 1)));
      peptide_protein_map[metCleavedSequence].push_back(protein_idx);
    }
    PercolatorCrux::Peptide::free(peptide);
    ++numPeptides;
  }
  num_peptides_per_protein[protein_idx] = numPeptides;
}

//!
//! fills \p fragment_protein_map which contains groups of proteins with same or
//! subset peptides 
//! 
//! @param[in] db fasta database representation
//! @param[in] protein_idx protein index of protein to be digested
//! @param[in] peptide_constraint digestion parameters
//! @param[in] num_peptides_per_protein helps to find out which proteins have 
//!   identical sets and which ones are proper subsets.
//! @param[out] fragment_protein_map groups of proteins with same or subset peptides
//!
void PickedProteinCaller::findFragmentProteins(Database& db, 
    size_t protein_idx, PeptideConstraint& peptide_constraint,
    std::map<std::string, std::vector<size_t> >& peptide_protein_map,
    std::map<size_t, size_t>& num_peptides_per_protein,
    std::map<size_t, std::vector<size_t> >& fragment_protein_map) {
  PercolatorCrux::Protein* protein = db.getProteinAtIdx(static_cast<unsigned int>(protein_idx));
  
  ProteinPeptideIterator cur_protein_peptide_iterator(protein, &peptide_constraint);
  
  // find all proteins of which the current protein's peptides are a subset (possibly identical)
  bool is_first = true;
  std::vector<size_t> protein_idx_intersection;
  while (cur_protein_peptide_iterator.hasNext()) {
    PercolatorCrux::Peptide* peptide = cur_protein_peptide_iterator.next();
    std::string sequence(peptide->getSequencePointer(), peptide->getLength());   
    
    if (is_first) {
      protein_idx_intersection = peptide_protein_map[sequence]; // proteins sharing the first peptide
      is_first = false;
    } else {
      std::vector<size_t> newprotein_idx_intersection(protein_idx_intersection.size());
      std::vector<size_t>::iterator it = std::set_intersection(
          protein_idx_intersection.begin(), protein_idx_intersection.end(), 
          peptide_protein_map[sequence].begin(), peptide_protein_map[sequence].end(), 
          newprotein_idx_intersection.begin());
      newprotein_idx_intersection.resize(static_cast<std::size_t>(it - newprotein_idx_intersection.begin()));
      protein_idx_intersection = newprotein_idx_intersection;
    }
    
    if (sequence[0] == 'M' && peptide->getNTermFlankingAA() == '-'
          && peptide->getLength() - 1 >= min_peptide_length_) {
      std::string metCleavedSequence(sequence.substr(1, static_cast<std::size_t>(peptide->getLength() - 1)));
      std::vector<size_t> newprotein_idx_intersection(protein_idx_intersection.size());
      std::vector<size_t>::iterator it = std::set_intersection(
          protein_idx_intersection.begin(), protein_idx_intersection.end(), 
          peptide_protein_map[metCleavedSequence].begin(), peptide_protein_map[metCleavedSequence].end(), 
          newprotein_idx_intersection.begin());
      newprotein_idx_intersection.resize(static_cast<std::size_t>(it - newprotein_idx_intersection.begin()));
      protein_idx_intersection = newprotein_idx_intersection;
    }
     
    PercolatorCrux::Peptide::free(peptide);
    if (protein_idx_intersection.size() < 2) break;
  }
  
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

void PickedProteinCaller::findFragmentsAndDuplicates(Database& db,
    std::map<size_t, std::vector<size_t> >& fragment_protein_map,
    std::map<size_t, size_t>& num_peptides_per_protein,
    std::map<std::string, std::string>& fragment_map, 
    std::map<std::string, std::string>& duplicate_map) {
  std::map<size_t, std::vector<size_t> >::iterator it;
  for (it = fragment_protein_map.begin(); it != fragment_protein_map.end(); ++it) {
    size_t i = it->first;
    std::string this_protein_id(db.getProteinAtIdx(static_cast<unsigned int>(i))->getIdPointer());
    
    for (std::vector<size_t>::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
      size_t j = *it2;
      
      if (i != j) {
        std::string that_protein_id(db.getProteinAtIdx(static_cast<unsigned int>(j))->getIdPointer());
        
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

//! creates a local peptide->protein map for each candidate protein group that 
//! which was based on a full digest with max_len <= 50 and max_misclv <= 2
void PickedProteinCaller::findFragmentsAndDuplicatesExtraDigest(Database& db, 
    PeptideConstraint& peptide_constraint,
    std::map<size_t, std::vector<size_t> >& fragment_protein_map,
    std::map<std::string, std::string>& fragment_map, 
    std::map<std::string, std::string>& duplicate_map) {
  bool reverseProteinSeqs = false;
  std::map<size_t, std::vector<size_t> >::iterator it;
  for (it = fragment_protein_map.begin(); it != fragment_protein_map.end(); ++it) {
    size_t i = it->first;
    if (std::find(it->second.begin(), it->second.end(), i) == it->second.end()) {
      it->second.push_back(i);
    }
    std::sort(it->second.begin(), it->second.end());
    
    std::map<std::string, std::vector<size_t> > peptide_protein_map;
    std::map<size_t, size_t> num_peptides_per_protein_local;
    for (std::vector<size_t>::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
      addToPeptideProteinMap(db, *it2, peptide_constraint, 
          peptide_protein_map, num_peptides_per_protein_local, reverseProteinSeqs);
    }
    
    std::map<size_t, std::vector<size_t> > fragment_protein_map_local;
    for (std::vector<size_t>::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
      findFragmentProteins(db, *it2, peptide_constraint, 
          peptide_protein_map, num_peptides_per_protein_local, 
          fragment_protein_map_local);
    }
    
    findFragmentsAndDuplicates(db, fragment_protein_map_local, 
        num_peptides_per_protein_local, fragment_map, duplicate_map);
  }
}

void PickedProteinCaller::findFragmentsAndDuplicatesNonSpecificDigest(
    Database& db, 
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
    
    std::vector<std::string> sequences;
    std::map<size_t, size_t> num_peptides_per_protein_local;
    for (std::vector<size_t>::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
      size_t protein_idx = *it2;
      PercolatorCrux::Protein* protein = db.getProteinAtIdx(static_cast<unsigned int>(protein_idx));
      std::string sequence(protein->getSequencePointer(), protein->getLength());
      sequences.push_back(sequence);
      num_peptides_per_protein_local[protein_idx] = protein->getLength();
    }
    
    // In a non-specific digest the only possibility for one protein to be subset
    // of another is if the entire string is contained (except for some very
    // unlikely cases where a region longer than max_len is repeated more than twice)
    std::map<size_t, std::vector<size_t> > fragment_protein_map_local;
    std::vector<std::string>::iterator sit2 = sequences.begin();
    for (std::vector<size_t>::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2, ++sit2) {
      size_t protein_idx = *it2;
      
      std::vector<size_t> protein_idx_intersection;
      std::vector<size_t>::iterator it3 = it->second.begin();
      for (std::vector<std::string>::iterator sit3 = sequences.begin(); sit3 != sequences.end(); ++sit3, ++it3) {
        if (sit3->find(*sit2) != std::string::npos) {
          protein_idx_intersection.push_back(*it3);
        }
      }
      
      if (protein_idx_intersection.size() > 1) {
        addToFragmentProteinMap(protein_idx, protein_idx_intersection, 
            num_peptides_per_protein_local, fragment_protein_map_local);
      }
    }
    
    findFragmentsAndDuplicates(db, fragment_protein_map_local, 
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
  
  bool is_memmap = false;
  Database db(protein_db_file_.c_str(), is_memmap);
  
  if (!db.parse()) {
    std::cerr << "Failed to parse database, cannot create index for " 
              << protein_db_file_ << std::endl;
    return EXIT_FAILURE;
  }
  
  if (VERB > 3) reportProgress("Creating database", startTime, startClock);
  
  // First do a "basic" digest to get candidates for protein grouping*
  //
  // * there are some rare cases in which fragment proteins have
  // distinct subsets in a full digest, but form a proper subset in a
  // partial or non-enzymatic digest (e.g. if one would add a single
  // amino acid to the beginning of a protein sequence). It's too hard 
  // to check for these at the moment...
  //
  PeptideConstraint peptide_constraint(enzyme_, FULL_DIGEST, 
      min_peptide_length_, (std::min)(50, max_peptide_length_), 
      (std::min)(2, max_miscleavages_) );
  std::map<std::string, std::vector<size_t> > peptide_protein_map;
  std::map<size_t, size_t> num_peptides_per_protein;
  for (size_t protein_idx = 0; protein_idx < db.getNumProteins(); 
       ++protein_idx) {
    addToPeptideProteinMap(db, protein_idx, peptide_constraint, 
        peptide_protein_map, num_peptides_per_protein, reverseProteinSeqs);
  }
  
  if (VERB > 3) {
    reportProgress("Creating protein peptide map", startTime, startClock);
  }
  
  // Find all proteins whose peptides form a subset (possibly identical) 
  // of another protein
  std::map<size_t, std::vector<size_t> > fragment_protein_map;
  for (size_t protein_idx = 0; protein_idx < db.getNumProteins(); 
       ++protein_idx) {
    findFragmentProteins(db, protein_idx, peptide_constraint,
        peptide_protein_map, num_peptides_per_protein, fragment_protein_map);
  }
  
  if (VERB > 3) {
    reportProgress("Creating fragment protein map", startTime, startClock);
  }
  
  // If the actual digest was more permissive than the basic digest (see above), 
  // check validity of each candidate protein group
  if (digestion_ == FULL_DIGEST && max_miscleavages_ <= 2 
          && max_peptide_length_ <= 50) {
    findFragmentsAndDuplicates(db, fragment_protein_map, 
                  num_peptides_per_protein, fragment_map, duplicate_map);
  } else if (digestion_ != NON_SPECIFIC_DIGEST) {
    PeptideConstraint peptide_constraint_extra_digest(
        enzyme_, digestion_, min_peptide_length_, max_peptide_length_, 
        max_miscleavages_);
    findFragmentsAndDuplicatesExtraDigest(db, 
                  peptide_constraint_extra_digest, fragment_protein_map, 
                  fragment_map, duplicate_map);
  } else {
    findFragmentsAndDuplicatesNonSpecificDigest(db, 
                  fragment_protein_map, fragment_map, duplicate_map);
  }
  
  if (VERB > 2) {
    reportProgress("Resolving protein fragments and duplicates", 
                      startTime, startClock);
  }
  
  return EXIT_SUCCESS;
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
