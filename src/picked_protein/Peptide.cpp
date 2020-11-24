/*******************************************
 *** adapted from crux/src/c/Peptide.cpp ***
 *******************************************/

/*************************************************************************//**
 * \file Peptide.cpp
 * \brief Object for representing a single peptide.
 ****************************************************************************/
#include "Peptide.h"
#include "PeptideSrc.h"
#include <string.h>

#include <set>
#include <vector>
#include <algorithm>
#include <sstream>

using namespace std;
using namespace PercolatorCrux;
/*
  TABLE OF CONTENTS
  Global variables
  Private data types
  Private functions
  Public functions
    Allocators/deallocators
    Getters and Setters
      mass-related
      source-related
      sequence-related
      getters requiring calculation

    Comparisons for sorting // MOVEME
    Printing / parsing      // MOVEME
      text
      binary

    Iterators
      residue iterator
      source iterator

      
 */


// the struct to be printed and read; possible fix for adding a field to the peptide
struct PRINT_PEPTIDE_T {
  unsigned char length; ///< The length of the peptide
  vector<PeptideSrc*>* peptide_src; ///< a vector of peptide_src  
};

/**
 * \struct peptide_src_iterator
 * \brief Object to iterate over the peptide_srcs Vector in a peptide
 */
struct peptide_src_iterator{
  Peptide*  peptide; ///< The peptide whose peptide_srcs to iterate over.
  PeptideSrc* current; ///< the current peptide_srcs
};

/* Private functions */

/**
 * Initializes an (empty) peptide object
 */
void Peptide::init() {
  length_ = 0;
  peptide_srcs_.clear();
}

/* Public functions--Allocators/Deallocators */



/**
 * \returns An (empty) peptide object.
 */
Peptide::Peptide() {
  init();
}

Peptide* Peptide::copyPtr() {

  pointer_count_++;
  return this;
}

void Peptide::free(Peptide* peptide) {

  peptide->pointer_count_--;
  if (peptide->pointer_count_ <= 0) {
    delete peptide;
  }
}


/**
 *\returns the protein struct size, value of sizeof function
 */
int Peptide::getSizeOf(){
  return sizeof(Peptide);
}

// FIXME association part might be need to change
/**
 * \returns A new peptide object, populated with the user specified parameters.
 */
  Peptide::Peptide(
  unsigned char length,     ///< The length of the peptide -in
  Protein* parent_protein, ///< the parent_protein of this peptide -in
  int start_idx ///< the start index of this peptide in the protein sequence -in
  ) {

  init();
  setLength(length);
  // FIXME: find the level of digest for this specific protein
  peptide_srcs_.push_back(new PeptideSrc(NON_SPECIFIC_DIGEST, parent_protein, start_idx ));
}

  
/**
 * \brief Allocates a new peptide giving it the values of the source
 * peptide.
 * \returns A newly allocated peptide identical to the source.
 */
Peptide::Peptide(
  Peptide* src ///< source peptide -in
){
  
  init();

  if( src == NULL ){
    //carp(CARP_ERROR, "Cannot copy null peptide!");
  } else {
    length_ = src->length_;
    PeptideSrc::copy(src->peptide_srcs_,peptide_srcs_);
  }
}

/**
 * Merge two identical peptides, copy all peptide_src into one of the peptide
 * peptide_dest, peptide_bye must have at least one peptide src
 * frees the peptide_bye, once the peptide_src are re-linked to the peptide_dest
 * \returns true if merge is successful else false
 */
bool Peptide::mergePeptides(
  Peptide* peptide_dest, ///< the peptide to merge into  -out
  Peptide* peptide_bye ///< the peptide to be merged  -in
  ) {

  vector<PeptideSrc*>& dest_srcs = peptide_dest->peptide_srcs_;
  vector<PeptideSrc*>& bye_srcs = peptide_bye->peptide_srcs_;
  // do both peptides have at least one peptide_src?
  if(dest_srcs.empty() || bye_srcs.empty()){
    //carp(CARP_ERROR, "failed to merge two peptides");
    return false;
  }

  for(vector<PeptideSrc*>::iterator iter=bye_srcs.begin();
    iter!=bye_srcs.end();
    iter++
  ){
   
    dest_srcs.push_back(*iter);
  }

  // peptide_dest now points to the src's, don't delete them
  bye_srcs.clear();
  
  delete peptide_bye;

  return true;
}

/**
 * Merges two identical peptides by adding the peptide_src of the
 * second to the first.  The second peptide remains unchanged.
 * Does not comfirm identity of peptides.
 * \returns true if merge is successfull.
 */
bool Peptide::mergePeptidesCopySrc(
  Peptide* peptide_dest,
  Peptide* peptide_giver
  ){

  vector<PeptideSrc*>& dest_srcs = peptide_dest->peptide_srcs_;
  vector<PeptideSrc*>& giver_srcs = peptide_giver->peptide_srcs_;

  // do both peptides have at least one peptide_src?
  if(dest_srcs.empty() || giver_srcs.empty()){
    //carp(CARP_ERROR, "failed to merge two peptides");
    return false;
  }
 
  for(vector<PeptideSrc*>::iterator iter=giver_srcs.begin();
    iter != giver_srcs.end();
    iter++
  ){
    PeptideSrc* new_src = new PeptideSrc(*(*iter));
    dest_srcs.push_back(new_src);
  }

  return true;
}

/**
 * Frees an allocated peptide object.
 * Depending on peptide_src implementation determines how to free srcs
 */
Peptide::~Peptide() {
  
   PeptideSrc::free(peptide_srcs_);
}

// Public functions--Getters and Setters 

/* source-related getters and setters */

/**
 * sets the peptide_src field in the peptide
 * this method should be ONLY used when the peptide has no existing list of peptide_src
 * must pass on a heap allocated peptide_srcs_ object
 * does not copy in the object, just the pointer to the object.
 */
void Peptide::setPeptideSrc(
  PeptideSrc*  new_association ///< new peptide_src -in
  ) {

  //assert(peptide_srcs_.empty());
  peptide_srcs_.push_back(new_association);
}

/**
 * this method adds the new_association to the end of the existing peptide's 
 * if no prior existing list, adds it at the front
 * must pass on a heap allocated peptide_srcs_ object
 * does not copy in the object, just the pointer to the object.
 */
void Peptide::addPeptideSrc(
  PeptideSrc* new_association ///< new peptide_src -in
  ) {

  peptide_srcs_.push_back(new_association);

}


// TODO: why do we need both of these?
/**
 * returns a pointer to the peptide_protein_association field of the peptide
 */
PeptideSrc* Peptide::getPeptideSrc() {

  
  return peptide_srcs_.at(0);
}
/**
 *returns the pepide_srcs_
 */

vector<PeptideSrc*>& Peptide::getPeptideSrcVector() {
  return peptide_srcs_;
}

/**
 *Return the begining of the peptide_srcs_
 */
PeptideSrcIterator Peptide::getPeptideSrcBegin() {

  return peptide_srcs_.begin();
}
/**
 *Return the end of the peptide_srcs_
 */

PeptideSrcIterator Peptide::getPeptideSrcEnd() {
  return peptide_srcs_.end();
}


/**
 * \returns The number of peptide sources (i.e. proteins) the peptide has.
 */
int Peptide::getNumPeptideSrc(){

  return peptide_srcs_.size();
}

/**
 * get the peptide->first peptide_src->parent protein->database
 */
Database* Peptide::getFirstSrcDatabase() {

  return peptide_srcs_.at(0)->getParentProtein()->getDatabase();
}

// set by peptide_src?
/**
 * returns a pointer to the peptide's first parent protein field of the peptide
 */
Protein* Peptide::getParentProtein() {

  return peptide_srcs_.at(0)->getParentProtein();
}

/**
 * sets the sequence length of the peptide
 * length maximum of 255
 */
void Peptide::setLength(
  unsigned char length  ///< the length of sequence -in
  ) {

  length_ = length;
}

/* sequence-related getters and setters */
/**
 *\returns the sequence length of the peptide
 */
unsigned char Peptide::getLength() {

  return length_;
}

/**
 * \brief Get a string representation of the peptide sequence with no
 * added modification symbols.
 * Returns decoy sequence for decoy peptides.
 * \returns The newly-allocated sequence of peptide
 */
char* Peptide::getSequence() {

  //Check that the peptide has parent(s)
  if(peptide_srcs_.empty()   ){
    //carp(CARP_ERROR, "Cannot get sequence from peptide with no peptide src.");
    return NULL;
  }

  char* seq_copy = NULL;

  seq_copy = getUnshuffledSequence();
 
  return seq_copy; 
}

/**
 * \brief Get a string representation of the target (unshuffled)
 * peptide sequence with no added modification symbols.
 * For target peptides, returns the same as get_peptide_sequence.
 * \returns The newly-allocated sequence of peptide
 */
char* Peptide::getUnshuffledSequence() {

  //check  parent(s) protein of the peptide not be empty 
  if(peptide_srcs_.empty()){
    //carp(CARP_ERROR, "Cannot get sequence from peptide with no peptide src.");
    return NULL;
  }
  int start_idx = getPeptideSrc()->getStartIdx();
  char* parent_sequence = 
    getPeptideSrc()->getParentProtein()->getSequencePointer(start_idx-1);
  

  char* copy_sequence = copy_string_part(parent_sequence,
                                         length_);
 
  return copy_sequence; 
}

/**
 * \brief Get a pointer to the peptide sequence that is NOT null
 * terminated.
 * USE WITH CAUTION.  Pointer is to the parent protein sequence and
 * thus is not null-terminated until the end of the protein.  Parent
 * protein is taken from the first protein source.
 * 
 * \returns A pointer to an existing peptide sequence.
 */
char* Peptide::getSequencePointer() {

  if(peptide_srcs_.empty()){
    //carp(CARP_FATAL, "ERROR: no peptide_src to retrieve peptide sequence pointer\n");
  }

  int start_idx = getPeptideSrc()->getStartIdx();
  char* parent_sequence = 
    getPeptideSrc()->getParentProtein()->getSequencePointer(start_idx-1);

  char* pointer_peptide_sequence = parent_sequence;
  
  return pointer_peptide_sequence;
}

/**
 * \brief Return a char for the amino acid n-terminal to the peptide
 * in the peptide src at the given index.
 *
 * \returns A char (A-Z) or - if peptide is the first in the protein.
 */
char Peptide::getNTermFlankingAA() {
  PeptideSrc* src = getPeptideSrc();
  if (src == NULL) {
    return '-';
  }

  // get protein seq
  Protein* protein = src->getParentProtein();

  // get peptide start idx, protein index starts at 1
  int start_index = src->getStartIdx();

  char* protein_seq = protein->getSequencePointer();

  char aa = '-';
  // if not at beginning, return char
  if( start_index > 1 ){
    aa = protein_seq[start_index - 2]; // -1 for 1-based shift
                                       // -1 for aa before start
  }
  return aa;
}

/**
 * Examines the peptide sequence and counts how many tryptic missed
 * cleavage sites exist. 
 *\returns the number of missed cleavage sites in the peptide
 */
int Peptide::getMissedCleavageSites() {

  int missed_count = 0;
  int aa_idx = 0;
  char* sequence = getSequencePointer();

  // count the missed cleavage sites
  for(; aa_idx < length_-1; ++aa_idx){
    if(sequence[aa_idx] == 'K' ||
       sequence[aa_idx] == 'R'){
      
      // skip one that are followed by a P
      if(sequence[aa_idx+1] == 'P'){
        continue;
      }
      else{
        ++missed_count;
      }      
    } 
  }
  
  return missed_count;
}

int Peptide::getMissedCleavageSites(
  set<int> skip ///< skip these amino acid indices.
  ) {

  int missed_count = 0;
  char* sequence = getSequencePointer();

  // count the missed cleavage sites
  for(int aa_idx=0; aa_idx < length_-1; ++aa_idx){
    
    if (skip.find(aa_idx) != skip.end()) {
      continue;
    }

    bool cleavage_prevented = true;
    if(sequence[aa_idx] == 'K' ||
       sequence[aa_idx] == 'R'){

      cleavage_prevented = false;
      // skip one that are followed by a P
      if(sequence[aa_idx+1] == 'P'){
        cleavage_prevented = true;
      }
    }

    if (!cleavage_prevented) {
      ++missed_count;
    }
  }
  
  return missed_count;
}

/* Comparisons for sorting */

/**
 * Compare peptide sequence
 * \returns true if peptide sequence is identical else false
 */ 
bool Peptide::compareSequence(
  Peptide* peptide_one,  ///< the peptide sequence to compare  -out
  Peptide* peptide_two  ///< the peptide sequence to compare  -out
  )
{
  // is length identical?
  if(peptide_one->length_ != peptide_two->length_){
    return false;
  }
  else{
    int current_idx = 0;
    char* start_one = peptide_one->getPeptideSrc()->getSequencePointer();
    char* start_two = peptide_two->getPeptideSrc()->getSequencePointer();
    
    while(current_idx < peptide_one->length_){
      if(start_one[current_idx] != start_two[current_idx]){
        return false;
      }
      ++current_idx;
    } 
  }
  return true;
}

/**
 * Compare two peptide sequences.
 * \returns Zero (0) if the sequences are identical, -1 if the first
 * sequence is less than the first and 1 if the first sequence is
 * greater than teh first.
 */
int Peptide::triCompareSequence(
  Peptide* peptide_one,  ///< the peptide sequence to compare  -out
  Peptide* peptide_two  ///< the peptide sequence to compare  -out
  )
{
  // find the shorter peptide
  int short_len = 0;
  if( peptide_one->length_ < peptide_two->length_ ){
    short_len = peptide_one->length_;
  } else {
    short_len = peptide_two->length_;
  }

  char* seq_one = peptide_one->getPeptideSrc()->getSequencePointer();
  char* seq_two = peptide_two->getPeptideSrc()->getSequencePointer();
    
  // stop comparing as soon as they differ
  int pep_idx = 0;
  for(pep_idx = 0; pep_idx < short_len; pep_idx++ ){
      if(seq_one[pep_idx] != seq_two[pep_idx]){
        pep_idx++; // stop pointing one after the difference
        break;
      }
  }

  // move index back one to the last letter compared
  pep_idx--;

  // if the seqs are the same up to this point, then compare the lengths
  if( seq_one[pep_idx] == seq_two[pep_idx] ){
    if( peptide_one->length_ == peptide_two->length_ ){ // same seq
      return 0;
    } else if ( peptide_one->length_ < peptide_two->length_ ){ 
      return -1;
    } else {
      return 1;
    }
  }

  // else, the seqs are different
  if(seq_one[pep_idx] < seq_two[pep_idx]){
    return -1;
  } else {
    return 1;
  }

}
/**
 * Compare the sequence of two peptides and return true if the first
 * petpide sequence is less than (in a lexical sort) the second
 * peptide.  Return false if they are idential peptides.
 */
bool Peptide::lessThan(
  Peptide* peptide_one,
  Peptide* peptide_two
  ){
  // find the shorter peptide
  int short_len = 0;
  if( peptide_one->length_ < peptide_two->length_ ){
    short_len = peptide_one->length_;
  } else {
    short_len = peptide_two->length_;
  }

  char* seq_one = peptide_one->getPeptideSrc()->getSequencePointer();
  char* seq_two = peptide_two->getPeptideSrc()->getSequencePointer();
    
  // stop comparing as soon as they differ
  int pep_idx = 0;
  for(pep_idx = 0; pep_idx < short_len; pep_idx++ ){
      if(seq_one[pep_idx] != seq_two[pep_idx]){
        break;
      }
  }
  if( seq_one[pep_idx] == seq_two[pep_idx] ){
    return (peptide_one->length_ < peptide_two->length_);
  } 

  return (seq_one[pep_idx] < seq_two[pep_idx]);
}

/**
 * compares two peptides with the lexical sort type
 * for qsort
 * /returns 1 if peptide_one has lower priority, 0 if equal, -1 if greater priority
 */
int Peptide::compareLexicalQSort(
  Peptide** peptide_one, ///< peptide to compare one -in
  Peptide** peptide_two ///< peptide to compare two -in
  )
{
  // convert the protein to heavy if needed
  /* uncomment if needed, to use light heavy protein
  protein_to_heavy(get_peptide_parent_protein(peptide_one));
  protein_to_heavy(get_peptide_parent_protein(peptide_two));
  */
  Peptide* peptide_1 = *peptide_one;
  Peptide* peptide_2 = *peptide_two;
  char* peptide_one_sequence = peptide_1->getSequencePointer();
  char* peptide_two_sequence = peptide_2->getSequencePointer();
  int peptide_one_length = peptide_1->length_;
  int peptide_two_length = peptide_2->length_;
  int current_idx = 0;
  
  // check if all alphabetically identical
  while(current_idx < peptide_one_length &&
        current_idx < peptide_two_length){
    if(peptide_one_sequence[current_idx] > peptide_two_sequence[current_idx]){        
      return 1;
    }
    else if(peptide_one_sequence[current_idx] < peptide_two_sequence[current_idx]){
      return -1;
    }
    ++current_idx;
  }
  
  // alphabetically identical, check if same length
  if(peptide_one_length > peptide_two_length){
    return 1;
  }
  else if(peptide_one_length < peptide_two_length){
    return -1;
  }
  else{
    return 0;
  }
}

/**
 * compares two peptides with the length sort type
 * /returns 1 if peptide_one has lower priority, 0 if equal, -1 if greater priority
 */
int Peptide::compareLengthQSort(
  Peptide** peptide_one, ///< peptide to compare one -in
  Peptide** peptide_two ///< peptide to compare two -in
  )
{
  int peptide_one_length = (*peptide_one)->length_;
  int peptide_two_length = (*peptide_two)->length_;
  
  // alphabetically identical, check if same length
  if(peptide_one_length > peptide_two_length){
    return 1;
  }
  else if(peptide_one_length < peptide_two_length){
    return -1;
  }
  else{
    return 0;
  }
}

/* Public functions--Printing / parsing */

/* text writing */

/**
 * \brief Prints a peptide object in text to file.
 *
 * Used by crux-generate-peptides. Prints the peptide once for all
 * peptide_src objects associated with it.  Optional fields are
 * determined by arguments. Tab-delimited format is:
 * mass protein-id peptide-start peptide-length <peptide-trypticity>
 * <peptide-sequence>
 * Peptide start begins with 1.
 */
void Peptide::printInFormat(
  bool flag_out, ///< print peptide sequence? -in
  FILE* file  ///< the out put stream -out
  ) {

  Protein* parent = NULL;
  char* id = NULL;
  int start_idx = 0;
  char* sequence = NULL;

  // obtain peptide sequence
  if(flag_out){
    sequence = getSequence();
  }

  for(vector<PeptideSrc*>::iterator iter=peptide_srcs_.begin();
      iter!=peptide_srcs_.end(); 
      ++iter
      ){
    PeptideSrc* next_src= *iter; 
    parent=next_src->getParentProtein();
     id = parent->getIdPointer();
    start_idx = next_src->getStartIdx();
    fprintf(file, "\t%s\t%d\t%d", id, start_idx, length_);
    // print peptide sequence?
    if(flag_out){
      fprintf(file, "\t%s\n", sequence);
    }
    else{
      fprintf(file, "\n");
    }
 }

  // free sequence if allocated
  if(flag_out){
    std::free(sequence);
  }
}

// TODO: this should be merged with other print, flags for optional
// fields and filtering should be taken from parameter.c and format
// adjusted accordingly
/**
 * Prints a peptide object to file.
 * ONLY prints peptide_src that match the peptide_src
 * mass \\t protein-id \\t peptide-start \\t peptide-length <\\t peptide-sequence> \n
 *      \\t protein-id \\t peptide-start \\t peptide-length <\\t peptide-sequence> \n
 * prints in correct format for generate_peptide
 */
void Peptide::printFilteredInFormat(
  bool flag_out, ///< print peptide sequence? -in
  FILE* file  ///< the out put stream -out
  ) {

  char* sequence = NULL;
  // bool light = false;

  // obtain peptide sequence
  if(flag_out){
    //parent = next_src->getParentProtein();
    
    // covnert to heavy protein
/*    
    FIXME, IF use light heavy put back
    if(get_protein_is_light(parent)){
      protein_to_heavy(parent);
      light = true;
    }
*/
    sequence = getSequence();
  }

  // iterate over all peptide src
/*
  while(next_src != NULL){
    if(peptide_type == ANY_TRYPTIC ||
       peptide_type == get_peptide_src_peptide_type(next_src) ||
       (peptide_type == PARTIALLY_TRYPTIC && 
        (get_peptide_src_peptide_type(next_src) == N_TRYPTIC ||
         get_peptide_src_peptide_type(next_src) == C_TRYPTIC)) ){
      
      // if(!light){
      parent = get_peptide_src_parent_protein(next_src);
        
      // covnert to heavy protein
      FIXME, IF use light heavy put back
      if(get_protein_is_light(parent)){
        protein_to_heavy(parent);
        light = true;
      }
        // }
      
      id = get_protein_id_pointer(parent);
      start_idx = get_peptide_src_start_idx(next_src);
      
      fprintf(file, "\t%s\t%d\t%d", id, start_idx, peptide->length);
      
      // print peptide sequence?
      if(flag_out){
        fprintf(file, "\t%s\n", sequence);
      }
      else{
        fprintf(file, "\n");
      }
    
       * uncomment this code if you want to restore a protein to 
       * light after converted to heavy
      // convert back to light
      if(light){
        protein_to_light(parent);
        light = false;
      }
    }
    next_src = get_peptide_src_next_association(next_src);
  }
*/

  // free sequence if allocated
  if(flag_out){
    std::free(sequence);
  }
}


/* binary read/writing (serialization) */

/**
 * Fills the given vectors with the names and descriptions of all
 * proteins containing this peptide.  Makes the descriptions
 * xml-friendly by swapping double for single quotes and angled braces
 * for square. Returned in the same order as getFlankingAAs().  Clears
 * any existing values in the vectors.
 * Adapted from Match::get_information_of_proteins()
 * \returns The number of proteins.
 */
int Peptide::getProteinInfo(vector<string>& protein_ids,
                            vector<string>& protein_descriptions){

  protein_ids.clear();
  protein_descriptions.clear();

  for(PeptideSrcIterator iter = getPeptideSrcBegin(); 
      iter!=getPeptideSrcEnd();
      ++iter
   ){
    PeptideSrc* peptide_src = *iter; 
    Protein* protein = peptide_src->getParentProtein();
    protein_ids.push_back(protein->getIdPointer());

    string description = "";
    const char* annotation_pointer = protein->getAnnotationPointer();
    if (annotation_pointer != NULL) {
      description = protein->getAnnotationPointer();
      // replace double quotes with single quotes
      replace(description.begin(), description.end(), '"', '\'');
      // remove any xml tags in the description by replacing <> with []
      replace(description.begin(), description.end(), '<', '[');
      replace(description.begin(), description.end(), '>', ']');
    }
    protein_descriptions.push_back(description);

  } 

  return protein_ids.size();
}



/* Public functions--Iterators */

/**
 * Instantiates a new residue_iterator from a peptide.
 * \returns a RESIDUE_ITERATOR_T object.
 */
RESIDUE_ITERATOR_T* new_residue_iterator(
  Peptide* peptide ///< peptide sequence to iterate -in
  )
{
  RESIDUE_ITERATOR_T* residue_iterator =
    (RESIDUE_ITERATOR_T*)calloc(1, sizeof(RESIDUE_ITERATOR_T));
  
  residue_iterator->peptide =  peptide;
  residue_iterator->residue_idx = 0;
  residue_iterator->sequence = peptide->getSequence();
  return residue_iterator;
}        

/**
 * Frees an allocated residue_iterator object.
 */
void free_residue_iterator(
  RESIDUE_ITERATOR_T* residue_iterator ///< free this object -in
  )
{
  free(residue_iterator->sequence);
  free(residue_iterator);
}

/**
 * The basic iterator functions.
 * \returns true if there are additional residues to iterate over, false if not.
 */
bool residue_iterator_has_next(
  RESIDUE_ITERATOR_T* residue_iterator ///< the query iterator -in
  )
{
  return (residue_iterator->residue_idx < residue_iterator->peptide->getLength());
}

/**
 * \returns The next residue (a character) in the peptide.
 */
char residue_iterator_next(
  RESIDUE_ITERATOR_T* residue_iterator  ///< the query iterator -in
  )
{
  ++residue_iterator->residue_idx;
  return residue_iterator->sequence[residue_iterator->residue_idx - 1];
}

/**
 * \brief Builds a comma delimited string listing the 
 * protein id(peptide start index) for the sources of 
 * a peptide
 *
 * \returns a string of the protein sources for this peptide
 */
string Peptide::getProteinIdsLocations() {

  set<string> protein_ids_locations;

 
  std::ostringstream protein_field_stream;
  if (!peptide_srcs_.empty()) {
    for( PeptideSrcIterator iter = getPeptideSrcBegin();
      iter!=getPeptideSrcEnd();++iter){
      
      PeptideSrc* peptide_src =*iter;
      Protein* protein = peptide_src->getParentProtein();
      char* protein_id = protein->getId();
      std::ostringstream protein_loc_stream;
      protein_loc_stream << protein_id;

      if (!protein->isPostProcess()) {
        int peptide_loc = peptide_src->getStartIdx();
        
        protein_loc_stream << "(" << peptide_loc << ")";
      } else if (peptide_src->getStartIdxOriginal() > 0) {
        protein_loc_stream << "(" << peptide_src->getStartIdxOriginal() << ")";
      }
      std::free(protein_id);
      protein_ids_locations.insert(protein_loc_stream.str());
    }
  }


  set<string>::iterator result_iter = protein_ids_locations.begin();
  string protein_field_string = *result_iter;

  while(++result_iter != protein_ids_locations.end()) {
    protein_field_string += "," + *result_iter;
  }

  return protein_field_string;
}

 
/**
 * \brief Builds a comma delimited string listing the protein ids
 * for the sources of a peptide.
 *
 * \returns a pointer to the string. Caller is responsible for freeing memeory.
 * If peptide has no sources returns NULL.
 */
char* Peptide::getProteinIds() {

  char *protein_field = NULL;

  if (!peptide_srcs_.empty()) {

    // Peptide has at least one parent.
    
    PeptideSrcIterator iter = getPeptideSrcBegin();
    PeptideSrc* peptide_src = *iter;
    Protein* protein = peptide_src->getParentProtein();

    const int allocation_factor = 1;
    char* protein_id = protein->getId();
    size_t protein_id_len = strlen(protein_id);
    size_t protein_field_len = allocation_factor * protein_id_len + 1; // Buffer size
    protein_field = (char*)malloc(sizeof(char) * protein_field_len);
    size_t protein_field_free = protein_field_len;  // Remaining free buffer space
    char *protein_field_tail = protein_field;
    *protein_field = 0;

    // First protein id in list doesn't have leading ','

    strncpy(protein_field_tail, protein_id, protein_field_free);
    protein_field_tail += protein_id_len;
    protein_field_free -= protein_id_len;
    delete protein_id;

    // Following proteins in list have leading ','
    for (;iter != getPeptideSrcEnd(); ++iter) {
      peptide_src = *iter;
      protein = peptide_src->getParentProtein();
      protein_id = protein->getId();
      protein_id_len = strlen(protein_id);

      // Allocate more memory if needed, allow space for comma and null
      if (protein_field_free < (protein_id_len + 2)) {
        size_t tail_offset = static_cast<std::size_t>(protein_field_tail - protein_field);
        protein_field = (char*)realloc(
          protein_field, 
          sizeof(char) * ((allocation_factor * (protein_id_len + 1)) + protein_field_len)
        );
        protein_field_len += allocation_factor * (protein_id_len + 1);
        protein_field_free += allocation_factor * (protein_id_len + 1);
        protein_field_tail = protein_field + tail_offset;
      }
    
      *protein_field_tail = ',';
      ++protein_field_tail;
      --protein_field_free;
      strncpy(protein_field_tail, protein_id, protein_field_free);
      protein_field_tail += protein_id_len;
      protein_field_free -= protein_id_len;
      delete protein_id;
    }
  }
  //carp(CARP_DEBUG, "Peptide::getProteinIds(): %s", protein_field);
  return protein_field;
}


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

