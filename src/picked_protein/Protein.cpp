/*******************************************
 *** adapted from crux/src/c/Protein.cpp ***
 *******************************************/

/*************************************************************************//**
 * \file protein.cpp
 * \brief Object for representing a single protein.
 ****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>
#include <vector>

#include "Protein.h"
#include "Peptide.h"
#include "PeptideSrc.h"
#include "PeptideConstraint.h"

using namespace std;
using namespace PercolatorCrux;
/**
 * Constants
 */
static const int PROTEIN_ID_LENGTH = 100;
static const int PROTEIN_SEQUENCE_LENGTH = 40000;
static const int PROTEIN_ANNOTATION_LENGTH = 100;
static const int LONGEST_LINE  = PROTEIN_ID_LENGTH + PROTEIN_ID_LENGTH;
static const int FASTA_LINE = 50;
static const int SMALLEST_MASS = 57;
static const int LARGEST_MASS = 190;

/**
 * Set member variables to default values.
 */
void Protein::init() {
  database_ = NULL;
  offset_ = 0;
  protein_idx_ = 0;
  is_light_ = false;
  is_memmap_ = false;
  id_ = NULL;
  sequence_ = NULL;
  length_ = 0;
  annotation_ = NULL;
}

/**
 * \returns An (empty) protein object.
 */
Protein::Protein() {
  init();
}

/**
 * \returns A new protein object(heavy).
 * The protein is does not constain a database, users must provide one.
 */
Protein::Protein(
  const char*         id, ///< The protein sequence id. -in
  const char*   sequence, ///< The protein sequence. -in
  unsigned int length, ///< The length of the protein sequence. -in
  const char* annotation,  ///< Optional protein annotation.  -in
  unsigned long int offset, 
  ///< The file location in the source file in the database -in
  unsigned int protein_idx, ///< The index of the protein in it's database.-in  
  Database* database ///< the database of its origin
  )
{
  init();
  setId(id);
  setSequence(sequence);
  setLength(length);
  setOffset(offset);
  setProteinIdx(protein_idx);
  setIsLight(false);
  database_ = Database::copyPtr(database); 
  is_memmap_ = false;
}         




/**
 * \returns A new light protein object.
 */
Protein* Protein::newLightProtein(
  unsigned long int offset, 
  ///< The file location in the source file in the database -in
  unsigned int protein_idx ///< The index of the protein in it's database. -in
  )
{
  Protein* protein = new Protein();
  protein->setIsLight(true);
  protein->setOffset(offset);
  protein->setProteinIdx(protein_idx);
  return protein;
}

/**
 * Frees an allocated protein object.
 */
Protein::~Protein() 
{
  // FIXME what is the point of checking this?
  if(!is_memmap_ && !is_light_){ 
    if (id_ != NULL){
      free(id_);
    }
    if (sequence_ != NULL){
      free(sequence_);
    }
    if (annotation_ != NULL){
      free(annotation_);
    }
  }
}

/**
 * Prints a protein object to file.
 * if light protein coverts it to heavy protein
 */
void Protein::print(
  FILE* file ///< output stream -out
  )
{
  int   sequence_index;
  int   sequence_length = getLength();
  char* sequence = getSequence();
  char* id = getId();
  char* annotation = getAnnotation();
  
  fprintf(file, ">%s %s\n", id, annotation);

  sequence_index = 0;
  while (sequence_length - sequence_index > FASTA_LINE) {
    fprintf(file, "%.*s\n", FASTA_LINE, &(sequence[sequence_index]));
    sequence_index += FASTA_LINE;
  }
  fprintf(file, "%s\n\n", &(sequence[sequence_index]));

  free(sequence);
  free(id);
  free(annotation);
}

/**
 * \returns the starting location of the sequence in a protein.  If not found, returns -1
 */
int Protein::findStart(
  string peptide_sequence, ///< the sequence to find
  string prev_aa, ///< the previous amino acid for the sequence
  string next_aa ///< the next amino acid for the sequence
  ) {

  if (getSequencePointer() == NULL) {
    return -1;
  }

  if (prev_aa == "-") {
    return 1;
  } else if (next_aa == "-") {
    return getLength() - peptide_sequence.length();
  } else {
    //use the flanking amino acids to further constrain our search in the sequence
    size_t pos = string::npos; 
    string seq = prev_aa + peptide_sequence + next_aa;
    string protein_seq = getSequencePointer();
    pos = protein_seq.find(seq);
    if (pos == string::npos) {
      //carp(CARP_DEBUG, "could not find %s in protein %s\n%s", seq.c_str(), getIdPointer(), protein_seq.c_str());
      //finding the sequence with the flanks failed, try finding without the flanks.
      seq = peptide_sequence;
      pos = protein_seq.find(seq);
      if (pos == string::npos) {
        //carp(CARP_ERROR, "could not %s in protein %s\n%s", seq.c_str(), getIdPointer(), protein_seq.c_str());
        return -1;
      }
      return (pos+1);
    }
    return (pos+2);
  }

}

bool Protein::isPostProcess() {

  return false;

}

/**
 * prints a binary representation of the protein
 * 
 * FORMAT (no line break)
 * <int: id length><char: id><int: annotation length>
    <char: annotation><int: sequence length><char: sequence>
 *
 * Make sure when reading the binary data, add one to the length so
 * that it will read in the terminating char as well.
 */
void Protein::serialize(
  FILE* file ///< output stream -out
  )
{
  
  int id_length = strlen(id_);
  int annotation_length = strlen(annotation_);

  // write the protein id length
  fwrite(&id_length, sizeof(int), 1, file);
  
  // write the protein id 
 // include "/0"
  fwrite(id_, sizeof(char), id_length+1, file);

  // write the protein annotation length
  fwrite(&annotation_length, sizeof(int), 1, file);

  // write the protein annotation
  // include "/0"
  fwrite(annotation_, sizeof(char), annotation_length+1, file);
  
  // write the protein sequence length
  fwrite(&length_, sizeof(unsigned int), 1, file);
  
  // write the protein sequence
  // include "/0"
  fwrite(sequence_, sizeof(char), length_+1, file);
}


/**
 * Copies protein object src to dest.
 * assumes that the protein is heavy
 * dest must be a heap allocated object 
 */
void Protein::copy(
  Protein* src,///< protein to copy -in
  Protein* dest ///< protein to copy to -out
  )
{
  char* id = src->getId();
  char* sequence = src->getSequence();
  char* annotation = src->getAnnotation();
  
  dest->setId(id);
  dest->setSequence(sequence);
  dest->setLength(src->getLength());
  dest->setOffset(src->offset_);
  dest->setProteinIdx(src->protein_idx_);
  dest->setIsLight(src->is_light_);
  dest->database_ = src->database_;
  
  free(id);
  free(sequence);
  free(annotation);
}


/**
 * Parses a protein from an memory mapped binary fasta file
 * the protein_idx field of the protein must be added before or after
 * you parse the protein 
 * protein must be a heap allocated
 * 
 * Assume memmap pointer is set at beginning of protein
 * Assume protein binary format (no line break)
 * <int: id length><char: id><int: annotation length>
     <char: annotation><int: sequence length><char: sequence>
 *
 * modifies the *memmap pointer!
 * \returns TRUE if success. FALSE is failure.
 */
bool Protein::parseProteinBinaryMemmap(
  char** memmap 
  ///< a pointer to a pointer to the memory mapped binary fasta file -in
  )
{
  int id_length = 0;
  int annotation_length = 0;
  int sequence_length = 0;

  /* FIXME, maybe use this to check if still within file
  if(*memmap_as_char[0] == EOF){
    //carp(CARP_ERROR, "end of file");
  }
  */

  /***set protein ID***/

  // read id length
  id_length = *((int *) *memmap);

  // reset pointer to start of id
  *memmap += sizeof(int);

  // set protein id to mem mapped id
  id_ = *memmap;

  // reset pointer to move to annotation_length
  *memmap += (id_length + 1);


  /***set protein annotation***/

  // read annotation length
  annotation_length = *((int *) *memmap);

  // reset pointer to start of annotation
  *memmap += sizeof(int);

  // set protein annotation to mem mapped annotation
  annotation_ = *memmap;

  // reset pointer to move to sequence_length
  *memmap += (annotation_length + 1);


  /***set protein sequence***/
  
  // read sequence length
  sequence_length = *((int *) *memmap);
  length_ = sequence_length;

  // reset pointer to start of sequence
  *memmap += sizeof(int);

  // set protein annotation to mem mapped sequence
  sequence_ = *memmap;

  // reset pointer to move to start of next protein
  *memmap += (sequence_length + 1);
  
  // now this protein has been created from memory mapped!
  is_memmap_ = true;

  return true;
}

/**
 * Parses a protein from an open (FASTA) file.
 * the protein_idx field of the protein must be added before or after
 * you parse the protein  
 * \returns TRUE if success. FALSE is failure.
 * protein must be a heap allocated
 */
bool Protein::parseProteinFastaFile(
  FILE* file ///< fasta file -in
  )
{
  static char name[LONGEST_LINE];    ///< Just the sequence ID.
  static char desc[LONGEST_LINE];    ///< Just the comment field.
  static char buffer[PROTEIN_SEQUENCE_LENGTH];///> The sequence to read in.
  static unsigned int sequence_length; // the sequence length

  // Read the title line.
  if (!readTitleLine(file, name, desc)) {
    return(false);
  }
  
  buffer[0] = '\0';

  // Read the sequence.
  if (!readRawSequence(file, name, PROTEIN_SEQUENCE_LENGTH, buffer, &sequence_length)) {
    //carp(CARP_FATAL, "Sequence %s is too long.\n", name);
    exit(1);
  }

  // update the protein object.
  setLength(sequence_length);
  setId(name);
  setSequence(buffer);

  return(true);

}

/**************************************************/

/**
 * FASTA file parsing code
 * AUTHOR: William Stafford Noble
 * modified by Chris Park
 */

/**
 * Find the beginning of the next sequence, and read the sequence ID
 * and the comment.
 */
bool Protein::readTitleLine
(FILE* fasta_file, ///< file to read -in
 char* name, ///< write protein name here -out
 char* description) ///< write description here -out
{
  static char id_line[LONGEST_LINE];  // Line containing the ID and comment.
  int a_char;                         // The most recently read character.

  // Read until the first occurrence of ">".
  while ((a_char = getc(fasta_file)) != '>') {
    // If we hit the end of the file, return FALSE.
    if (a_char == EOF) {
      return(false);
    }  
  }
  // set protein offset                   FIXME: might not need to "-1" -CHRIS
  offset_ = ftell(fasta_file) - 1;

  /**
   * chris edited, added this block to make sure all of comment line
   * is read although might not be stored, to ensure the file* is at
   * start of the sequence
   */

  {
    char* new_line = NULL;
    int line_length;
    size_t buf_length = 0;

    if((line_length =  getline(&new_line, &buf_length, fasta_file)) == -1){
      //carp(CARP_FATAL, "Error reading Fasta file.\n");
    }
    strncpy(id_line, new_line, LONGEST_LINE-1);
    free(new_line);
  }

  // Remove EOL.
  id_line[strlen(id_line) - 1] = '\0';

  // Extract the ID from the beginning of the line.
  if (sscanf(id_line, "%s", name) != 1) {
    //carp(CARP_FATAL, "Error reading sequence ID.\n%s\n", id_line);
  }

  // Store the rest of the line as the comment.
  strcpy(description, &(id_line[strlen(name)+1]));

  return(true);
}


/****************************************************************************
 * Read raw sequence until a '>' is encountered or too many letters
 * are read.  The new sequence is appended to the end of the given
 * sequence.
 *
 * Return: Was the sequence read completely?
 ****************************************************************************/
bool Protein::readRawSequence
  (FILE* fasta_file,   // Input Fasta file.
   char* name,         // Sequence ID (used in error messages).
   unsigned int   max_chars,    // Maximum number of allowed characters.
   char* raw_sequence, // Pre-allocated sequence.
   unsigned int* sequence_length // the sequence length -chris added
   )
{
  int a_char;
  unsigned int i_seq;
  bool return_value = true;

  // Start at the end of the given sequence.
  i_seq = strlen(raw_sequence);
  assert((unsigned int)strlen(raw_sequence) < max_chars);

  // Read character by character.
  while ((a_char = getc(fasta_file)) != EOF) {

    // Check for the beginning of the next sequence.
    if (a_char == '>') {
      // Put the ">" back onto the stream for the next call to find.
      ungetc(a_char, fasta_file);
      break;
    }

    // Skip non-alphabetic characters.
    if (!isalpha((int)a_char)) {
      if ((a_char != ' ') && (a_char != '\t') && (a_char != '\n') && (a_char != '\r')) {
        //carp(CARP_WARNING,"Skipping character %c in sequence %s.", a_char, name);
      }

    } else {

      // Convert invalid characters to X.
      a_char = toupper((int)a_char);

      /**
       * Check the ASCII code.  If the char is above or below the
       * A(65)~Z(90)range, convert the character to an 'X'.
       */
      if ( (int)a_char < 65 || (int)a_char  > 90 ) {
        //carp(CARP_WARNING, "Converting illegal character %c to X ", a_char);
        //carp(CARP_WARNING, "in sequence %s.", name);
        a_char = 'X';
      }
      
      raw_sequence[i_seq] = a_char;
      i_seq++;
    }
    if (i_seq >= max_chars) {
      return_value = false;
      break;
    }
  }
  raw_sequence[i_seq] = '\0';
  *sequence_length = i_seq; // chris added

  return(return_value);
}


/**
 * end of FASTA parsing
 * Thanks Bill!
 */

/**
 * Change the sequence of a protein to be a randomized version of
 * itself.  The method of randomization is dependant on the
 * decoy_type (shuffle or reverse).  The name of the protein is also
 * changed by prefixing with "decoy-prefix"
 */
void Protein::shuffle(
  DECOY_TYPE_T decoy_type){ ///< method for shuffling

  switch(decoy_type){
  case NO_DECOYS:
    return;

  case PROTEIN_REVERSE_DECOYS:
    reverse(sequence_, sequence_ + strlen(sequence_));
    break;

  case PROTEIN_SHUFFLE_DECOYS:
  case PEPTIDE_SHUFFLE_DECOYS:
  case INVALID_DECOY_TYPE:
  case NUMBER_DECOY_TYPES:
    //carp(CARP_FATAL, "Illegal decoy type for shuffling protein.");
    break;
  default:
    exit(EXIT_FAILURE);
  }
}

/** 
 * Access routines of the form get_<object>_<field> and set_<object>_<field>. 
 */

/**
 * Additional get and set methods
 */

/**
 *\returns the id of the protein
 * returns a heap allocated new copy of the id
 * user must free the return id
 * assumes that the protein is heavy
 */
char* Protein::getId()
{

  if(is_light_){
    //carp(CARP_FATAL, "Cannot get ID from light protein.");
  }
  
  int id_length = strlen(id_) +1; // +\0
  char* copy_id = 
    (char *)malloc(sizeof(char)*id_length);
  
  strncpy(copy_id, id_, id_length); 

  return copy_id;
}

/**
 *\returns a pointer to the id of the protein
 * assumes that the protein is heavy
 */
char* Protein::getIdPointer()
{
  if(is_light_){
    //carp(CARP_FATAL, "Cannot get ID pointer from light protein.");
  }
  return id_; 
}

/**
 * sets the id of the protein
 */
void Protein::setId(
  const char* id ///< the sequence to add -in
  )
{
  int id_length = strlen(id) +1; // +\0
  free(id_);
  char* copy_id = 
    (char *)malloc(sizeof(char)*id_length);
  id_ =
    strncpy(copy_id, id, id_length);  
}

/**
 *\returns the sequence of the protein
 * returns a heap allocated new copy of the sequence
 * user must free the return sequence 
 * assumes that the protein is heavy
 */
char* Protein::getSequence(
  int offset
) {
  if(is_light_){
    //carp(CARP_FATAL, "Cannot get sequence from light protein.");
  }
  unsigned int sequence_length = strlen(sequence_) +1-offset; // +\0
  char * copy_sequence = 
    (char *)malloc(sizeof(char)*sequence_length);
  return strncpy(copy_sequence, sequence_+offset, sequence_length);  
}

/**
 *\returns a pointer to the sequence of the protein
 * assumes that the protein is heavy
 */
char* Protein::getSequencePointer(
  int offset
)
{
  if(is_light_){
    //carp(CARP_FATAL, "Cannot get sequence pointer from light protein.");
  }
  return sequence_+offset;
}

/**
 * sets the sequence of the protein
 */
void Protein::setSequence(
  const char* sequence ///< the sequence to add -in
  )
{

  unsigned int sequence_length = strlen(sequence) +1; // +\0
  free(sequence_);
  char * copy_sequence = 
    (char *)malloc(sizeof(char)*sequence_length);
  sequence_ =
    strncpy(copy_sequence, sequence, sequence_length);  
}

/**
 *\returns the length of the protein
 * assumes that the protein is heavy
 */
unsigned int Protein::getLength()
{
  return length_;
}

/**
 * sets the id of the protein
 */
void Protein::setLength(
  unsigned int length ///< the length to add -in
  )
{
  length_ = length;
}

/**
 *\returns the annotation of the protein
 * returns a heap allocated new copy of the annotation
 * user must free the return annotation
 * assumes that the protein is heavy
 */
char* Protein::getAnnotation()
{
  if(is_light_){
    //carp(CARP_FATAL, "Cannot get annotation from light protein.");
  }
  int annotation_length = strlen(annotation_) +1; // +\0
  char * copy_annotation = 
    (char *)malloc(sizeof(char)*annotation_length);
  return strncpy(copy_annotation, annotation_, annotation_length);  
}

/**
 *\returns A const pointer to the annotation of the protein.
 */
const char* Protein::getAnnotationPointer(){
  return annotation_;
}

/**
 * sets the offset of the protein in the fasta file
 */
void Protein::setOffset(
  unsigned long int offset 
  ///< The file location in the source file in the database -in
  )
{
  offset_ = offset;
}

/**
 *\returns the offset the protein
 */
unsigned long int Protein::getOffset()
{
  return offset_;
}

/**
 * sets the protein_idx (if, idx=n, nth protein in the fasta file)
 */
void Protein::setProteinIdx(
  unsigned int protein_idx ///< The index of the protein in it's database. -in
  )
{
  // //carp(CARP_DETAILED_DEBUG, "set protein idx = %i", protein_idx);
  protein_idx_ = protein_idx;
}

/**
 *\returns the protein_idx field
 */
unsigned int Protein::getProteinIdx()
{
  return protein_idx_;
}

/**
 * sets the is_light field (is the protein a light protein?)
 */
void Protein::setIsLight(
  bool is_light ///< is the protein a light protein? -in
  )
{
  is_light_ = is_light;
}

/**
 *\returns TRUE if the protein is light protein
 */
bool Protein::getIsLight()
{
  return is_light_;
}

/**
 * sets the database for protein
 */
void Protein::setDatabase(
  Database*  database ///< Which database is this protein part of -in
  )
{
  database_ = Database::copyPtr(database);
}

/**
 *\returns Which database is this protein part of
 */
Database* Protein::getDatabase()
{
  return database_;
}

/** 
 * Comparison function for sorting proteins by protein id.
 */
bool PercolatorCrux::protein_id_less_than(Protein* protein_one, Protein* protein_two){
  int compare = strcmp(protein_one->getIdPointer(),
                       protein_two->getIdPointer());
  return (compare > 0);
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

