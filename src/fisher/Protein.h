/*****************************************
 *** adapted from crux/src/c/Protein.h ***
 *****************************************/

/**
 * \file Protein.h 
 * $Revision: 1.25 $
 * \brief Object for representing one protein sequence.
 *****************************************************************************/
#ifndef PROTEIN__H 
#define PROTEIN__H

#include <stdio.h>
#include <string>
#include <cstring>
#include <algorithm>
//#include "utils.h"
#include "objects.h"
//#include "Peptide.h"
//#include "PeptideSrc.h"
//#include "carp.h"
//#include "PeptideConstraint.h"

namespace Crux {

class Protein {
 protected:
  Database*  database_; ///< Which database is this protein part of
  unsigned long int offset_; ///< The file location in the database source file
  unsigned int protein_idx_; ///< The index of the protein in it's database.
  bool    is_light_; ///< is the protein a light protein?
  bool    is_memmap_; ///< is the protein produced from memory mapped file
  char*              id_; ///< The protein sequence id.
  char*        sequence_; ///< The protein sequence.
  unsigned int   length_; ///< The length of the protein sequence.
  char*      annotation_; ///< Optional protein annotation.

  /**
   * Find the beginning of the next sequence, and read the sequence ID
   * and the comment.
   */
  bool readTitleLine
    (FILE* fasta_file,///< file to read
     char* name,      ///< put protein name here
     char* description);///< put protein description here

  
  /**
   * Read raw sequence until a '>' is encountered or too many letters
   * are read.  The new sequence is appended to the end of the given
   * sequence.
   *
   * Return: Was the sequence read completely?
   **/
  static bool readRawSequence
    (FILE* fasta_file,   // Input Fasta file.
     char* name,         // Sequence ID (used in error messages).
     unsigned int   max_chars,    // Maximum number of allowed characters.
     char* raw_sequence, // Pre-allocated sequence.
     unsigned int* sequence_length // the sequence length -chris added
   );
 
  /**
   * Rearrange the sequence_ between cleavage sites, keeping residues
   * on either side of a cleavage in place.  Get enzyme from
   * parameters.  Same behavior for full and partial digest, min/max
   * length/mass and missed cleavages, i.e. shuffle between every
   * cleavage site.
   */
  void peptideShuffleSequence();

  /**
   * Shuffle the region of the sequence between start and end, leaving
   * start and end residues in place.  Repeat up to three times if the
   * shuffled sequence doesn't change.
   */
  void shuffleRegion(
    int start, ///< index of peptide start
    int end);  ///< index of last residue in peptide
  
  /**
   * given two strings return a concatenated third string
   * \returns a heap allocated string that concatenates the two inputs
   */
  char* cat_string(const char* string_one, const char* string_two){
    int len_one = strlen(string_one);
    int len_two = strlen(string_two);
    
    char* result = (char*)calloc(len_one + len_two + 1, sizeof(char));
    strncpy(result, string_one, len_one);
    strncpy(&result[len_one], string_two, len_two);
    return result;
  }
 public:

  /**
   * Initialize member variables to default values.
   */
  void init();

  /**
   * \returns An (empty) protein object.
   */
  Protein();

  /**
   * \returns A new protein object(heavy).
   */
  Protein(
    const char*         id, ///< The protein sequence id.
    const char*   sequence, ///< The protein sequence.
    unsigned int length, ///< The length of the protein sequence.
    const char* annotation,  ///< Optional protein annotation.  -in
    unsigned long int offset, 
    ///< The file location in the source file in the database -in
    unsigned int protein_idx,///< The index of the protein in its database. -in
    Database* database ///< the database of its origin
  );         

  /**
   * \returns A new light protein object.
   */
  static Protein* newLightProtein(
    unsigned long int offset, 
    ///< The file location in the source file in the database -in
    unsigned int protein_idx ///< The index of the protein in its database. -in
  );

  /**
   * convert light protein to heavy, by parsing all the sequence from fasta file
   * \returns TRUE if successfully converts the protein to heavy 
   */
  bool toHeavy();

  /**
   * covert heavy protein back to light
   * \returns TRUE if successfully converts the protein to light
   */
  bool toLight();

  /**
   * Frees an allocated protein object.
   */
  
  virtual ~Protein();

  /**
   * Prints a protein object to file.
   */
  void print(
    FILE* file ///< output stream -out
  );

  /**
   * Copies protein object src to dest.
   * dest must be a heap allocated object 
   */
  static void copy(
    Protein* src,///< protein to copy -in
    Protein* dest ///< protein to copy to -out
  );

  /**
   * Parses a protein from an open (FASTA) file.
   * \returns TRUE if success. FALSE is failure.
   */
  bool parseProteinFastaFile(
    FILE* file ///< fasta file -in
  );

  /**
   * Parses a protein from an memory mapped binary fasta file
   * the protein_idx field of the protein must be added before or
   * after you parse the protein.
   * \returns TRUE if success. FALSE is failure.
   * protein must be a heap allocated
   * 
   * Assume memmap pointer is set at beginning of protein
   * Assume protein binary format (no line break)
   * <int: id length><char: id><int: annotation length>
     <char: annotation><int: sequence length><char: sequence>
   *
   * modifies the *memmap pointer!
   */
  bool parseProteinBinaryMemmap(
    char** memmap 
    ///< a pointer to a pointer to the memory mapped binary fasta file -in
  );

  /**
   * Change the sequence of a protein to be a randomized version of
   * itself.  The method of randomization is dependant on the
   * decoy_type (shuffle or reverse).  The name of the protein is also
   * changed by prefixing with reverse_ or rand_, depending on how it
   * was randomized. 
   */
  void shuffle(DECOY_TYPE_T decoy_type);

  /**
   * Additional get and set methods
   */

  /**
   *\returns the id of the protein
   * returns a heap allocated new copy of the id
   * user must free the return id
   */
  char* getId();

  /**
   *\returns a pointer to the id of the protein
   */
  char* getIdPointer();

  /**
   * sets the id of the protein
   */
  void setId(
    const char* id ///< the sequence to add -in
    );

  /**
   *\returns the sequence of the protein
   * returns a heap allocated new copy of the sequence
   * user must free the return sequence 
   */
  virtual char* getSequence(
    int offset=0
  );

  /**
   *\returns a pointer to the sequence of the protein
   */
  virtual char* getSequencePointer(
    int offset=0
  );

  /**
   * sets the sequence of the protein
   */
  void setSequence(
    const char* sequence ///< the sequence to add -in
  );

  /**
   *\returns the length of the protein
   */
  virtual unsigned int getLength();

  /**
   * sets the id of the protein
   */
  void setLength(
    unsigned int length ///< the length to add -in
  );

  /**
   *\returns the annotation of the protein
   * returns a heap allocated new copy of the annotation
   * user must free the return annotation
   */
  char* getAnnotation();

  /**
   *\returns A const pointer to the annotation of the protein.
   */
  const char* getAnnotationPointer();

   /**
   * sets the annotation of the protein
   */
  void setAnnotation(
    const char* annotation ///< the sequence to add -in
  );


  /**
   * sets the offset of the protein in the fasta file
   */
  void setOffset(
    unsigned long int offset 
    ///< The file location in the source file in the database -in
    );

  /**
   *\returns the offset the protein
   */
  unsigned long int getOffset();

  /**
   * sets the protein_idx (if, idx=n, nth protein in the fasta file)
   */
  void setProteinIdx(
    unsigned int protein_idx ///< The index of the protein in its database. -in
  );

  /**
   *\returns the protein_idx field
   */
  unsigned int getProteinIdx();

  /**
   * sets the is_light field (is the protein a light protein?)
   */
  void setIsLight(
    bool is_light ///< is the protein a light protein? -in
    );

  /**
   *\returns TRUE if the protein is light protein
   */
  bool getIsLight();

  /**
   * sets the database for protein
   */
  void setDatabase(
    Database*  database ///< Which database is this protein part of -in
    );

  /**
   *\returns Which database is this protein part of
   */
  virtual Database* getDatabase();

  /**
   * prints a binary representation of the protein
   * 
   * FORMAT (no line break)
   * <int: id length><char: id><int: annotation length>
     <char: annotation><int: sequence length><char: sequence>
   *
   * make sure when rading the binary data, add one to the length so
   * that it will read in the terminating char as well.
   */
  void serialize(
    FILE* file ///< output stream -out
    );

  /**
   * \returns the start index of the sequence in the protein
   * tries to use the previous aa and next aa to refine where
   * this sequence exits
   */
  virtual int findStart(
    std::string sequence,
    std::string prev_aa,
    std::string next_aa
  );

  virtual bool isPostProcess();

};

/** 
 * Comparison function for sorting proteins by protein id.
 */
bool protein_id_less_than(Protein* protein_one, Protein* protein_two);

};

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
