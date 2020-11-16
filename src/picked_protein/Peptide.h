/*****************************************
 *** adapted from crux/src/c/Peptide.h ***
 *****************************************/
 
 /**
 * \file peptide.h 
 * $Revision: 1.52 $
 * \brief Object for representing one peptide.
 */
#ifndef CRUX_PEPTIDE_H 
#define CRUX_PEPTIDE_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <set>
#include <string>

#include "objects.h"
#include "Protein.h"
#include "Database.h"

namespace PercolatorCrux {

//these may be elsewhere
static const int MAX_PEPTIDE_LENGTH = 255;

/**
 * \class peptide
 * \brief A subsequence of a protein.
 */

class Peptide {

 protected:

  /**
   * static global variable
   * determines if the peptide src are created by link lists or array
   * if true, peptides are implented with link list peptide src, else array
   */

  /* Private data types */

  int pointer_count_;
  unsigned char length_; ///< The length of the peptide
  std:: vector<PeptideSrc*> peptide_srcs_; ///< a vector of peptide_srcs_

  /**
   * Initializes an (empty) peptide object
   */
  void init();
  
  // MT: helper function copied from crux/src/c/crux-utils.cpp
  /**
   * Returns copy of the src string upto the specified length.
   * Includes a null terminating character.
   * The string is heap allocated; thus, user must free.
   */
  char* copy_string_part(const char* src, int length){
    char* copy = (char*)calloc(length+1, sizeof(char));
    strncpy(copy, src, length);
    copy[length] = '\0';
    return copy;
  }
 public:

  /*  Allocators/deallocators  */
  
  /**
   * \returns An (empty) peptide object.
   */
  Peptide();

  /**
   * \returns A new peptide object, populated with the user specified
   * parameters.
   */
  Peptide(
    unsigned char length,     ///< The length of the peptide -in
    PercolatorCrux::Protein* parent_protein, ///< The parent_protein of this peptide -in
    int start_idx ///< Start index of peptide in the protein sequence -in
    );

  /**
   * \brief Allocates a new peptide giving it the values of the source
   * peptide.
   * \returns A newly allocated peptide identical to the source.
   */
  Peptide(
    Peptide* src ///< source peptide -in
  );

  /**
   *\returns the protein struct size, value of sizeof function
   */
  int getSizeOf(void);
 
  /**
   * Merge to identical peptides, copy all peptide_src into one of the peptide
   * peptide_dest, peptide_bye must have at least one peptide src
   * frees the peptide_bye, once the peptide_src are re-linked to the peptide_dest
   * Assumes that both peptides use linklist implemenation for peptide_src
   * \returns true if merge is successful else false
   */
  static bool mergePeptides(
    Peptide* peptide_dest,
    Peptide* peptide_bye
    );
                             
  /**
   * Merges two identical peptides by adding the peptide_src of the
   * second to the first.  The second peptide remains unchanged.
   * Does not comfirm identity of peptides.
   * \returns true if merge is successfull.
   */
  static bool mergePeptidesCopySrc(
    Peptide* peptide_dest,
    Peptide* peptide_giver
    );

  /**
   * Frees an allocated peptide object.
   * Depending on peptide_src implementation determines how to free srcs
   * This decision is made by global variable PEPTIDE_SRC_USE_LINK_LIST
   */
  ~Peptide();

  static void free(PercolatorCrux::Peptide* peptide);
  PercolatorCrux::Peptide* copyPtr();

  /*  Getters and Setters  */
  
  /*  Get-set:  source */
  
  /**
   * sets the peptide_src field in the peptide
   * must pass on a heap allocated peptide_src object
   * does not copy in the object, just the pointer to the object.
   */
  void setPeptideSrc(
    PeptideSrc* new_association ///< new peptide_src -in
    );

  /**
   * this method adds the new_association to the end of the existing peptide's 
   * linklist of peptide_srcs
   * must pass on a heap allocated peptide_src object
   * does not copy in the object, just the pointer to the object.
   */
  void addPeptideSrc(
    PeptideSrc* new_association ///< new peptide_src -in
    );

  /**
   * this method adds the peptide src array to an EMPTY peptide
   * only used in index.c, when the peptide src count for  peptide is known
   * Any existing peptide_src will lose it's reference
   */
  void addPeptideSrcArray(
    PeptideSrc* peptide_src_array ///< new peptide_src -in
    );

  /**
   * returns a pointer to the first PeptideSrc object of the peptide
   */
  PeptideSrc* getPeptideSrc();

  /**
   * returns a point to the peptide_protein_association field of the peptide
   */
  std::vector<PeptideSrc*>& getPeptideSrcVector();

  PeptideSrcIterator getPeptideSrcBegin();
  PeptideSrcIterator getPeptideSrcEnd();

  /**
   * get the peptide->first peptide_src->parent protein->database
   */
  Database* getFirstSrcDatabase();

  /**
   * \returns The number of peptide sources (i.e. proteins) the peptide has.
   */
  int getNumPeptideSrc();

  /**
   * returns a pointer to the peptide's first parent protein field of the peptide
   */
  PercolatorCrux::Protein* getParentProtein();

  /*  Get-set:  sequence */

  /**
   * sets the sequence length of the peptide
   */
  void setLength(
    unsigned char length  ///< the length of sequence -in
    );

  /**
   *\returns the sequence length of the peptide
   */
  unsigned char getLength();

  /**
   * \brief Get the sequence of a peptide.
   * Goes to the first peptide_src to gain sequence, thus must have at
   * least one peptide src 
   * \returns A newly allocated copy of the sequence.
   */
  char* getSequence();

  /**
   * \brief Get a string representation of the target (unshuffled)
   * peptide sequence with no added modification symbols.
   * For target peptides, returns the same as get_peptide_sequence.
   * \returns The newly-allocated sequence of peptide
   */
  char* getUnshuffledSequence();

  /**
   * \returns a pointer to the start of peptide sequence with in it's protein parent sequence, 
   * thus does not have terminating signe until end of parent protein
   * goes to the first peptide_src to find the location of start, thus must have at least one peptide src
   * should not print, will result in printing the entire protein sequence
   */
  char* getSequencePointer();

  /**
   * \brief Formats the sequence of the peptide with each flanking AA.
   * 
   * Format is "X.peptide_sequence.X", where "X" is a flanking amino acid.
   * "X", is printed as "-" if there is no flanking sequence.
   * Goes to the first peptide_src to gain sequence, thus must have at
   * least one peptide src 
   * \returns A newly allocated string with the sqt-formated peptide sequence.
   */
  char* getSequenceSqt();

  /**
   * \brief Formats the sequence of the peptide from a particular
   * peptide_src.
   *
   * Is called by get_peptide_sequence_sqt()
   * Format is "X.peptide_sequence.X", where "X" is a flanking amino acid.
   * "X", is printed as "-" if there is no flanking sequence.
   * Goes to the first peptide_src to gain sequence, thus must have at
   * least one peptide src 
   *
   * \returns A newly allocated string with the sqt-formated peptide sequence.
   */
  char* getSequenceFromPeptideSrcSqt(
    PeptideSrc* peptide_src ///< peptide_src -in 
   );

  /**
   * \brief Return a char for the amino acid c-terminal to the peptide
   * in the peptide src at the given index.
   *
   * \returns A char (A-Z) or - if peptide is the first in the protein.
   */
  char getCTermFlankingAA();

  /**
   * \brief Return a char for the amino acid n-terminal to the peptide
   * in the peptide src at the given index.
   *
   * \returns A char (A-Z) or - if peptide is the last in the protein.
   */
  char getNTermFlankingAA();

  bool isModified();

  bool isDecoy();

  /*  Getters requiring calculation */
  
  /**
   * Examines the peptide sequence and counts how many tryptic missed
   * cleavage sites exist. 
   *\returns the number of missed cleavage sites in the peptide
   */
  int getMissedCleavageSites();

  int getMissedCleavageSites(
    std::set<int> skip //skip these amino acid indices.
  );

  /**
   * Creates a heap allocated hash_value for the peptide that should
   * uniquely identify the peptide
   *\returns the string of "<first src protein idx><start idx><length>"
   */
  char* getHashValue();

  /**
   * Change the given target peptide into a decoy by randomizing its sequence.
   * Uses settings in parameter.c to decide between shuffling and
   * reversing the sequence.  Any modifications that exist will be
   * maintained on the same amino acids whose position will move.
   */
  void transformToDecoy();

  /*  Comparisons for sorting  */
  
  /**
   * Compare peptide sequence
   * \returns true if peptide sequence is identical else false
   */
  static bool compareSequence(
    Peptide* peptide_one,
    Peptide* peptide_two
  );

  /**
   * Compare two peptide sequences.
   * \returns Zero (0) if the sequences are identical, -1 if the first
   * sequence is less than the first and 1 if the first sequence is
   * greater than teh first.
   */
  static int triCompareSequence(
    Peptide* peptide_one,  ///< the peptide sequence to compare  -out
    Peptide* peptide_two  ///< the peptide sequence to compare  -out
    );

  /**
   * Compare the sequence of two peptides and return true if the first
   * petpide sequence is less than (in a lexical sort) the second peptide.
   */
  static bool lessThan(
    Peptide* peptide_one,
    Peptide* peptide_two
    );

  /**
   * compares two peptides with the lexical sort type
   * for qsort
   * /returns 1 if peptide_one has lower priority, 0 if equal, -1 if greater priority
   */
  static int compareLexicalQSort(
    Peptide** peptide_one, ///< peptide to compare one -in
    Peptide** peptide_two ///< peptide to compare two -in
    );

  /**
   * compares two peptides with the length sort type
   * /returns 1 if peptide_one has lower priority, 0 if equal, -1 if greater priority
   */
  static int compareLengthQSort(
    Peptide** peptide_one, ///< peptide to compare one -in
    Peptide** peptide_two ///< peptide to compare two -in
    );

  /*  Printing / parsing       */

  /**
   * Prints a peptide object to file.
   * prints all peptide_src object it's associated 
   * mass \\t protein-id \\t peptide-start \\t peptide-length <\\t peptide-trypticity> <\\t peptide-sequence> \n
   *      \\t protein-id \\t peptide-start \\t peptide-length <\\t peptide-trypticity> <\\t peptide-sequence> \n
   * prints in correct format for generate_peptide
   */
  void printInFormat(
    bool flag_out, ///< print peptide sequence? -in
    FILE* file  ///< the out put stream -out
    );

  /**
   * Prints a peptide object to file.
   * ONLY prints peptide_src that match the peptide_src
   * mass \\t protein-id \\t peptide-start \\t peptide-length <\\t peptide-sequence> \n
   *      \\t protein-id \\t peptide-start \\t peptide-length <\\t peptide-sequence> \n
   * prints in correct format for generate_peptide
   */
  void printFilteredInFormat(
    bool flag_out, ///< print peptide sequence? -in
    FILE* file  ///< the out put stream -out
    );

  /**
   * Serialize a peptide to a FILE in binary
   * \returns true if serialization is successful, else false
   *
   * The peptide serialization format looks like this:
   *
   *<Peptide: peptide struct><int: number of peptide_src>[<int: protein index><PEPTIDE_TYPE_T: peptide_type><int: peptide start index>]+
   * the bracket peptide src information repeats for the number of peptide src listed before the bracket
   * the protein index is the index of the parent protein in the database Database
   *
   */
  bool serialize(
    FILE* file,
    FILE* text_file
    );
 
  /**
   * \brief Builds a comma delimited string listing the 
   * protein id(peptide start index) for the sources of 
   * a peptide
   *
   * \returns a string of the protein sources for this peptide
   */
  std::string getProteinIdsLocations();

  /**
   * \brief Builds a comma delimited string listing the protein ids
   * for the sources of a peptide.
   *
   * \returns a pointer to the string. Caller is responsible for freeing memeory.
   * If peptide has no sources returns NULL.
   */
  char *getProteinIds();

  /**
   * \brief Builds a comma delimited string listing the flanking amino acids
   * for the sources of a peptide.
   *
   * \returns a pointer to the string. Caller is responsible for freeing memeory.
   * If peptide has no sources returns NULL.
   */
  char* getFlankingAAs();

  /**
   * Fills the given vectors with the names and descriptions of all
   * proteins containing this peptide.  Returned in the same order as
   * getFlankingAAs().  Clears any existing data in the vectors.
   * \returns The number of proteins.
   */
  int getProteinInfo(std::vector<std::string>& protein_ids,
                     std::vector<std::string>& protein_descriptions);

};  // class Peptide

/*  Iterators */

/**
 * Instantiates a new residue_iterator from a peptide.
 * \returns a RESIDUE_ITERATOR_T object.
 */
RESIDUE_ITERATOR_T* new_residue_iterator(
  PercolatorCrux::Peptide* peptide ///< peptide sequence to iterate -in
  );

/**
 * Frees an allocated residue_iterator object.
 */
void free_residue_iterator(
  RESIDUE_ITERATOR_T* residue_iterator ///< free this object -in
  );

/**
 * The basic iterator functions.
 * \returns true if there are additional residues to iterate over, false if not.
 */
bool residue_iterator_has_next(
  RESIDUE_ITERATOR_T* residue_iterator ///< the query iterator -in
  );

/**
 * \returns The next residue (a character) in the peptide.
 */
char residue_iterator_next(
  RESIDUE_ITERATOR_T* residue_iterator  ///< the query iterator -in
  );

} // end namespace PercolatorCrux

#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
