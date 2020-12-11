/********************************************
 *** adapted from crux/src/c/PeptideSrc.h ***
 ********************************************/

/**
 * \file PeptideSrch
 * $Revision: 1.12 $
 * \brief Object for mapping a peptide to it's parent protein.
 */
#ifndef PEPTIDE_SRC_H
#define PEPTIDE_SRC_H
#include <map>
#include <stdio.h>
#include <string>

#include "objects.h"

namespace PercolatorCrux {

class PeptideSrc {

 protected:
  DIGEST_T digestion_; ///< how specific the ends are relative to the enzyme
  PercolatorCrux::Protein* parent_protein_; ///< the parent of this preptide

  /*
   * start_idx_ tracks the index of the peptide within the corresponding protein
   * start_idx_original_ is for keeping track of the index of the peptide within
   * the original protein, which we may not have the full sequence of
   */
  int start_idx_; ///< start index of the peptide in the protein sequence, first residue is 1 
  int start_idx_original_;  ///< start index of the peptide in the original protein sequence
  static std::map<std::string, PercolatorCrux::Peptide* > sequence_to_peptide_; ///< Maps a sequence to a peptide object
  static std::map<std::string, PercolatorCrux::Peptide* > decoy_sequence_to_peptide_; ///< Maps a decoy sequence to a peptide object

  /**
   * \brief fills the sequence_to_peptide_ member variable for use in parseTabDelimited
   * used when the tab delimited file doesn't provide a protein id, but we have
   * sequences and access to the database.
   */
  static void fillPeptides(
    Database* database, ///< the protein database 
    Database* decoy_database ///< the decoy database
    );

 public:

  /**
   * \returns An (empty) peptide_src object.
   */
  PeptideSrc();

  /**
   * \returns a PeptideSrc object, populated with user
   * specified parameters 
   */
  PeptideSrc(
    DIGEST_T digest,
    PercolatorCrux::Protein* parent_protein, ///< the parent of this preptide -in
    int start_idx ///< peptide start index in protein sequence, first is 1 -in
    );

  /**
   * Frees the entire allocated peptide_src object
   */
  static void free(
    std::vector<PeptideSrc*>& peptide_srcs ///< object to free -in 
  );

  /**
   * Frees the an individual allocated peptide_src object
   */
  virtual ~PeptideSrc();


  /**
   * Prints a peptide object to file.
   */
  /*
  void print(
    FILE* file  ///< the output stream -out
    );
  */

  /**
   * Copies the entire vector of peptide_src object src to dest.
   * dest must be a heap allocated peptide_src
   */
  static void copy(
    std::vector<PeptideSrc*>& src, ///< source peptide_src -in
    std::vector<PeptideSrc*>& dest ///< destination peptide_src -out
    );

  /**
   * sets the level of digestion
   */
  void setDigest(
    DIGEST_T digest ///< the type of the peptide -in
    );

  /**
   * \returns the level of digestion
   */
  DIGEST_T getDigest();

  /**
   * sets the parent protein
   */
  void setParentProtein(
    PercolatorCrux::Protein* parent_protein ///< the parent of this preptide -in  
    );

  /**
   * \returns a pointer to the parent protein
   */
  PercolatorCrux::Protein* getParentProtein();

  /**
   * sets the start index of the peptide in the protein sequence
   */
  void setStartIdx(
    int start_idx ///< start index of the peptide in the protein sequence -in
    );

  /**
   * \returns the start index of the peptide in the protein sequence
   */
  int getStartIdx();

  /**
   * sets the original start index of the peptide in the protein sequence
   */
  void setStartIdxOriginal(
    int start_idx ///< start index of the peptide in the original protein sequence -in
  );

  /**
   * \returns the original start index of the peptide in the protein sequence
   */
  int getStartIdxOriginal();

  /**
   * \returns a pointer to the start of the peptide with in it's parent protein sequence
   */
  char* getSequencePointer();

  /**
   *\returns the peptide_src strct size, value of sizeof function
   */
  static int getSizeOf(void);

};
/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
 
} // end namespace PercolatorCrux
#endif
