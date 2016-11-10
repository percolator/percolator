/********************************************************
 *** adapted from crux/src/c/ProteinPeptideIterator.h ***
 ********************************************************/

/**
 * \file ProteinPeptideIterator.h 
 * $Revision: 1.25 $
 * \brief Object to iterate over the peptides within a protein in an
 * unspecified order. The peptides should satisfy the constraints specified
 * in the peptide_constraint object.
 *****************************************************************************/
#ifndef PROTEINPEPTIDEITERATOR_H
#define PROTEINPEPTIDEITERATOR_H

#include "Protein.h"
#include "Peptide.h"
#include "PeptideConstraint.h"

#include <iterator>
#include <vector>

namespace PercolatorCrux {

class ProteinPeptideIterator {

 protected:
  PercolatorCrux::Protein* protein_; ///< The protein whose peptides to iterate over. 
  unsigned short int cur_start_; ///< Start in protein of the current peptide.
  unsigned short int cur_length_; ///< The length of the current peptide.
  unsigned int peptide_idx_; ///< The index of the current peptide.
  PeptideConstraint* peptide_constraint_; ///< peptide type to iterate over.
  double* mass_array_; ///< stores all the peptides' masses
  std::vector<int>* nterm_cleavage_positions_; ///< nterm cleavages that satisfy 
                                        ///< constraint. 1st aa is 1.
  std::vector<int>* peptide_lengths_; ///< all the lengths of valid peptides
  std::vector<FLOAT_T>* peptide_masses_; ///< all the masses of valid peptides
  std::vector<int>* cumulative_cleavages_; ///< cumulative number of cleavages so far
  int current_cleavage_idx_; /// where are we in the cleavage positions?
  int num_cleavages_; /// how many cleavage positions?

  bool has_next_; ///< is there a next? 
  int num_mis_cleavage_; ///< The maximum mis cleavage of the peptide

  /*
   * Takes a cumulative distribution of peptide masses (the mass_array) and
   * the start index and end index and returns a peptide mass
   */
  FLOAT_T calculateSubsequenceMass(
    double* mass_array,
    int start_idx,
    int cur_length
  );

  /**
   * \brief Decide if a residue is in an inclusion list or is not in an
   * exclusion list. 
   *
   * For use with the user-specified enzyme digestion.  Takes an amino
   * acid, a list of amino acids, and a flag for if it is an inclusion
   * list or an exclusion list.  A cleavage can happen before/after the
   * given residue if it is either in the inclusion list or is not in
   * the exculsion list.
   * \returns TRUE if the residue is in the inclusion list or not in the
   * exclusion list.
   */
  static bool isResidueLegal(
    char aa, 
    char* aa_list, 
    int list_size, 
    bool for_inclusion
  );

  /**
   * \brief Adds cleavages to the protein peptide iterator that obey iterator
   * constraint.
   *
   * Uses the allowed cleavages arrays, and whether skipped cleavages
   * are allowed. 
   * A small inconsistency: 
   *  Allowed cleavages start at 0, while the output cleavages start at 1.
   */
  void selectPeptides(
    int* nterm_allowed_cleavages, 
    int  nterm_num_cleavages, 
    int* cterm_allowed_cleavages, 
    int  cterm_num_cleavages, 
    int  num_skip_cleavages);

  /**
   * Creates the data structures in the protein_peptide_iterator object necessary
   * for creating peptide objects.
   * - mass_array - cumulative distribution of masses. used to determine 
   *     the mass of any peptide subsequence.
   * - nterm_cleavage_positions - the nterm cleavage positions of the 
   *     peptides that satisfy the protein_peptide_iterator contraints
   * - peptide_lengths - the lengths of the peptides that satisfy the constraints
   * - cumulative_cleavages - cumulative distribution of cleavage positions
   *    used to determine if a cleavage location has been skipped
   */
  void prepare();

  /**
   * \brief Estimate the maximum number of peptides a protein can
   * produce.  Counts the number of subsequences of length
   * min_seq_length, min_seq_length + 1, ..., max_seq_length that can be
   * formed from a protein of the given length.  No enzyme specificity
   * assumed.  
   */
  unsigned int countMaxPeptides(
    unsigned int protein_length,   ///< length of protein
    unsigned int min_seq_length,   ///< min peptide length
    unsigned int max_seq_length);  ///< max peptide length

 public:
  /**
   * Instantiates a new peptide_iterator from a peptide.
   * \returns a PROTEIN_PEPTIDE_ITERATOR_T object.
   */
  ProteinPeptideIterator(
    PercolatorCrux::Protein* protein,
    PeptideConstraint* peptide_constraint
  );

  /**
   * Frees an allocated peptide_iterator object.
   */
  ~ProteinPeptideIterator();

  /**
   * The basic iterator functions.
   * \returns TRUE if there are additional peptides to iterate over, FALSE if not.
   */
  bool hasNext();

  /**
   * \returns The next peptide in the protein, in an unspecified order
   */
  PercolatorCrux::Peptide* next();

  /**
   *\returns the protein that the iterator was created on
   */
  PercolatorCrux::Protein* getProtein();

  /**
   * \returns The total number of peptides in this protein.
   */
  int getTotalPeptides();

  /**
   * Creates the data structures in the protein_peptide_iterator object necessary
   * for creating peptide objects.
   * - mass_array - cumulative distribution of masses. used to determine 
   *     the mass of any peptide subsequence.
   * - nterm_cleavage_positions - the nterm cleavage positions of the 
   *     peptides that satisfy the protein_peptide_iterator contraints
   * - peptide_lengths - the lengths of the peptides that satisfy the constraints
   * - cumulative_cleavages - cumulative distribution of cleavage positions
   *    used to determine if a cleavage location has been skipped
   */
  void prepareMc(
    int missed_cleavages
  );

  /**
   * Compares the first and second amino acids in the given sequence to
   * see if they conform to the cleavage rules of the given enzyme.  For
   * NO_ENZYME, always returns TRUE.
   *
   * \returns TRUE if this is a valid cleavage position for the given enzyme.
   */
  static bool validCleavagePosition(
    const char* sequence,
    ENZYME_T enzyme
  );

};

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

}; // end namespace PercolatorCrux

#endif
