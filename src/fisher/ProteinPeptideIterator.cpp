/**********************************************************
 *** adapted from crux/src/c/ProteinPeptideIterator.cpp ***
 **********************************************************/

/**
 * \file ProteinPeptideIterator.cpp
 * $Revision: 0.0 $
 */

#include "ProteinPeptideIterator.h"

using namespace std;
using namespace Crux;


/**
 * \brief Decide if a residue is in an inclusion list or is not in an
 * exclusion list. 
 *
 * For use with the user-specified enzyme digestion.  Takes an amino
 * acid, a list of amino acids, and a flag for if it is an inclusion
 * list or an exclusion list.  A cleavage can happen before/after the
 * given residue if it is either in the inclusion list or is not in
 * the exculsion list.
 * \returns true if the residue is in the inclusion list or not in the
 * exclusion list.
 */
bool ProteinPeptideIterator::isResidueLegal(char aa, 
                           char* aa_list, 
                           int list_size, 
                           bool for_inclusion){

  // The logic for returning for_inclusion:
  // For an inclusion list (true), once we find the aa it passes (true)
  // For an exclusion list (false), once we find the aa, it fails (false)
  int idx=0;
  for(idx=0; idx < list_size; idx++){
    if( aa == aa_list[idx] ){ return for_inclusion; }
  }
  // or if we got to the end of the list and didn't find a match
  // for inclusion, it fails (!true)
  // for exclusion, it passes (!false)
  return ! for_inclusion;
}

/**
 * Compares the first and second amino acids in the given sequence to
 * see if they conform to the cleavage rules of the given enzyme.  For
 * NO_ENZYME, always returns true.
 *
 * \returns true if this is a valid cleavage position for the given enzyme.
 */
bool ProteinPeptideIterator::validCleavagePosition(
   const char* sequence,
   ENZYME_T enzyme
){

  switch(enzyme){

  case TRYPSIN:
    return (sequence[0] == 'K' || sequence[0] == 'R') && sequence[1] != 'P';
  //trypsin/p k or r 
  case TRYPSINP:
    return sequence[0] == 'K' || sequence[0] == 'R';
  case ELASTASE:
    return (sequence[0] == 'A' || sequence[0] == 'L' ||
            sequence[0] == 'I' || sequence[0] == 'V') && sequence[1] != 'P';
  case CLOSTRIPAIN:
    return sequence[0] == 'R';
  case CYANOGEN_BROMIDE:
    return sequence[0] == 'M';
  case IODOSOBENZOATE:
    return sequence[0] == 'W';
  case PROLINE_ENDOPEPTIDASE:
    return sequence[0] == 'P';
  case STAPH_PROTEASE:
    return sequence[0] == 'E';
  case ASPN:
    return sequence[1] == 'D';
  case LYSC:
    return sequence[0] == 'K' && sequence[1] != 'P';
  case LYSN:
    return sequence[1] == 'K';
  case ARGC:
    return sequence[0] == 'R' && sequence [1] != 'P';
  case GLUC:
    return (sequence[0] == 'D' || sequence[0] == 'E') && sequence[1] != 'P';
  case PEPSINA:
    return (sequence[0] == 'F' || sequence[0] == 'L') && sequence[1] != 'P';
  case CHYMOTRYPSIN:
    return (sequence[0] == 'F' || sequence[0] == 'L' ||
            sequence[0] == 'W' || sequence[0] == 'Y') && sequence[1] != 'P';
  case ELASTASE_TRYPSIN_CHYMOTRYPSIN:
    return (sequence[0] == 'A' || sequence[0] == 'L' ||
            sequence[0] == 'I' || sequence[0] == 'V' ||
            sequence[0] == 'K' || sequence[0] == 'R' ||
            sequence[0] == 'W' || sequence[0] == 'F' ||
            sequence[0] == 'Y' ) && sequence[1] != 'P';
  case CUSTOM_ENZYME:
    //carp(CARP_FATAL, "The custom enzyme is not yet implmented.");
    return false;
    /*
    return ( isResidueLegal(sequence[0], 
                              pre_cleavage_list,
                              pre_list_size, 
                              pre_for_inclusion)
             && 
             isResidueLegal(sequence[1], 
                              post_cleavage_list,
                              post_list_size, 
                              post_for_inclusion) );
                              */
    break;

  case NO_ENZYME:
    return true;
    break;

  case INVALID_ENZYME:
  case NUMBER_ENZYME_TYPES:
    //carp(CARP_FATAL, "Cannot generate peptides with invalid enzyme.");
    break;

  }// end switch

  return false;
}

/**
 * \brief Adds peptides to the iterator based on our constraint and
 * the given possible cleavage positions.
 * 
 * The allowed cleavages on either end of the peptide are specified
 * separately so that the ends can obey different cleavage rules
 * (e.g. tryptic and non).  Use the member variable
 * cumulative_cleavages_ to keep track of skipped enzyme cleavage
 * sites since the cterm allowed cleavages may be non-tryptic for a partially
 * tryptic search.  Add a peptide to the iterator for each
 * pair of n- and c-term cleavages that obey all of the peptide
 * constraints: correct length, mass, and number of internal cleavage
 * positions.
 * A small inconsistency: 
 *  Allowed cleavages start at 0, while the output cleavages start at 1.
 */
void ProteinPeptideIterator::selectPeptides(
    int* nterm_allowed_cleavages, 
    int  nterm_num_cleavages, 
    int* cterm_allowed_cleavages, 
    int  cterm_num_cleavages, 
    int  int_num_skip_cleavages){

  // to avoid checking a lot of C-term before our current N-term cleavage
  int previous_cterm_cleavage_start= 0;

  PeptideConstraint* constraint = peptide_constraint_;
  int nterm_idx, cterm_idx;

  // for each possible n-term (start) position...
  for (nterm_idx=0; nterm_idx < nterm_num_cleavages; nterm_idx++){

    // check all possible c-term (end) positions

    int next_cterm_cleavage_start = previous_cterm_cleavage_start;
    bool no_new_cterm_cleavage_start = true;
    for (cterm_idx = previous_cterm_cleavage_start; 
         cterm_idx < cterm_num_cleavages; cterm_idx++){

      if ((*cumulative_cleavages_)[cterm_allowed_cleavages[cterm_idx]-1] - \
          (*cumulative_cleavages_)[nterm_allowed_cleavages[nterm_idx]]  \
          > int_num_skip_cleavages) {
        break;
      }
      if (cterm_allowed_cleavages[cterm_idx] 
          <= nterm_allowed_cleavages[nterm_idx]){
        continue;
      }
      
      // check our length constraint
      int length = 
        cterm_allowed_cleavages[cterm_idx] - nterm_allowed_cleavages[nterm_idx];
      
      // if too short, try next cterm position
      if (length < constraint->getMinLength()){
        continue;
        // if too long, go to next nterm (start) position
      } else if (length > constraint->getMaxLength()){
        break;
      } else if (no_new_cterm_cleavage_start){
        next_cterm_cleavage_start = cterm_idx;
        no_new_cterm_cleavage_start = false;
      }
     
      // we have found a peptide
      nterm_cleavage_positions_->push_back(nterm_allowed_cleavages[nterm_idx] + 1);

      peptide_lengths_->push_back(length);
      //carp(CARP_DETAILED_DEBUG, "New pep: %i (%i)", nterm_allowed_cleavages[nterm_idx], length);

      num_cleavages_++;
    }
    previous_cterm_cleavage_start = next_cterm_cleavage_start;

  }
}

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
void ProteinPeptideIterator::prepare()
{
  // TODO: use global variable for this
  prepareMc(0);
}

void ProteinPeptideIterator::prepareMc(
    int missed_cleavages)
{
  Protein* protein = protein_;

  ENZYME_T enzyme = peptide_constraint_->getEnzyme();

  // initialize mass matrix and enzyme cleavage positions
  int* cleavage_positions = (int*) calloc(protein->getLength()+1, sizeof(int));
  int* non_cleavage_positions = (int*)calloc(protein->getLength()+1, sizeof(int));
  int* all_positions = (int*) calloc(protein->getLength()+1, sizeof(int));

  // initialize first value in all array except non_cleavage_positions
  unsigned int start_idx = 0;
  int cleavage_position_idx = 0;
  int non_cleavage_position_idx = 0;
  cleavage_positions[cleavage_position_idx++] = 0;

  // calculate our cleavage positions and masses
  for(start_idx = 1; start_idx < protein->getLength()+1; start_idx++){
    int sequence_idx = start_idx - 1;
    char amino_acid = protein->getSequencePointer()[sequence_idx];

    if( amino_acid == 'B' || amino_acid == 'X' || amino_acid == 'Z' ){
      //carp_once(CARP_WARNING, "Ignoring peptides with ambiguous amino acids (B, X, Z).");
    } 

    // increment cumulative cleavages before we check if current position
    // is a cleavage site because cleavages come *after* the current amino acid
    cumulative_cleavages_->push_back(cleavage_position_idx);

    //if (valid_cleavage_position(protein->sequence + sequence_idx)){ 
    if (validCleavagePosition(protein->getSequencePointer() + sequence_idx, enzyme)){ 
      cleavage_positions[cleavage_position_idx++] = sequence_idx + 1;
    } else {
      non_cleavage_positions[non_cleavage_position_idx++] = sequence_idx + 1;
    }

    all_positions[sequence_idx] = sequence_idx;
  }

  // put in the implicit cleavage at end of protein
  if (cleavage_positions[cleavage_position_idx-1] != (int)protein->getLength()){
    cleavage_positions[cleavage_position_idx++] = protein->getLength(); 
  }

  all_positions[protein->getLength()] = (int)protein->getLength();

  int num_cleavage_positions = cleavage_position_idx;
  int num_non_cleavage_positions = non_cleavage_position_idx;

  //carp(CARP_DETAILED_DEBUG, "num_cleavage_positions = %i", num_cleavage_positions);

  // now determine the cleavage positions that actually match our constraints

  DIGEST_T digestion = 
    peptide_constraint_->getDigest();

  switch (digestion){

  case FULL_DIGEST:
      this->selectPeptides(
        cleavage_positions, num_cleavage_positions-1,
        cleavage_positions+1, num_cleavage_positions-1, 
        missed_cleavages);

      break;

  case PARTIAL_DIGEST:
      // add the C-term tryptic cleavage positions.
      this->selectPeptides(
        all_positions, protein->getLength(),
        cleavage_positions+1, num_cleavage_positions-1, 
        missed_cleavages);

      // add the N-term tryptic cleavage positions.
      // no +1 below for non_cleavage_positions below 
      // because it does not include sequence beginning. it is *special*
      this->selectPeptides(
        cleavage_positions, num_cleavage_positions-1,
        non_cleavage_positions, num_non_cleavage_positions-1,
        missed_cleavages);

      break;

  case NON_SPECIFIC_DIGEST:
      this->selectPeptides(
        all_positions, protein->getLength(),
        all_positions+1, protein->getLength(), // len-1?
        500); // for unspecific ends, allow internal cleavage sites
      break;

  case INVALID_DIGEST:
  case NUMBER_DIGEST_TYPES:
    //carp(CARP_FATAL, "Invalid digestion type in protein peptide iterator.");
    break;
  }

/*
  int idx;
  for (idx=0; idx < iterator->num_cleavages; idx++){
    //carp(CARP_DETAILED_DEBUG, "%i->%i", 
         iterator->nterm_cleavage_positions[idx], 
         //iterator->peptide_lengths[idx], 
         //iterator->peptide_lengths[idx], 
       iterator->protein->sequence[iterator->nterm_cleavage_positions[idx]-1]);
  }
*/
  if (num_cleavages_ > 0){
    has_next_ = true;
  } else { 
    has_next_ = false;
  }

  free(cleavage_positions);
  free(non_cleavage_positions);
  free(all_positions);
}

/**
 * \brief Estimate the maximum number of peptides a protein can
 * produce.  Counts the number of subsequences of length
 * min_seq_length, min_seq_length + 1, ..., max_seq_length that can be
 * formed from a protein of the given length.  No enzyme specificity
 * assumed.  
 */
unsigned int ProteinPeptideIterator::countMaxPeptides(
 unsigned int protein_length,   ///< length of protein
 unsigned int min_seq_length,   ///< min peptide length
 unsigned int max_seq_length)  ///< max peptide length
{
  if( max_seq_length > protein_length ){
    max_seq_length = protein_length;
  }

  unsigned int total_peptides = 0;
  for(unsigned int len = min_seq_length; len <= max_seq_length; len++){
    total_peptides += protein_length + 1 - len;
  }
  return total_peptides;
}

/**
 * Instantiates a new peptide_iterator from a protein.
 * \returns a PROTEIN_PEPTIDE_ITERATOR_T object.
 * assumes that the protein is heavy
 */
ProteinPeptideIterator::ProteinPeptideIterator(
  Protein* protein, ///< the protein's peptide to iterate -in
  PeptideConstraint* peptide_constraint ///< the peptide constraints -in
  )
{

  // initialize iterator
  protein_ = NULL;
  cur_start_ = 0;
  cur_length_ = 1;
  peptide_idx_ = 0;
  peptide_constraint_ = NULL;
  current_cleavage_idx_ = 0;
  has_next_ = false;

  peptide_idx_ = 0;
  peptide_constraint_ 
    = PeptideConstraint::copyPtr(peptide_constraint);
  cur_start_ = 0; 
  cur_length_ = 1;  
  num_mis_cleavage_ 
    = peptide_constraint_->getNumMisCleavage();
  protein_ = protein;

  nterm_cleavage_positions_ = new vector<int>();
  peptide_lengths_ = new vector<int>();
  cumulative_cleavages_ = new vector<int>();

  // estimate array size and reserve space to avoid resizing vector
  // TODO: create global variable for these
  int max_peptides = countMaxPeptides(protein->getLength(), 6, 50);
  nterm_cleavage_positions_->reserve(max_peptides); 
  peptide_lengths_->reserve(max_peptides);
  cumulative_cleavages_->reserve(max_peptides);
  num_cleavages_ = 0;

  // prepare the iterator data structures
  prepare();
}


/**
 * Frees an allocated peptide_iterator object.
 */
ProteinPeptideIterator::~ProteinPeptideIterator() 
{
  PeptideConstraint::free(peptide_constraint_);
  delete nterm_cleavage_positions_; 
  delete peptide_lengths_; 
  delete cumulative_cleavages_; 
}

/**
 * The basic iterator functions.
 * \returns true if there are additional peptides, false if not.
 */
bool ProteinPeptideIterator::hasNext()
{
  return has_next_;
}

/**
 * \returns The next peptide in the protein, in an unspecified order
 * the Peptide is new heap allocated object, user must free it
 */
Crux::Peptide* ProteinPeptideIterator::next()
{
  if( !has_next_){
    //carp(CARP_DEBUG, "Returning null");
    return NULL;
  }

  int cleavage_idx = current_cleavage_idx_;
  int current_start = (*nterm_cleavage_positions_)[cleavage_idx];
  int current_length = (*peptide_lengths_)[cleavage_idx];

  // create new peptide
  Peptide* peptide = new Peptide(current_length, protein_, current_start);//, peptide_type);
  // update position of iterator
  ++current_cleavage_idx_;

  // update has_next field
  if (current_cleavage_idx_ == num_cleavages_){
    has_next_ = false;
  } else {
    has_next_ = true;
  }
  return peptide;
}

/**
 *\returns the protein that the iterator was created on
 */
Protein* ProteinPeptideIterator::getProtein()
{
  return protein_;
}

/**
 * \returns The total number of peptides in this protein.
 */
int ProteinPeptideIterator::getTotalPeptides(){
  return num_cleavages_;
}


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
