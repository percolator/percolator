/**********************************************
 *** adapted from crux/src/c/PeptideSrc.cpp ***
 **********************************************/

/*************************************************************************//**
 * \file PeptideSrc.cpp
 * \brief Object for mapping a peptide to its parent protein.
 ****************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <string>

#include "objects.h"
#include "Peptide.h"
#include "Protein.h"
#include "PeptideSrc.h"
#include "PeptideConstraint.h"

using namespace std;
using namespace PercolatorCrux;

/**
 * Static variable definitions
 */
map<string, Peptide* > PeptideSrc::sequence_to_peptide_; ///< Maps a sequence to a peptide object
map<string, Peptide* > PeptideSrc::decoy_sequence_to_peptide_; ///< Maps a decoy sequence to a peptide object


/**
 * \returns An (empty) peptide_src object.
 */
PeptideSrc::PeptideSrc() {
  digestion_ = (DIGEST_T)0;
  parent_protein_ = NULL;
  start_idx_ = 0;
  start_idx_original_ = 0;
}

/**
 *\returns a PeptideSrc object, populated with user specified parameters
 */
PeptideSrc::PeptideSrc(
  DIGEST_T digest,
  Protein* parent_protein, ///< the parent of this peptide -in
  int start_idx ///< start index of the peptide in the protein sequence -in
  ) {

  setDigest(digest);
  setParentProtein(parent_protein);
  setStartIdx(start_idx);
}


/**
 * Frees the entire allocated peptide_srcs object
 */
void PeptideSrc::free(vector<PeptideSrc*>& peptide_srcs) {


  for(vector<PeptideSrc*>::iterator iter = peptide_srcs.begin();
     iter != peptide_srcs.end();
     ++iter) {
  delete *iter; 

  }
  peptide_srcs.clear();
}

/**
 * Frees the an individual allocated peptide_src object
 * assumes that new_association pointer is NULL or some other pointer exist for the rest of the linklist 
 */
PeptideSrc::~PeptideSrc() {
  ;
}

// FIXME might need to change how this is printed
/**
 * Prints a peptide object to file.
 */
/*
void print_peptide_src(
  PEPTIDE_SRC_T* peptide_src, ///< object to print -in 
  FILE* file  ///< the out put stream -out
  )
{
  char* sequence = get_protein_sequence(peptide_src->parent_protein);
  fprintf(file, "parent protein:%s\n", sequence);

  fprintf(file, "peptide start: %d\n", peptide_src->start_idx);

  if(peptide_type == TRYPTIC){
    fprintf(file, "peptide type:%s\n", "TRYPTIC");
  }
  else if(peptide_type == PARTIALLY_TRYPTIC){
    fprintf(file, "peptide type:%s\n", "PARTIALLY_TRYPTIC");
  }
  else if(peptide_type == N_TRYPTIC){
    fprintf(file, "%s", "N_TRYPTIC");
  }
  else if(peptide_type == C_TRYPTIC){
    fprintf(file, "%s", "C_TRYPTIC");
  }
  else if(peptide_type == NOT_TRYPTIC){
    fprintf(file, "peptide type:%s\n", "NOT_TRYPTIC");
  }
  else if(peptide_type == ANY_TRYPTIC){
    fprintf(file, "peptide type:%s\n", "ANY_TRYPTIC");
  }
  free(sequence);
}
*/

/**
 * Copies the entire linklist of peptide_src object src to dest.
 * dest must be a heap allocated peptide_src
 */
void PeptideSrc::copy(
  vector<PeptideSrc*>& src, ///< source peptide_src -in
  vector<PeptideSrc*>& dest ///< destination peptide_src -out
  )
{

  for (vector<PeptideSrc*>::iterator iter = src.begin();
       iter != src.end();
       ++iter) {
    PeptideSrc* dest_src = new PeptideSrc(*(*iter));
    dest.push_back(dest_src);

  }
}

/**
 * sets the level of digestion
 */
void PeptideSrc::setDigest(
  DIGEST_T digest ///< the type of the peptide -in
  ){

  digestion_ = digest;
}

/**
 * \returns the level of digestion
 */
DIGEST_T PeptideSrc::getDigest() {

  return digestion_;
}

/**
 * sets the parent protein
 */
void PeptideSrc::setParentProtein(
  Protein* parent_protein ///< the parent of this preptide -in  
  ) {

  parent_protein_ = parent_protein;
}

/**
 * \returns a pointer to the parent protein
 */
Protein* PeptideSrc::getParentProtein() {

  return parent_protein_;
}
/**
 * sets the start index of the peptide in the protein sequence
 */
void PeptideSrc::setStartIdx(
  int start_idx ///< start index of the peptide in the protein sequence -in
  ) {

  start_idx_ = start_idx;
}

/**
 * \returns the start index of the peptide in the protein sequence
 */
int PeptideSrc::getStartIdx() {

  return start_idx_;
}

/**
 * \sets the original start index of the peptide in the protein sequence
 */
void PeptideSrc::setStartIdxOriginal(
  int start_idx ///< start index of the peptide in the original protein sequence -in
) {
  start_idx_original_ = start_idx;
}

/**
 * \returns the original start index of the peptide in the protein sequence
 */
int PeptideSrc::getStartIdxOriginal() {
  return start_idx_original_;
}

/**
 * \returns a pointer to the start of the peptide with in it's parent protein sequence
 */
char* PeptideSrc::getSequencePointer() {

  return parent_protein_->getSequencePointer(start_idx_ - 1);

}

/**
 *\returns the peptide_src strct size, value of sizeof function
 */
int PeptideSrc::getSizeOf(){
  return sizeof(PeptideSrc);
}


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
