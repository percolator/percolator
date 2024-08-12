/*****************************************************
 *** adapted from crux/src/c/PeptideConstraint.cpp ***
 *****************************************************/

/*************************************************************************//**
 * \file PeptideConstraint.cpp
 * \brief Object for holding the peptide constraint information.
 ****************************************************************************/
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "PeptideConstraint.h"

using namespace PercolatorCrux;

/**
 * Allocates a new (empty) peptide_constraint object.
 * \returns An allocated PEPTIDE_CONSTRAINT_T object.
 */
void PeptideConstraint::init() {
  //carp(CARP_DETAILED_DEBUG, "Initializing peptide constraint");
  enzyme_ = (ENZYME_T)0;
  digestion_ = (DIGEST_T)0;
  min_length_ = 0;
  max_length_ = 0;
  num_mis_cleavage_ = 0;
  num_pointers_ = 1;

}

PeptideConstraint::PeptideConstraint() {
  init();
}

/**
 * Instantiates a new peptide_constraint object.
 * \returns An allocated PEPTIDE_CONSTRAINT_T object.
 */
PeptideConstraint::PeptideConstraint(
  ENZYME_T enzyme, 
  DIGEST_T digest,
  int min_length, ///< the minimum length of peptide -in
  int max_length,  ///< the maximum lenth of peptide(max limit = 255) -in
  int num_mis_cleavage ///< The maximum mis cleavage of the peptide -in
  )
{
  // max length must be less or equal than 255 becuase of the unsigned char limit of 255
  if(max_length > 255){
    //carp(CARP_FATAL, "ERROR: cannot set max length higer than 255");
  }
  
  init();

  setEnzyme(enzyme);
  setDigest(digest);
  setMinLength(min_length);
  setMaxLength(max_length);
  setNumMisCleavage(num_mis_cleavage);
}


/**
 * Copy peptide pointer and increment pointer count
 */
PeptideConstraint* PeptideConstraint::copyPtr(
  PeptideConstraint* constraint
  ) {

  constraint->num_pointers_++;
  return constraint;
}

// FIXME check the association..as long as there is one tryptic parent then true
// num_miss_cleavage is not implemented..add if needed
/** 
 * Determines if a peptide satisfies a peptide_constraint.
 * \returns TRUE if the constraint is satisified. FALSE if not.
 */
bool PeptideConstraint::isSatisfied(
  Peptide* peptide ///< the query peptide -in
  ) {

  return (peptide->getLength() <= getMaxLength() &&
     peptide->getLength() >= getMinLength()
     );
}

/**
 * Frees an allocated peptide_constraint object.
 */
void PeptideConstraint::free(
  PeptideConstraint* peptide_constraint ///< object to free -in 
  ) {

  peptide_constraint->num_pointers_--;
  if (peptide_constraint->num_pointers_ <= 0) {
    //carp(CARP_DETAILED_DEBUG, "Final free of peptide constraint");
    delete peptide_constraint;
  }

}

PeptideConstraint::~PeptideConstraint() {

}

/**
 * Sets the enzyme used for the in silicos digestion
 * of the protein sequence into peptides.
 */
void PeptideConstraint::setEnzyme(
  ENZYME_T enzyme
){

  enzyme_ = enzyme;
}

/**
 * \returns The enzyme for this peptide constraint.
 */
ENZYME_T PeptideConstraint::getEnzyme()
{
  return enzyme_;
}

/**
 * Sets the level of digestion for the peptide constraint.
 */
void PeptideConstraint::setDigest(
  DIGEST_T digest
  ){

  digestion_ = digest;
}

/**
 * \returns The level of digestion for the peptide constraint.
 */
DIGEST_T PeptideConstraint::getDigest() {

  return digestion_;
}

/**
 * sets the min length of the peptide_constraint
 */
void PeptideConstraint::setMinLength(
  int min_length  ///< the min length of the peptide constraint - in
  ) {

  min_length_ = min_length;
}

/**
 * \returns the min length of the peptide_constraint
 */
int PeptideConstraint::getMinLength() {

  return min_length_;
}

/**
 * sets the max length of the peptide_constraint
 * maximum limit 255
 */
void PeptideConstraint::setMaxLength(
  int max_length  ///< the max length of the peptide constraint - in
  ) {

  // check if maximum length is with in range <= 255
  if(max_length > 255){
    //carp(CARP_FATAL, "maximum length:%d over limit 255.", max_length);
  }
  
  max_length_ = max_length;
}

/**
 * \returns the max length of the peptide_constraint
 */
int PeptideConstraint::getMaxLength(
  ) {

  return max_length_;
}


/**
 * sets the num_mis_cleavage of the peptide_constraint
 */
void PeptideConstraint::setNumMisCleavage(
  int num_mis_cleavage ///< The maximum mis cleavage of the peptide -in
  ) {

  num_mis_cleavage_ = num_mis_cleavage;
}

/**
 * \returns the num_mis_cleavage of the peptide_constraint
 */
int PeptideConstraint::getNumMisCleavage() {

  return num_mis_cleavage_;
}
