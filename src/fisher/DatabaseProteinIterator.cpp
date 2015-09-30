/**
 * \file DatabaseProteinIterator.cpp
 * $Revision: 1.27 $
 * \brief Object to iterator over the proteins within a database
 *****************************************************************************/

#include "DatabaseProteinIterator.h"
using namespace Crux;

/**
 * Instantiates a new database_protein_iterator from a database.
 * 
 * \returns a DATABASE_PROTEIN_ITERATOR_T object.
 */
DatabaseProteinIterator::DatabaseProteinIterator(
  Database* database ///< the database to create a protein iterator -in
  )
{
  // if database is parsed, if not do so..
  if(!database->getIsParsed()){
    // failed to parse database
    if(!database->parse()){
      //carp(CARP_FATAL, "Failed to parse database, cannot create iterator");
    }
  }

  database_ = Database::copyPtr(database);
  cur_protein_ = 0;
}        


/**
 * Frees an allocated database_protein_iterator object.
 */
DatabaseProteinIterator::~DatabaseProteinIterator() {

  // subtract pointer count
  Database::freeDatabase(database_);
}

/**
 * The basic iterator functions.
 * \returns true if there are additional proteins to iterate over, false if not.
 */
bool DatabaseProteinIterator::hasNext()
{
  return (cur_protein_ < database_->getNumProteins());
}


/**
 * \returns The next protein in the database.
 */
Protein* DatabaseProteinIterator::next()
{
  ++cur_protein_;

  // print number of protein generated to STDERR for every 500 protein reached
  if(cur_protein_ % 500 == 0){
    //carp(CARP_DETAILED_DEBUG, "Reached protein %d out of %d", cur_protein_, database_->getNumProteins());
  }
  
  return database_->getProteinAtIdx(cur_protein_-1);
}

/**
 * \returns the protein database for this iterator
 */
Database* DatabaseProteinIterator::getDatabase() {

  return database_;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
