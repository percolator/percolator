/**
 * \file DatabaseProteinIterator.h 
 * $Revision: 1.27 $
 * \brief Object to iterator over the proteins within a database
 *****************************************************************************/

#ifndef DATABASEPROTEINITERATOR_H
#define DATABASEPROTEINITERATOR_H

#include "Database.h"

class DatabaseProteinIterator {
  
 protected:
  Database* database_;  ///< The database whose proteins to iterate over 
  unsigned int cur_protein_;      ///< The index of the current protein

 public:
  
  /**
   * Instantiates a new database_protein_iterator from a database.
   * \returns a DATABASE_PROTEIN_ITERATOR_T object.
   */
  DatabaseProteinIterator(
    Database* database ///< the database to create a protein iterator -in
    );        

  /**
   * Frees an allocated database_protein_iterator object.
   */
  virtual ~DatabaseProteinIterator();
  
  /**
   * The basic iterator functions.
   * \returns TRUE if there are additional proteins to iterate over, FALSE if not.
   */
  bool hasNext();
  
  /**
   * \returns The next protein in the database.
   */
  Crux::Protein* next();
  
  /**
   * \returns the protein to the corresponding protein_idx in the database.
   */
  Crux::Protein* getProtein(
    unsigned int protein_idx ///< protein_idx to which protein to return -in
    );

  /**
   * \returns the protein database for this iterator
   */
  Database* getDatabase();

};

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
