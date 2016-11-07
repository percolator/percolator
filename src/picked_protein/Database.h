/******************************************
 *** adapted from crux/src/c/Database.h ***
 ******************************************/

/**
 * \file Database.h 
 * $Revision: 1.27 $
 * \brief Object for representing a database of protein sequences.
 *****************************************************************************/
#ifndef DATABASE_H
#define DATABASE_H

#include <stdio.h>
#include "objects.h"
#include "Protein.h"
#include <string>
#include <cstring>
#include <map>
#include <vector>

namespace PercolatorCrux {

#ifdef _MSC_VER
//#include "WinCrux.h"
/*********************************************************
 This function replaces the GNU extension of the same name.
 Reads a line from the given stream.
 *********************************************************/
int getline(char **lineptr, size_t *n, FILE *stream);
#endif

//Comparator function for c type strings.
struct cmp_str {

  bool operator()(char const *a, char const *b) {
    return strcmp(a, b) < 0;
  }
};

class Database {
 protected:
  std::string fasta_filename_; ///< Name of the text file.
  std::string binary_filename_;///< Full path to the binary protein sequence file.
  FILE*        file_;     ///< Open filehandle for this database.
                         ///  A database has only one associated file.
  bool is_parsed_;  ///< Has this database been parsed yet.
  std::vector<PercolatorCrux::Protein*>* proteins_; ///< Proteins in this database.
  std::map<char*, PercolatorCrux::Protein*, cmp_str>* protein_map_; //map for proteins 
  bool is_hashed_; //Indicator of whether the database has been hashed/mapped.
  unsigned long int size_; ///< The size of the database in bytes (convenience)
  bool use_light_protein_; ///< should I use the light/heavy protein option
  bool is_memmap_; ///< Are we using a memory mapped fasta file? 
  void* data_address_; ///< pointer to the beginning of the memory mapped data, 
  unsigned int pointer_count_; ///< number of pointers referencing this database. 
  long file_size_; ///< the size of the binary fasta file, when memory mapping
  DECOY_TYPE_T decoys_; ///< the type of decoys, none if target db
  bool binary_is_temp_; ///< should we delete the binary fasta in destructor

  /**
   * Parses a database from the text based fasta file in the filename
   * member variable
   * reads in all proteins in the fasta file and creates a protein object
   * and adds them to the database protein array
   * total proteins in fasta file must not exceed MAX_PROTEIN constant
   * IF using light_protein functionality will not read in the sequence or id.
   * \returns true if success. false if failure.
   */
  bool parseTextFasta();

  /**
   * memory maps the binary fasta file for the database
   *\return true if successfully memory map binary fasta file, else false
   */
  bool memoryMap(
    int file_d  ///<  file descriptor -in
    );

  /**
   * Assumes that there is a 1 at the very end after all the proteins in binary file
   *\return true successfully populates the proteins from memory mapped binary fasta file, else false
   */
  bool populateProteinsFromMemmap();

  /**
   * \brief Parses a database from the binary fasta file in the filename
   * member variable.
   *
   * Memory maps the binary fasta file into memory. The protein
   * sequences are not copied, but just pointed to the memory mapped
   * location. 
   * \returns true if success. false if failure.
   */
  bool parseMemmapBinary();

  /**
   * intializes a database object
   */
  void init();

 public:
  /**
   * The suffix on binary and text fasta files.
   */
  static const std::string binary_suffix;
  static const std::string decoy_binary_suffix;
  static const std::string decoy_fasta_suffix;

  /**
   * \returns An (empty) database object.
   */

  Database();

  /**
   * \returns A new database object.
   */
  Database(
    const char* filename, ///< The file from which to parse the database. either text fasta file or binary fasta file -in
    bool is_memmap, ///< are we using a memory mapped binary fasta file, thus proteins are all memory mapped -in
    DECOY_TYPE_T decoys = NO_DECOYS ///< is this to be a decoy database
    );         

  void addProtein(
    PercolatorCrux::Protein* protein
  );

  /**
   * Frees an allocated protein object.
   */
  static void freeDatabase(
    Database* database ///< An allocated database -in
    );

  virtual ~Database();

  /**
   * Prints a database object to file.
   */
  void print(
    FILE* file    ///< output file stream -out             
    );

  /**
   * Parses a database from the file in the filename member variable
   * reads in all proteins in the fasta file and creates a protein object
   * and adds them to the database protein array
   * total proteins in fasta file must not exceed MAX_PROTEIN constant
   * \returns TRUE if success. FALSE if failure.
   */
  bool parse();

  /**
   * \brief Changes a database from one that reads from a fasta file to
   * one that reads from a binary/memmory mapped protein file.
   *
   * If database already has binary source (i.e. is_memmap == TRUE), 
   * returns TRUE.  
   * Opens the fasta file pointed to by filename for reading.  Creates a
   * file with the name given.  Reads in each protein from the text file
   * and serializes it to the output file.  Closes both files.  Changes
   * filename to point to new output file and sets is_memmap to true.
   * Parses the database.
   * \returns TRUE if all processes succeed, else FALSE.
   */
  bool transformTextToMemmap(
    const char* binary_protein_filename,
    bool binary_is_temp
    );

  /**
   * \returns FALSE if nth protein cannot be parsed or does not exist 
   */
  /**
  bool get_database_protein_at_idx(
      DATABASE_T* database, ///< A parsed database object -in
      unsigned int protein_idx, ///< The index of the protein to retrieve -in
      Protein** protein   ///< A pointer to a pointer to a PROTEIN object -out
      );
  **/

  /**
   * Using the fasta file the Database was instantiated with, write a
   * binary protein file in the given directory to use for memory
   * mapping.  If is_temp, delete the file on destruction.  Warns if
   * Database was not opened with a text file.
   */
  void createBinaryFasta(const char* directory, bool is_temp = false);

  /** 
   * Access routines of the form get_<object>_<field> and set_<object>_<field>. 
   */

  /**
   * Additional get and set methods
   */

  /**
   *\returns the filename of the database
   * returns a heap allocated new copy of the filename
   * user must free the return filename
   */
  char* getFilename();

  /**
   *\returns the pointer to the filename of the database
   * user must not free or change the filename
   */
  const char* getFilenamePointer();

  /**
   * sets the filename of the database
   * protein->sequence must been initiailized
   */
  void setFilename(
    const char* filename ///< the filename to add -in
    );

  /**
   *\returns TRUE|FALSE whether the database has been parsed?
   */
  bool getIsParsed();

  void setIsParsed(
    bool is_parsed
  );
  
  /**
   * \returns The type of shuffling used on the proteins in this database
   */
  DECOY_TYPE_T getDecoyType();

  /**
   *\returns the total number of proteins of the database
   */
  unsigned int getNumProteins();

  /**
   *\returns the src FILE* of the database
   */
  FILE* getFile();

  /**
   * sets the src FILE* of the database
   */
  void setFile(
    FILE* file ///< the src file to add -in
    );

  /**
   *\returns the nth protein of the database
   */
  PercolatorCrux::Protein* getProteinAtIdx(
    unsigned int protein_idx ///< The index of the protein to retrieve -in
    );

  /**
   *\returns the protein designated by protein id of the database
   */
  PercolatorCrux::Protein* getProteinByIdString(
    const char* protein_id ///< The id string for this protein -in
    );

  /**
   * sets the use_light_protein of the database
   */
  void setUseLightProtein(
    bool use ///< should I use the light/heavy functionality?
    );
  
  /**
   *\returns TRUE|FALSE whether the database uses light/heavy
   */
  bool getUseLightProtein();

  /**
   *sets TRUE,FALSE whether the database uses memory mapped
   */
  void setMemmap(
    bool is_memmap  ///< is the database memory mapped?
    );

  /**
   * increase the pointer_count produced by this database.
   * \returns database pointer
   */
  static Database* copyPtr(
    Database* database ///< the query database -in/out
  );

};

/**
 * Frees an allocated database_peptide_iterator object.
 */
void void_free_database_peptide_iterator(
  void* database_peptide_iterator ///< the iterator to free -in
  );

/**
 * The basic iterator functions.
 * \returns TRUE if there are additional peptides to iterate over, FALSE if not.
 */
bool void_database_peptide_iterator_has_next(
  void* database_peptide_iterator ///< the iterator of interest -in
  );

/**
 * \returns The next peptide in the database.
 */
PercolatorCrux::Peptide* void_database_peptide_iterator_next(
  void* database_peptide_iterator ///< the iterator of interest -in
  );

/***********************************
 * database sorted peptide iterator
 ***********************************/

/**********************************************************************
 * wrapper, for generate_peptides_iterator, cast back to original type
 ***********************************************************************/

/**
 * Frees an allocated database_sorted_peptide_iterator object.
 */
void void_free_database_sorted_peptide_iterator(
  void* database_peptide_iterator ///< the iterator to free -in
  );

/**
 * The basic iterator functions.
 * \returns TRUE if there are additional peptides to iterate over, FALSE if not.
 */
bool void_database_sorted_peptide_iterator_has_next(
  void* database_peptide_iterator ///< the iterator of interest -in
  );

/**
 * returns each peptide in sorted order
 * \returns The next peptide in the database.
 */
PercolatorCrux::Peptide* void_database_sorted_peptide_iterator_next(
  void* database_peptide_iterator ///< the iterator of interest -in
  );


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

}; // end namespace PercolatorCrux

#endif
