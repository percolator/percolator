/********************************************
 *** adapted from crux/src/c/Database.cpp ***
 ********************************************/
 

/*************************************************************************//**
 * \file Database.cpp
 * \brief Object for representing a database of protein sequences.
 ****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>
#include <sys/types.h>
#include <fcntl.h>
#ifndef _MSC_VER
#include <sys/mman.h>
#include <unistd.h>
#include <errno.h>
#endif
#include "Database.h"

#include <map>
#include <vector>
#include <iostream>

#ifdef _MSC_VER
/*********************************************************
 This function replaces the GNU extension of the same name.
 Reads a line from the given stream.
 *********************************************************/
int getline(char **lineptr, size_t *n, FILE *stream) {

  const size_t BUFFSIZE = 100;

  // Check the input values.
  if (lineptr == NULL || stream == NULL) {
    errno = EINVAL;
    return -1;
  }

  // Read the first character from the stream.
  size_t index = 0;
  int c = fgetc(stream);
  int e = ferror(stream);
  if (c == EOF || e) {
    return -1;
  }

  // Allocate an buffer if needed.
  if (*lineptr == NULL) {
    *lineptr = (char *) malloc((*n + BUFFSIZE) * sizeof(char));
    *n += BUFFSIZE;
  }

  // Copy from the stream to the buffer until we find a line end.
  while(c != '\n' && c != EOF && !e) {
    (*lineptr)[index++] = c;
    if (index > *n - 1) {
      // Out of space, expand the buffer.
      *lineptr = (char *) realloc(*lineptr, *n + BUFFSIZE);
      *n += BUFFSIZE;
    }
    c = fgetc(stream);
    e = ferror(stream);
  }

  // Reached end of line, end of file, or read error.
  if (!e) {

    if (c != EOF) {
      (*lineptr)[index++] = c;
      if (index > (*n - 1)) {
        *lineptr = (char *) realloc(*lineptr, *n + 1);
        (*n)++;
      }
    }
    
    // Terminate the string.
    (*lineptr)[index] = 0;
    
    // Return the length of the string
    // without the terminating null.
    return index;
  }
  else {
    // Some sort of read error
    return -1;
  }
}
#endif

using namespace std;
using namespace PercolatorCrux;
const string Database::binary_suffix = "-binary-fasta";
const string Database::decoy_binary_suffix = "-binary-fasta-decoy";
const string Database::decoy_fasta_suffix = "-random.fasta";

/**
 * intializes a database object
 */
void Database::init(){
  file_ = NULL;
  is_parsed_ = false;
  size_ = 0; 
  use_light_protein_ = false; 
  is_memmap_ = false;
  data_address_ = NULL;
  pointer_count_ = 1;
  file_size_ = 0;
  is_hashed_ = false;
  proteins_ = new vector<Protein*>();
  protein_map_ = new map<char*, Protein*, cmp_str>();
  decoys_ = NO_DECOYS;
  binary_is_temp_ = false;
}

/**
 * \returns An (empty) database object.
 */
Database::Database() {
  init();
}


/**
 * \returns A new database object.
 */
Database::Database(
  const char*         filename, ///< The file from which to parse the database. 
  ///< either text fasta file or binary fasta file -in
  bool is_memmap, ///< are we using a memory mapped binary fasta file? 
  ///< If so, all proteins are memory mapped -in
  DECOY_TYPE_T decoys ///< is this a decoy database
  )         
{
  //carp(CARP_DEBUG, "Creating new database from '%s'", filename);
  init();
  is_memmap_ = is_memmap;
  if( is_memmap_ ){
    binary_filename_ = filename;
  } else {
    fasta_filename_ = filename;
  }
  decoys_ = decoys;
}  

/**
 * Frees an allocated protein object.
 */
void Database::freeDatabase(
  Database* database ///< An allocated database -in
  )
{
  if(database == NULL){
    return;
  }

  // decrement database pointer counter
  --database->pointer_count_;
  //carp(CARP_DETAILED_DEBUG, "Database pointer count %i", database->pointer_count_);

  // DEBUG show the databse pointer count
  // printf("Free: After free: %s: %d\n", database->pointer_count);

  // only free up memory when remaining pointers are from the proteins
  if((size_t)database->pointer_count_ > database->proteins_->size() ){ 
    return;
  }

  delete database;
}

Database::~Database() {

  // only free proteins if been parsed and file has been opened
  if(is_parsed_){
    //carp(CARP_DEBUG, "Freeing database.");
    
    // free each protein in the array
    unsigned int protein_idx;
    for(protein_idx=0; protein_idx < proteins_->size(); ++protein_idx){
      delete ((*proteins_)[protein_idx]);
    }
    delete proteins_;
    delete protein_map_; // contents already deleted
    
    if (file_ != NULL) {
      // close file handle
      //carp(CARP_DEBUG, "Closing database filehandle");
      fclose(file_);
    }
  }
}

/**
 * Prints a database object to file.
 */
void Database::print(
  FILE* file    ///< output file stream -out             
  )
{
  Protein* protein = NULL;

  fprintf(file, "filename:%s\n", fasta_filename_.c_str());
}

void Database::addProtein(
  Protein* protein
  ) {

  protein->setDatabase(this);
      
  // add protein to database
  proteins_->push_back(protein);
  
  protein->setProteinIdx(proteins_->size()-1);

  if (is_hashed_) {
    char* id = protein->getIdPointer();
    protein_map_->insert(make_pair(id, protein));
  }
}

/**
 * Parses a database from the text based fasta file in the filename
 * member variable
 * reads in all proteins in the fasta file and creates a protein object
 * and adds them to the database protein array
 * total proteins in fasta file must not exceed MAX_PROTEIN constant
 * IF using light_protein functionality will not read in the sequence or id.
 * \returns true if success. false if failure.
 */
bool Database::parseTextFasta()
{
  unsigned long working_index;
  FILE* file = NULL;
  char* new_line = NULL;
  int line_length;
  size_t buf_length = 0;
  Protein* new_protein;
  unsigned int protein_idx;

  //carp(CARP_DEBUG, "Parsing text fasta file '%s'", fasta_filename_.c_str());
  // check if already parsed
  if(is_parsed_){
    return true;
  }
  
  // open file and 
  file = fopen(fasta_filename_.c_str(), "rb");
  
  // check if succesfully opened file
  if(file == NULL){
    //carp(CARP_ERROR, "Failed to open fasta file %s", fasta_filename_.c_str());
    return false;
  }
   
  working_index = ftell(file);
  // check each line until reach '>' line
  while((line_length =  getline(&new_line, &buf_length, file)) != -1){
    if(new_line[0] == '>'){
      // the new protein to be added
      new_protein = new Protein();
      
      // do not parse the protein sequence if using light/heavy functionality
      if(use_light_protein_){
        // set light and offset
        new_protein->setOffset(working_index);
        new_protein->setIsLight(true);
      }
      else{
        // rewind to the beginning of the protein to include ">" line
        fseek(file, working_index, SEEK_SET);
        
        // failed to parse the protein from fasta file
        // protein offset is set in the parse_protein_fasta_file method
        if(!new_protein->parseProteinFastaFile(file)){
          fclose(file);
          delete new_protein;
          for(protein_idx=0;protein_idx<proteins_->size();protein_idx++){
            delete (proteins_->at(protein_idx));
          }
          proteins_->clear();
          //carp(CARP_ERROR, "Failed to parse fasta file");
          return false;
        }
        new_protein->setIsLight(false);
      }
      // add protein to database
      proteins_->push_back(new_protein);
      // set protein index, database
      new_protein->setProteinIdx(proteins_->size()-1);
      new_protein->setDatabase(this);
    }
    working_index = ftell(file);
  }
  free(new_line);
  
  // yes the database is paresed now..!!
  is_parsed_ = true;
  file_ = file;
  return true;
}


/**
 * Parses a database from the file in the filename member variable
 * The is_memmap field in the database struct determines whether the
 * input file is a binary fasta file or normal text fasta file.
 *
 * Uses the traditional text fasta file which
 * it parses out the various peptides for each protein. Only when
 * using text fasta file can you use light/heavy protein, in which 
 * if using light_protein functionality will not read in the sequence
 * or id. Will parse sequence if protein  
 * is needed, lazy parsing.
 *
 * Reads in all proteins in file and creates a protein object
 * and adds them to the database protein array
 * total proteins in fasta file must not exceed MAX_PROTEIN constant
 *
 * \returns true if success. false if failure.
 */
bool Database::parse()
{
  return parseTextFasta();
}

/**
 *\returns the pointer to the filename of the database
 * user must not free or change the filename
 */
const char* Database::getFilenamePointer()
{
  return fasta_filename_.c_str();
}

/**
 *\returns true|false whether the database has been parsed?
 */
bool Database::getIsParsed()
{
  return is_parsed_;
}

void Database::setIsParsed(
  bool is_parsed
  ) {
  is_parsed_ = is_parsed;
}

/**
 * \returns The type of shuffling used on the proteins in this database
 */
DECOY_TYPE_T Database::getDecoyType(){
  return decoys_;
}

/**
 *\returns the total number of proteins of the database
 */
unsigned int Database::getNumProteins()
{
  return proteins_->size();
}

/**
 *\returns the src FILE* of the database
 */
FILE* Database::getFile()
{
  return file_;
}

/**
 * sets the src FILE* of the database
 */
void Database::setFile(
  FILE* file ///< the src file to add -in
  )
{
  file_ = file;
}

/**
 * \returns the nth protein of the database
 * 
 */
Protein* Database::getProteinAtIdx(
  unsigned int protein_idx ///< The index of the protein to retrieve -in
  )
{
  ////carp(CARP_DETAILED_DEBUG, "Getting db protein idx = %i, num proteins %i", 
  //     protein_idx, database->proteins.size());
  if( protein_idx >= proteins_->size()){
    //carp(CARP_FATAL, "Protein index %i out of bounds.  %i proteins in the database", protein_idx, proteins_->size());
  }

  return proteins_->at(protein_idx);
}

/**
 *\returns the protein designated by protein id of the database
 */
Protein* Database::getProteinByIdString(
  const char* protein_id ///< The id string for this protein -in
  ) {

  //TODO - Implement as a hashtable rather than a map to make 
  //this even faster if needed.
  Protein* protein = NULL;
  if (is_hashed_) {
    map<char*, Protein*, cmp_str>::iterator find_iter;
    find_iter = protein_map_->find((char*)protein_id);

    if (find_iter != protein_map_->end()) {
      protein = find_iter->second;
    }
  } else {
    //create the hashtable of protein ids
    for (unsigned int protein_idx = 0;
      protein_idx < proteins_->size();
      protein_idx++) {

      Protein* current_protein = proteins_->at(protein_idx);
      char* current_id = current_protein->getIdPointer();
      protein_map_->insert(make_pair(current_id, current_protein));

      if (strcmp(current_id, protein_id)==0) {
        protein = current_protein;
      }
        
    }
    is_hashed_ = true;
  }
  return protein;
}

/**
 * increase the pointer_count produced by this database.
 * \returns database pointer
 */
Database* Database::copyPtr(
  Database* database ///< the query database -in/out
  )
{
  if( database == NULL ){
    return NULL;
  }
  ++database->pointer_count_;
  return database;
}


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

