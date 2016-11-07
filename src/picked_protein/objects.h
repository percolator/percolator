/*****************************************
 *** adapted from crux/src/c/objects.h ***
 *****************************************/


/**
 * \file objects.h 
 * $Revision: 1.62 $
 * \brief The defined objects
 *****************************************************************************/
#ifndef OBJECTS_H 
#define OBJECTS_H

#include <vector>
#include <cstdlib>

// Macro allowing us to build using floats or double
#ifdef USE_DOUBLES
typedef double FLOAT_T;
#else
typedef float FLOAT_T;
#endif

/** 
  * MT: temporary namespace to allow automatic Crux-Percolator build
  */
namespace PercolatorCrux {
  
  /**
   * \class Protein
   * \brief A protein sequence
   */
  class Protein;
  
  /**
   * \class Peptide
   * \brief A peptide subsequence of a protein
   */
  class Peptide;
  
  /**
   * \class PeptideSrc
   * \brief object for mapping a peptide to it's parent protein.
   */
  class PeptideSrc;

  /**
   *\typedef  PeptideSrcIterator 
   *\brief An object to iterate over the PeptideSrc in a peptide  
   */
  typedef std::vector<PeptideSrc*>::iterator  PeptideSrcIterator; 
  
  /**
   * \struct residue_iterator
   * \brief Object to iterate over the residues in a peptide, starting at the
   * first residue of the peptide, and proceeding in order.
   */
  struct residue_iterator {
    Peptide*  peptide; ///< The peptide whose residues to iterate over.
    char*   sequence;    ///< The peptide sequence
    int     residue_idx; ///< The index of the current peak
  };
  
  /**
   * \typedef RESIDUE_ITERATOR_T 
   * \brief An object to iterate over the residues in a peptide
   */
  typedef struct residue_iterator RESIDUE_ITERATOR_T;

  /**
   * \class ProteinPeptideIterator
   * \brief An object to iterate over the peptides in a protein sequence
   */
  class ProteinPeptideIterator;

  /**
   * \class Database
   * \brief A database of protein sequences.
   */
  class Database;

  /**
   * \enum DECOY_TYPE_T
   */
  enum DECOY_TYPE_T {
    INVALID_DECOY_TYPE,
    NO_DECOYS,
    PROTEIN_REVERSE_DECOYS,
    PROTEIN_SHUFFLE_DECOYS,
    PEPTIDE_SHUFFLE_DECOYS,
    PEPTIDE_REVERSE_DECOYS,
    NUMBER_DECOY_TYPES
  };


  /**
   * \enum _digest_type
   * The rule governing how a peptide was cleaved from its source
   * protein sequence. 
   */
  enum _digest_type {
    INVALID_DIGEST,      ///< required invalid value for the enum
    FULL_DIGEST,         ///< c- AND n-term specific to ENZYME_T
    PARTIAL_DIGEST,      ///< c- OR n-term specific to ENZYME_T
    NON_SPECIFIC_DIGEST, ///< not specific to any enzyme cleavage rules
    NUMBER_DIGEST_TYPES  ///< keep last, number of types
  };

  /**
   * \typedef DIGEST_T
   * \brief The rule governing how a peptide was digested.  Used in
   * conjunction with ENZYME_T to define how peptides are generated.
   */
  typedef enum _digest_type DIGEST_T;

  /**
   * \enum _enzyme_type
   */
  enum _enzyme_type {
    INVALID_ENZYME,        ///< required invalid value for the enum
    NO_ENZYME,             ///< cleave anywhere
    TRYPSIN,               ///< cleave after K or R, not before P
    TRYPSINP,               ///< cleave after K or R
    CHYMOTRYPSIN,          ///< cleave after FWYL, not before P
    ELASTASE,              ///< cleave after ALIV, not before P
    CLOSTRIPAIN,           ///< cleave after R
    CYANOGEN_BROMIDE,      ///< cleave after M
    IODOSOBENZOATE,        ///< cleave after W
    PROLINE_ENDOPEPTIDASE, ///< cleave after P
    STAPH_PROTEASE,        ///< cleave after E
    ASPN,                  ///< cleave before D
    LYSC,                  ///< cleave after K , not befor P 
    LYSN,                  ///< cleave before K 
    ARGC,                  ///< cleave after R, not before P
    GLUC,                  ///< cleave after D or E, not before P
    PEPSINA,               ///< cleave after FL, not before P 
    ELASTASE_TRYPSIN_CHYMOTRYPSIN, ///< cleave after ALIVKRWFY, not before P
    CUSTOM_ENZYME,         ///< cleave after/before user-defined residues
    NUMBER_ENZYME_TYPES    ///< leave last, number of types
  };

  /**
   * \typedef ENZYME_T
   * \brief The enzyme with which a peptide was digested.  Used in
   * conjunction with DIGEST_T to define how peptides are generated.
   */
  typedef enum _enzyme_type ENZYME_T;
};



#endif
