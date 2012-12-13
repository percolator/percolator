#ifndef MISC_H_
#define MISC_H_

#include <limits.h> // for the bounds of types
#include <stddef.h> // for size_t
#include <string.h> // for memcpy
#include <math.h>
#include <algorithm>

// A macro to disallow the copy constructor and operator= functions
// This should be used in the private: declarations for a class
#define DISALLOW_COPY_AND_ASSIGN(TypeName) \
  TypeName(const TypeName&);               \
  void operator=(const TypeName&)

namespace phos_loc {

const int LENGTH_AMINO_ACID_TABLE = 26;
const double PROTON = 1.007276466812; // 1.00727647012
const double ELECTRON = 0.00054857990946; // 0.00054857990313
const double NEUTRON = 1.00866491600; // 1.00866490414
const double C13ISOTOPE = 1.0033548378; // 1.003354826
const double HYDROGEN[2] = {1.00794, 1.00782503207}; // 1.007825035
const double CARBON[2] = {12.0107, 12.0000000};
const double NITROGEN[2] = {14.0067, 14.0030740048}; // 14.003074002
const double OXYGEN[2] = {15.9994, 15.99491461956}; // 15.99491463
const double H2O[2] = {18.01528, 18.0105646837}; // 18.0105647
const double NH3[2] = {17.03052, 17.02654910101}; // 17.026549107
const double PHOSPHO_LOSS[2] = {97.9952, 97.976896};
const double NOT_VALID_MASS = -32768.0;
const unsigned short UNIMOD_PHOSPHO_ID = 21;

// maximum number of variable modification types that can be considered
const int MAX_NUM_VAR_MOD_TYPES = 20;
// maximum number of variable modifications per peptide
const int MAX_NUM_VAR_MODS_PER_PEPTIDE = 4;

enum MassType {
  AVERAGE = 0,
  MONO,
  NUM_MASS_TYPES
};

// TODO: define other properties here
const double AMINO_ACID_MASS[LENGTH_AMINO_ACID_TABLE][NUM_MASS_TYPES] = {
  { 71.0788,  71.037113805},  // A, Ala, Alanine
  {114.5962, 114.534935268},  // B, Asx, Avg_Asn_Asp
  {103.1388, 103.009184505},  // C, Cys, Cysteine
  {115.0886, 115.026943065},  // D, Asp, Aspartic Acid
  {129.1155, 129.042593135},  // E, Glu, Glutamic Acid
  {147.1766, 147.068413945},  // F, Phe, Phenylalanine
  { 57.0519,  57.021463735},  // G, Gly, Glycine
  {137.1411, 137.058911875},  // H, His, Histidine
  {113.1594, 113.084064015},  // I, Ile, Isoleucine
  {113.1594, 113.084064015},  // J, Leu_or_Ile
  {128.1741, 128.094963050},  // K, Lys, Lysine
  {113.1594, 113.084064015},  // L, Leu, Leucine
  {131.1926, 131.040484645},  // M, Met, Methionine
  {114.1038, 114.042927470},  // N, Asn, Asparagine
  {114.1472, 114.079312980},  // O, Ornithine
  { 97.1167,  97.052763875},  // P, Pro, Proline
  {128.1307, 128.058577540},  // Q, Gln, Glutamine
  {156.1875, 156.101111050},  // R, Arg, Arginine
  { 87.0782,  87.032028435},  // S, Ser, Serine
  {101.1051, 101.047678505},  // T, Thr, Threonine
  {150.0388, 150.953633405},  // U, SeC, Selenocysteine
  { 99.1326,  99.068413945},  // V, Val, Valine
  {186.2132, 186.079312980},  // W, Trp, Tryptophan
  {111.2137, 111.213700000},  // X, Xla, Average AA
  {163.1760, 163.063328575},  // Y, Tyr, Tyrosine
  {128.6231, 128.550585338}   // Z, Glx, Avg_Gln_Glu
};

enum IonCategory {
  NTERM_IONS = 0,
  CTERM_IONS,
  INTERNAL_IONS,
  IMMONIUM_IONS,
  SATELLITE_IONS,
  PRECURSOR_IONS,
  NUM_ION_CATEGORIES
};

enum IonTypeCode {
  IM_ITC = 0, // immonium ions
  A_ITC,
  A_NH3_ITC,  // a-NH3
  A_H2O_ITC,
  B_ITC,
  B_NH3_ITC,
  B_H2O_ITC,
  B_NL1_ITC,  // for neutral loss of PTM
  B_NL2_ITC,
  C_ITC,
  C_H_ITC,    // c-H
  X_ITC,
  Y_ITC,
  Y_NH3_ITC,
  Y_H2O_ITC,
  Y_NL1_ITC,  // for neutral loss of PTM
  Y_NL2_ITC,
  Z_ITC,
  Z__H_ITC,   // z+H
  Z__2H_ITC,  // z+2H
  IYA_ITC,   // internal ya
  IYB_ITC,
  NUM_ION_TYPE_CODES
};

enum ModificationPosition {
  ANYWHERE = 0,
  PEP_NTERM,
  PEP_CTERM,
  PROT_NTERM,
  PROT_CTERM,
  NUM_MODIFICATION_POSITIONS
};

enum DigestType {
  FULLY_SPECIFIC = 0,
  SEMI_SPECIFIC,
  NON_SPECIFIC,
  NUM_DIGEST_TYPES
};

enum ToleranceUnit {
  PPM_TOL = 0,    // fraction expressed as parts per million, abs_tol=MZ*tol*0.000001
  DA_TOL,         // absolute units of Da, abs_tol=tol(Th)
  MMU_TOL,        // absolute milli-mass units, abs_tol=tol*0.001 Da
  PERCENTAGE_TOL, // fraction expressed as a percentage, abs_tol=MZ*tol*0.01
  NUM_TOLERANCE_UNITS
};

struct Tolerance {
  double value;
  double drift;  // instrument drift, default=0.0
  ToleranceUnit unit;
  MassType mass_type;

  Tolerance()
      : value(0.5),
        drift(0),
        unit(DA_TOL),
        mass_type(MONO) {
  }
  Tolerance(double val, double dr, ToleranceUnit un, MassType mt)
      : value(val),
        drift(dr),
        unit(un),
        mass_type(mt) {
  }
  Tolerance(const Tolerance& tol)
      : value(tol.value),
        drift(tol.drift),
        unit(tol.unit),
        mass_type(tol.mass_type) {
  }
};

struct Precursor {
  char charge;
  double m_over_z;
  MassType mass_type;
  double intensity;

  Precursor()
      : charge(0),
        m_over_z(0),
        mass_type(MONO),
        intensity(0) {
  }
  Precursor(char chg, double mz, MassType mt, double inten = 0)
      : charge(chg),
        m_over_z(mz),
        mass_type(mt),
        intensity(inten) {
  }
  Precursor(const Precursor& prec)
      : charge(prec.charge),
        m_over_z(prec.m_over_z),
        mass_type(prec.mass_type),
        intensity(prec.intensity) {
  }
};

enum ActivationType {
  CAD_CID = 0,
  ECD_ETD,
  HCD,
  NUM_ACTIVATION_TYPES
};

enum ConsiderC13Precursor {  // for high accuracy data
  MONO_P0 = 0, // a) TOL > absolute(exp - calc)
  MONO_P1,     // a) and b) TOL > absolute(exp - calc - C13ISOTOPE)
  MONO_P2      // a) and b) and c) TOL > absolute(exp - calc - 2*C13ISOTOPE)
};

enum PreprocessMethod {
  TOPN = 0,     // retain the N most intense peaks
  THRESHOLD,    // keep all the peaks with intensity > THRESHOLD
  TOPN_WINDOW,  // top N intense peaks per window of width W(Th)
  NUM_PREPROCESS_METHODS
};

enum SpectraFormat {
  MGF = 0,
  MS2,
  DTA,
  MZML,
  MZXML,
  AGILENT_RAW,
  BRUKER_RAW,
  THERMO_RAW,
  WATERS_RAW,
  NUM_SPECTRA_FORMATS
};

const int MAX_LENGTH_SPECTRUM_LINE = 256;

enum InstrumentType {
  ESI_TRAP = 0,
  ESI_QTOF,
  MALDI_QTOF,
  MALDI_TOF_TOF,
  ETD_TRAP,
  NUM_INSTRUMENT_TYPES
};

// used to specify variable modifications in a peptide
struct LocationMod {
  unsigned char aa_idx; // index of the modified amino acid in a peptide
  unsigned short mod_id; // unimod id
  LocationMod(unsigned char aaidx, unsigned short modid)
      : aa_idx(aaidx),
        mod_id(modid) {
  }
};

} // namespace phos_loc

#endif // MISC_H_
