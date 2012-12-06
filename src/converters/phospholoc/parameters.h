#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include "misc.h"
#include "modification.h"
#include "ion_type.h"
#include <string>
#include <vector>

class Parameters {
 public:
  Parameters();
  virtual ~Parameters();

  void Print();

 private:
  void InitAminoAcidMass();

 public:
  std::string parameter_file_path_;
  std::string fasta_db_name_;
  std::string fasta_file_path_;
  std::vector<Modification> fixed_mods_;
  std::vector<Modification> variable_mods_;
  Tolerance pep_tolerance_;
  ConsiderC13Precursor consider_c13_precursor_;
  Tolerance frag_tolerance_;
  PreprocessMethod preproc_method_;
  double preproc_parameters_[3];
  SpectraFormat spectra_format_;
  InstrumentType instrument_;
  std::vector<IonType> specified_ion_types_;
  std::string spectra_path_;
  double amino_acid_mass_[LENGTH_AMINO_ACID_TABLE][NUM_MASS_TYPES];
};

#endif // PARAMETERS_H_
