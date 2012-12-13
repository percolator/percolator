#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include "phos_loc_misc.h"
#include "phos_loc_modification.h"
#include "phos_loc_ion_type.h"
#include <string>
#include <vector>
#include <map>

namespace phos_loc {

class Parameters {
 public:
  Parameters();
  virtual ~Parameters();

  void Print();

  void AddVariableModification(const Modification& mod);
  void SetFragmentIonTypes();
  void InitAminoAcidMass();

  std::string parameter_file_path_;
  std::string fasta_db_name_;
  std::string fasta_file_path_;
  std::vector<Modification> fixed_mods_;
  std::vector<Modification> variable_mods_;
  std::map<unsigned short, int> mod_id_map_;
  Tolerance pep_tolerance_;
  ConsiderC13Precursor consider_c13_precursor_;
  Tolerance frag_tolerance_;
  PreprocessMethod preproc_method_;
  double preproc_parameters_[3];
  SpectraFormat spectra_format_;
  InstrumentType instrument_;
  ActivationType activation_type_;
  std::vector<IonType> specified_ion_types_;
  std::string spectra_path_;
  double amino_acid_mass_[LENGTH_AMINO_ACID_TABLE][NUM_MASS_TYPES];
};

} // namespace phos_loc

#endif // PARAMETERS_H_
