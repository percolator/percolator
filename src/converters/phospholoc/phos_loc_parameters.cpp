#include "phos_loc_parameters.h"

namespace phos_loc {

Parameters::Parameters()
    : pep_tolerance_(3, 0, DA_TOL, MONO),
      consider_c13_precursor_(MONO_P0),
      frag_tolerance_(0.5, 0, DA_TOL, MONO),
      preproc_method_(TOPN_WINDOW),
      spectra_format_(DTA),
      instrument_(ESI_TRAP),
      activation_type_(HCD) {
  AddVariableModification(Modification("Phospho"));
  preproc_parameters_[0] = 50; // starting m/z, -1 means the lowest m/z; 50.0
  preproc_parameters_[1] = 100.0; // window width
  preproc_parameters_[2] = 10; // top-N per window
  SetFragmentIonTypes();
  InitAminoAcidMass();
}

Parameters::~Parameters() {
}

void Parameters::Print() {
  fprintf(stderr, ">Fixed modifications:\n");
  for (std::vector<Modification>::const_iterator it = fixed_mods_.begin();
       it != fixed_mods_.end(); ++it)
    fprintf(stderr, "%s %f\n", it->name().c_str(), it->mass_shift(MONO));

  fprintf(stderr, ">Variable modifications:\n");
  for (std::vector<Modification>::const_iterator it = variable_mods_.begin();
       it != variable_mods_.end(); ++it)
    fprintf(stderr, "%s %f\n", it->name().c_str(), it->mass_shift(MONO));

  fprintf(stderr, ">Selected ion types:\n");
  for (std::vector<IonType>::const_iterator it = specified_ion_types_.begin();
       it != specified_ion_types_.end(); ++it)
    fprintf(stderr, " %s", it->ion_name().c_str());
  fprintf(stderr, "\n");
}

void Parameters::AddVariableModification(const Modification& mod) {
  std::pair<std::map<unsigned short, int>::iterator, bool> res;
  res = mod_id_map_.insert(std::pair<unsigned short, int>(
          mod.unimod_id(), variable_mods_.size()));
  if (res.second)
    variable_mods_.push_back(mod);
}

void Parameters::SetFragmentIonTypes() {
  switch (activation_type_) {
    case CAD_CID:
      specified_ion_types_.push_back(std::string("b"));
      specified_ion_types_.push_back(std::string("y"));
      break;
    case ECD_ETD:
      specified_ion_types_.push_back(std::string("c"));
      specified_ion_types_.push_back(std::string("z"));
      specified_ion_types_.push_back(std::string("z+h"));
      break;
    case HCD:
      specified_ion_types_.push_back(std::string("b"));
      specified_ion_types_.push_back(std::string("b-phospho"));
      specified_ion_types_.push_back(std::string("y"));
      specified_ion_types_.push_back(std::string("y-phospho"));
      break;
    default:
      specified_ion_types_.push_back(std::string("b"));
      specified_ion_types_.push_back(std::string("y"));
      break;
  }
}

void Parameters::InitAminoAcidMass() {
  for (int i = 0; i < LENGTH_AMINO_ACID_TABLE; ++i) {
    amino_acid_mass_[i][AVERAGE] = AMINO_ACID_MASS[i][AVERAGE];
    amino_acid_mass_[i][MONO] = AMINO_ACID_MASS[i][MONO];
    for (std::vector<Modification>::size_type j = 0; j < fixed_mods_.size(); ++j) {
      if (fixed_mods_[j].modifying(i) && fixed_mods_[j].position() == ANYWHERE) {
        amino_acid_mass_[i][AVERAGE] += fixed_mods_[j].mass_shift(AVERAGE);
        amino_acid_mass_[i][MONO] += fixed_mods_[j].mass_shift(MONO);
      }
    }
  }
}

} // namespace phos_loc
