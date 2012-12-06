#include "modification.h"

Modification::Modification() {
  Init("Carbamidomethyl", 57.0513, 57.021464, 0, 0, "C", ANYWHERE, 0);
}

Modification::Modification(const std::string name, char max_num) {
  if (!name.compare("Carbamidomethyl"))
    Init(name, 57.0513, 57.021464, 0, 0, "C", ANYWHERE, max_num);
  else if (!name.compare("Oxidation"))
    Init(name, 15.9994, 15.994915, 0, 0, "M", ANYWHERE, max_num);
  else if (!name.compare("Phospho"))
    Init(name, 79.9799, 79.966331, 97.9952, 97.976896, "STY", ANYWHERE, max_num);
  else if (!name.compare("PhosphoST"))
    Init(name, 79.9799, 79.966331, 97.9952, 97.976896, "ST", ANYWHERE, max_num);
  else if (!name.compare("PhosphoY"))
    Init(name, 79.9799, 79.966331, 0, 0, "Y", ANYWHERE, max_num);
  else if (!name.compare("Acetyl"))
    Init(name, 42.0367, 42.010565, 0, 0, "*", PROT_NTERM, max_num);
  else
    fprintf(stderr, "- Please specify other fields of this modification: %s.\n",
            name.c_str());
}


Modification::Modification(const std::string name, double* mass_shift,
                           double* neutral_loss, const std::string mod_aas,
                           ModificationPosition position, char max_num) {
  Init(name, mass_shift[AVERAGE], mass_shift[MONO],
       neutral_loss[AVERAGE], neutral_loss[MONO], mod_aas, position, max_num);
}

Modification::Modification(const Modification& modification)
    : name_(modification.name_),
      position_(modification.position_),
      max_num_per_peptide_(modification.max_num_per_peptide_) {
  mass_shift_[AVERAGE] = modification.mass_shift_[AVERAGE];
  mass_shift_[MONO] = modification.mass_shift_[MONO];
  neutral_loss_[AVERAGE] = modification.neutral_loss_[AVERAGE];
  neutral_loss_[MONO] = modification.neutral_loss_[MONO];
  for (int i = 0; i < LENGTH_AMINO_ACID_TABLE; ++i)
    modifying_[i] = modification.modifying_[i];
}

Modification::~Modification() {
}

void Modification::Init(const std::string name, double mass_shift_aver,
                        double mass_shift_mono, double neutral_loss_aver,
                        double neutral_loss_mono, const std::string mod_aas,
                        ModificationPosition position, char max_num) {
  name_ = name;
  mass_shift_[AVERAGE] = mass_shift_aver;
  mass_shift_[MONO] = mass_shift_mono;
  neutral_loss_[AVERAGE] = neutral_loss_aver;
  neutral_loss_[MONO] = neutral_loss_mono;
  MapModAAs(mod_aas);
  position_ = position;
  max_num_per_peptide_ = max_num;
}

void Modification::MapModAAs(const std::string mod_aas) {
  memset(modifying_, 0, LENGTH_AMINO_ACID_TABLE * sizeof(bool));
  if (!mod_aas.compare("*")) {
    for (int i = 0; i < LENGTH_AMINO_ACID_TABLE; ++i)
      modifying_[i] = true;
    return;
  }
  for (std::string::const_iterator it = mod_aas.begin(); it != mod_aas.end(); ++it) {
    modifying_[*it-'A'] = true;
  }
}
