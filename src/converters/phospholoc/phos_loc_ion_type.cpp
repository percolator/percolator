#include "phos_loc_ion_type.h"
#include "phos_loc_peptide.h"

namespace phos_loc {

IonType::IonType()
    : ion_category_(NTERM_IONS),
      mass_type_(MONO),
      mass_shift_(0),
      charge_(0),
      condition_func_ptr_(NULL) {
}

IonType::IonType(const IonType& ion_type)
    : ion_name_(ion_type.ion_name_),
      ion_category_(ion_type.ion_category_),
      mass_type_(ion_type.mass_type_),
      mass_shift_(ion_type.mass_shift_),
      charge_(ion_type.charge_),
      condition_func_ptr_(ion_type.condition_func_ptr_) {
}

IonType::IonType(const std::string& ion_name, MassType mass_type)
    : mass_type_(mass_type) {
  InterpretIonName(ion_name);
}

IonType::~IonType() {
}

IonType& IonType::operator=(const IonType& ion_type) {
  if (this != &ion_type) {
    ion_name_ = ion_type.ion_name_;
    ion_category_ = ion_type.ion_category_;
    mass_type_ = ion_type.mass_type_;
    mass_shift_ = ion_type.mass_shift_;
    charge_ = ion_type.charge_;
    condition_func_ptr_ = ion_type.condition_func_ptr_;
  }
  return *this;
}

bool IonType::operator<(const IonType& ion_type) const {
  if (ion_name_ != ion_type.ion_name_)
    return ion_name_ < ion_type.ion_name_;
  else
    return charge_ < ion_type.charge_;
}

double IonType::GetMoverZ(const Peptide& peptide, int cleavage_site) const {
  double mass_neutral_residues = 0.0;
  if (ion_category_ == NTERM_IONS)
    mass_neutral_residues = peptide.nterm_ladder(cleavage_site);
  else if (ion_category_ == CTERM_IONS)
    mass_neutral_residues = peptide.cterm_ladder(cleavage_site);
  if (condition_func_ptr_ != NULL &&
      !((this->*condition_func_ptr_)(peptide, cleavage_site)))
    return NOT_VALID_MASS;
  return (mass_neutral_residues + mass_shift_ + PROTON * charge_) / charge_;
}

void IonType::Set(const std::string& ion_name, IonCategory ion_category,
                  double mass_shift, int charge) {
  ion_name_ = ion_name;
  ion_category_ = ion_category;
  mass_shift_ = mass_shift;
  charge_ = charge;
  condition_func_ptr_ = NULL;
}

void IonType::SetCondition() {
  if (!ion_name_.compare("b-phospho") ||
      !ion_name_.compare("y-phospho"))
    condition_func_ptr_ = &IonType::IsContainingPhospho;
}

void IonType::InterpretIonName(const std::string ion_name) {
  if (!ion_name.compare("b"))
    Set(ion_name, NTERM_IONS, 0, 1);
  else if (!ion_name.compare("y"))
    Set(ion_name, CTERM_IONS, H2O[mass_type_], 1);
  else if (!ion_name.compare("c"))
    Set(ion_name, NTERM_IONS, NH3[mass_type_], 1);
  else if (!ion_name.compare("c-h"))
    Set(ion_name, NTERM_IONS, NH3[mass_type_]-HYDROGEN[mass_type_], 1);
  else if (!ion_name.compare("z"))
    Set(ion_name, CTERM_IONS,
        -HYDROGEN[mass_type_]-NITROGEN[mass_type_]+OXYGEN[mass_type_], 1);
  else if (!ion_name.compare("z+h"))
    Set(ion_name, CTERM_IONS, -NITROGEN[mass_type_]+OXYGEN[mass_type_], 1);
  else if (!ion_name.compare("z+2h"))
    Set(ion_name, CTERM_IONS,
        HYDROGEN[mass_type_]-NITROGEN[mass_type_]+OXYGEN[mass_type_], 1);
  else if (!ion_name.compare("b-nh3"))
    Set(ion_name, NTERM_IONS, -NH3[mass_type_], 1);
  else if (!ion_name.compare("b-h2o"))
    Set(ion_name, NTERM_IONS, -H2O[mass_type_], 1);
  else if (!ion_name.compare("b-phospho")) {
    Set(ion_name, NTERM_IONS, -PHOSPHO_LOSS[mass_type_], 1);
    SetCondition();
  }
  else if (!ion_name.compare("y-nh3"))
    Set(ion_name, CTERM_IONS, H2O[mass_type_]-NH3[mass_type_], 1);
  else if (!ion_name.compare("y-h2o"))
    Set(ion_name, CTERM_IONS, 0, 1);
  else if (!ion_name.compare("y-phospho")) {
    Set(ion_name, CTERM_IONS, H2O[mass_type_]-PHOSPHO_LOSS[mass_type_], 1);
    SetCondition();
  }
  else if (!ion_name.compare("a"))
    Set(ion_name, NTERM_IONS, -CARBON[mass_type_]-OXYGEN[mass_type_], 1);
  else if (!ion_name.compare("a-nh3"))
    Set(ion_name, NTERM_IONS,
        -CARBON[mass_type_]-OXYGEN[mass_type_]-NH3[mass_type_], 1);
  else if (!ion_name.compare("a-h2o"))
    Set(ion_name, NTERM_IONS,
        -CARBON[mass_type_]-OXYGEN[mass_type_]-H2O[mass_type_], 1);
  else if (!ion_name.compare("im"))
    Set(ion_name, IMMONIUM_IONS, -CARBON[mass_type_]-OXYGEN[mass_type_], 1);
  else if (!ion_name.compare("iya"))
    Set(ion_name, INTERNAL_IONS, -CARBON[mass_type_]-OXYGEN[mass_type_], 1);
  else if (!ion_name.compare("iyb"))
    Set(ion_name, INTERNAL_IONS, 0, 1);
  else
    fprintf(stderr, "- Incorrect ion type: %s.\n", ion_name.c_str());
}

bool IonType::IsContainingPhospho(const Peptide& peptide,
                                  int cleavage_site) const {
  size_t start;
  size_t len = cleavage_site + 1;
  if (ion_category_ == NTERM_IONS)
    start = 0;
  else if (ion_category() == CTERM_IONS)
    start = peptide.sequence().size() - cleavage_site - 1;
  unsigned short mod_id = 21; // phosphorylation
  return peptide.IsContainingMod(mod_id, start, len);
}

bool IonType::IsContainingSTED(const Peptide& peptide,
                               int cleavage_site) const {
  size_t start;
  size_t len = cleavage_site + 1;
  if (ion_category_ == NTERM_IONS)
    start = 0;
  else if (ion_category_ == CTERM_IONS)
    start = peptide.sequence().size() - cleavage_site - 1;
  return (peptide.IsContainingAA('S', start, len) ||
          peptide.IsContainingAA('T', start, len) ||
          peptide.IsContainingAA('E', start, len) ||
          peptide.IsContainingAA('D', start, len));
}

int IonType::IndexToCleavageSite(int index, int pep_len) const {
  int clvg_site;
  if (ion_category_ == NTERM_IONS)
    clvg_site = index;
  else if (ion_category_ == CTERM_IONS)
    clvg_site = pep_len - index - 1;
  return clvg_site;
}

} // namespace phos_loc
