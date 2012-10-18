#include "ion_type.h"

IonType::IonType()
    : ion_category_(NTERM_IONS),
      mass_type_(MONO),
      mass_shift_(0),
      charge_(0) {
}

IonType::IonType(const IonType& ion_type)
    : ion_name_(ion_type.ion_name_),
      ion_category_(ion_type.ion_category_),
      mass_type_(ion_type.mass_type_),
      mass_shift_(ion_type.mass_shift_),
      charge_(ion_type.charge_) {
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
  }
  return *this;
}

bool IonType::operator<(const IonType& ion_type) const {
  if (ion_name_ != ion_type.ion_name_)
    return ion_name_ < ion_type.ion_name_;
  else
    return charge_ < ion_type.charge_;
}

double IonType::GetMoverZ(double molecular_mass_all_neutral_residues) const {
  return (molecular_mass_all_neutral_residues + mass_shift_ + PROTON * charge_) /
      charge_;
}

void IonType::Set(const std::string& ion_name, IonCategory ion_category,
                  double mass_shift, int charge) {
  ion_name_ = ion_name;
  ion_category_ = ion_category;
  mass_shift_ = mass_shift;
  charge_ = charge;
}

void IonType::InterpretIonName(const std::string ion_name) {
  if (!ion_name.compare("b"))
    Set("b", NTERM_IONS, 0, 1);
  else if (!ion_name.compare("y"))
    Set("y", CTERM_IONS, H2O[mass_type_], 1);
  else if (!ion_name.compare("c"))
    Set("c", NTERM_IONS, NH3[mass_type_], 1);
  else if (!ion_name.compare("c-h"))
    Set("c", NTERM_IONS, NH3[mass_type_]-HYDROGEN[mass_type_], 1);
  else if (!ion_name.compare("z"))
    Set("z", CTERM_IONS,
        -HYDROGEN[mass_type_]-NITROGEN[mass_type_]+OXYGEN[mass_type_], 1);
  else if (!ion_name.compare("z+h"))
    Set("z+h", CTERM_IONS, -NITROGEN[mass_type_]+OXYGEN[mass_type_], 1);
  else if (!ion_name.compare("z+2h"))
    Set("z+2h", CTERM_IONS,
        HYDROGEN[mass_type_]-NITROGEN[mass_type_]+OXYGEN[mass_type_], 1);
  else if (!ion_name.compare("b-nh3"))
    Set("b-nh3", NTERM_IONS, -NH3[mass_type_], 1);
  else if (!ion_name.compare("b-h2o"))
    Set("b-h2o", NTERM_IONS, -H2O[mass_type_], 1);
  else if (!ion_name.compare("y-nh3"))
    Set("y-nh3", CTERM_IONS, H2O[mass_type_]-NH3[mass_type_], 1);
  else if (!ion_name.compare("y-h2o"))
    Set("y-h2o", CTERM_IONS, 0, 1);
  else if (!ion_name.compare("a"))
    Set("a", NTERM_IONS, -CARBON[mass_type_]-OXYGEN[mass_type_], 1);
  else if (!ion_name.compare("a-nh3"))
    Set("a-nh3", NTERM_IONS,
        -CARBON[mass_type_]-OXYGEN[mass_type_]-NH3[mass_type_], 1);
  else if (!ion_name.compare("a-h2o"))
    Set("a-h2o", NTERM_IONS,
        -CARBON[mass_type_]-OXYGEN[mass_type_]-H2O[mass_type_], 1);
  else if (!ion_name.compare("im"))
    Set("im", IMMONIUM_IONS, -CARBON[mass_type_]-OXYGEN[mass_type_], 1);
  else if (!ion_name.compare("iya"))
    Set("iya", INTERNAL_IONS, -CARBON[mass_type_]-OXYGEN[mass_type_], 1);
  else if (!ion_name.compare("iyb"))
    Set("iya", INTERNAL_IONS, 0, 1);
  else
    fprintf(stderr, "- Incorrect ion type: %s.\n", ion_name.c_str());
}


