#ifndef IONTYPE_H_
#define IONTYPE_H_

#include "misc.h"
#include <cstdio>
#include <string>

class IonType {
 public:
  IonType();
  IonType(const IonType& ion_type);
  IonType(const std::string& ion_name, MassType mass_type = MONO);
  virtual ~IonType();

  IonType& operator=(const IonType& ion_type);
  bool operator<(const IonType& ion_type) const;

  std::string ion_name() const { return ion_name_; }
  IonCategory ion_category() const { return ion_category_; }
  double mass_shift() const { return mass_shift_; }
  char charge() const { return charge_; }

  void set_ion_name(const std::string& ion_name) { ion_name_ = ion_name; }
  void set_ion_category(IonCategory ion_category) { ion_category_ = ion_category; }
  void set_mass_shift(double mass_shift) { mass_shift_ = mass_shift; }
  void set_charge(int charge) { charge_ = charge; }

  double GetMoverZ(double molecular_mass_all_neutral_residues) const;
  void Set(const std::string& ion_name, IonCategory ion_category,
           double mass_shift, int charge);
  void InterpretIonName(const std::string ion_name);

 protected:
  std::string ion_name_;
  IonCategory ion_category_;
  MassType mass_type_;
  // mass shift in neutral state relative to total mass
  // of neutral residues in ions of a specific ion_category
  double mass_shift_;
  char charge_;
};

#endif // IONTYPE_H_
