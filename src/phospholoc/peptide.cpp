#include "peptide.h"

Peptide::Peptide()
    : nterm_aa_('\0'),
      cterm_aa_('\0'),
      molecular_mass_(0) {
}

Peptide::Peptide(const std::string& pep_seq, const Parameters& paras,
                 char nterm_aa, char cterm_aa)
    : sequence_(pep_seq),
      nterm_aa_(nterm_aa),
      cterm_aa_(cterm_aa) {
  InitMolecularMass(paras);
  InitNtermLadder(paras);
  InitCtermLadder(paras);
}

Peptide::Peptide(const std::string& pep_seq,
                 const std::vector<LocationMod>& loc_mods,
                 const Parameters& paras, char nterm_aa, char cterm_aa)
    : sequence_(pep_seq),
      location_mods_(loc_mods),
      nterm_aa_(nterm_aa),
      cterm_aa_(cterm_aa) {
  InitMolecularMass(paras);
  InitNtermLadder(paras);
  InitCtermLadder(paras);
}

Peptide::Peptide(const Peptide& peptide)
    : sequence_(peptide.sequence_),
      location_mods_(peptide.location_mods_),
      nterm_aa_(peptide.nterm_aa_),
      cterm_aa_(peptide.cterm_aa_),
      molecular_mass_(peptide.molecular_mass_),
      nterm_ladder_(peptide.nterm_ladder_),
      cterm_ladder_(peptide.cterm_ladder_) {
}

Peptide& Peptide::operator=(const Peptide& peptide) {
  if (this != &peptide) {
    sequence_ = peptide.sequence_;
    location_mods_.assign(peptide.location_mods_.begin(),
                          peptide.location_mods_.end());
    nterm_aa_ = peptide.nterm_aa_;
    cterm_aa_ = peptide.cterm_aa_;
    molecular_mass_ = peptide.molecular_mass_;
    nterm_ladder_.assign(peptide.nterm_ladder_.begin(), peptide.nterm_ladder_.end());
    cterm_ladder_.assign(peptide.cterm_ladder_.begin(), peptide.cterm_ladder_.end());
  }
  return *this;
}

Peptide::~Peptide() {
}

void Peptide::InitVarModMass(double* var_mod_mass, bool is_fragment,
                             const Parameters& paras) {
  if (is_fragment)
    for (std::vector<LocationMod>::const_iterator it = location_mods_.begin();
         it != location_mods_.end(); ++it)
      var_mod_mass[it->aa_idx] +=
          paras.variable_mods_[it->mod_id].mass_shift(
              paras.frag_tolerance_.mass_type);
  else
    for (std::vector<LocationMod>::const_iterator it = location_mods_.begin();
         it != location_mods_.end(); ++it)
      var_mod_mass[it->aa_idx] +=
          paras.variable_mods_[it->mod_id].mass_shift(
              paras.pep_tolerance_.mass_type);
}

void Peptide::InitMolecularMass(const Parameters& paras) {
  MassType mass_t = paras.pep_tolerance_.mass_type;
  molecular_mass_ = H2O[mass_t];
  double* var_mod_mass = new double[sequence_.size()]();
  InitVarModMass(var_mod_mass, false, paras);
  for (std::string::size_type i = 0; i < sequence_.size(); ++i) {
    molecular_mass_ +=
        paras.amino_acid_mass_[sequence_[i]-'A'][mass_t] + var_mod_mass[i];
  }
  delete[] var_mod_mass;
}

void Peptide::InitNtermLadder(const Parameters& paras) {
  MassType mass_t = paras.frag_tolerance_.mass_type;
  double mass = 0;
  double* var_mod_mass = new double[sequence_.size()]();
  InitVarModMass(var_mod_mass, true, paras);
  for (std::string::size_type i = 0; i < sequence_.size(); ++i) {
    mass += paras.amino_acid_mass_[sequence_[i]-'A'][mass_t] + var_mod_mass[i];
    nterm_ladder_.push_back(mass);
  }
  delete[] var_mod_mass;
}

void Peptide::InitCtermLadder(const Parameters& paras) {
  MassType mass_t = paras.frag_tolerance_.mass_type;
  double mass = 0;
  double* var_mod_mass = new double[sequence_.size()]();
  InitVarModMass(var_mod_mass, true, paras);
  for (std::string::size_type i = sequence_.size(); i > 0; --i) {
    // warning: size_type is unsigned, i>=0 always true
    mass += paras.amino_acid_mass_[sequence_[i-1]-'A'][mass_t] +
        var_mod_mass[i-1];
    cterm_ladder_.push_back(mass);
  }
  delete[] var_mod_mass;
}

