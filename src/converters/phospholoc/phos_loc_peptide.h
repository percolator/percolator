#ifndef PEPTIDE_H_
#define PEPTIDE_H_

#include "phos_loc_misc.h"
#include "phos_loc_parameters.h"
#include <string>
#include <vector>

namespace phos_loc {

class Peptide {
 public:
  Peptide();
  Peptide(const std::string& pep_seq, const Parameters& paras,
          char nterm_aa = '\0', char cterm_aa = '\0');
  Peptide(const std::string& pep_seq, const std::vector<LocationMod>& loc_mods,
          const Parameters& paras, char nterm_aa = '\0', char cterm_aa = '\0');
  Peptide(const Peptide& peptide);
  virtual ~Peptide();

  Peptide& operator=(const Peptide& peptide);

  const std::string& sequence() const { return sequence_; }
  const std::vector<LocationMod>& location_mods() const { return location_mods_; }
  char nterm_aa() const { return nterm_aa_; }
  char cterm_aa() const { return cterm_aa_; }
  double molecular_mass() const { return molecular_mass_; }
  double nterm_ladder(int clvg_site) const { return nterm_ladder_[clvg_site]; }
  double cterm_ladder(int clvg_site) const { return cterm_ladder_[clvg_site]; }

  void set_location_mods(const std::vector<LocationMod>& mods) {
    location_mods_.assign(mods.begin(), mods.end());
  }

  bool IsContainingAA(char aa, size_t start, size_t len) const;
  bool IsContainingMod(unsigned short mod_id, size_t start, size_t len) const;
  void ClearVarMods() { location_mods_.clear(); }
  void AddVarMod(unsigned char aa_idx, unsigned char mod_id) {
    location_mods_.push_back(LocationMod(aa_idx, mod_id));
  }

 private:
  void InitVarModMass(double* var_mod_mass, bool is_fragment,
                      const Parameters& paras);
  void InitMolecularMass(const Parameters& paras);
  void InitNtermLadder(const Parameters& paras);
  void InitCtermLadder(const Parameters& paras);

 private:
  std::string sequence_;
  std::vector<LocationMod> location_mods_; // locations of variable modifications
  char nterm_aa_;  // amino acid preceding the N-term cleavage site
  char cterm_aa_;
  double molecular_mass_;  // total mass of all AA residues + H2O
  std::vector<double> nterm_ladder_;  // accumulative N-term AA ladder masses
  std::vector<double> cterm_ladder_;
};

} // namespace phos_loc

#endif // PEPTIDE_H_
