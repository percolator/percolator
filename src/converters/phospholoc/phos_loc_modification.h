#ifndef MODIFICATION_H_
#define MODIFICATION_H_

#include "phos_loc_misc.h"
#include <cstdio>
#include <cstring>
#include <string>

namespace phos_loc {

class Modification {
 public:
  Modification();
  Modification(const std::string name, char max_num_per_pep = 3);
  Modification(unsigned short unimod_id, char max_num_per_pep = 3);
  Modification(const std::string name, unsigned short unimod_id,
               double* mass_shift, double* neutral_loss,
               const std::string mod_aas, ModificationPosition position,
               char max_num_per_pep = 3);
  Modification(const Modification& modification);
  virtual ~Modification();

  const std::string name() const { return name_; }
  unsigned short unimod_id() const { return unimod_id_; }
  double mass_shift(MassType mass_t) const { return mass_shift_[mass_t]; }
  double neutral_loss(MassType mass_t) const { return neutral_loss_[mass_t]; }
  bool modifying(int aa_order) const { return modifying_[aa_order]; }
  bool modifying(char aa) const { return modifying_[aa-'A']; }
  ModificationPosition position() const { return position_; }
  char max_num_per_peptide() const { return max_num_per_peptide_; }

 private:
  void Init(const std::string name, unsigned short unimod_id,
            double mass_shift_aver, double mass_shift_mono,
            double neutral_loss_aver, double neutral_loss_mono,
            const std::string mod_aas, ModificationPosition position,
            char max_num);
  void MapModAAs(const std::string mod_aas);

 protected:
  std::string name_;
  unsigned short unimod_id_;
  double mass_shift_[NUM_MASS_TYPES];
  double neutral_loss_[NUM_MASS_TYPES];
  bool modifying_[LENGTH_AMINO_ACID_TABLE];
  ModificationPosition position_;
  char max_num_per_peptide_; // maximum number of this var_mod per peptide
};

} // namespace phos_loc

#endif // MODIFICATION_H_
