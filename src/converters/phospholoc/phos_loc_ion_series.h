#ifndef IONSERIES_H_
#define IONSERIES_H_

#include "phos_loc_peptide.h"
#include "phos_loc_ion.h"
#include <vector>
#include <map>
#include <algorithm>

namespace phos_loc {

class IonSeries {
 public:
  IonSeries();
  IonSeries(const Peptide& peptide, const Parameters& paras,
            int precursor_charge);
  IonSeries(const std::string& pep_seq, const std::vector<LocationMod>& loc_mods,
            const Parameters& paras, int precursor_charge,
            char nterm_aa = '\0', char cterm_aa = '\0');
  IonSeries(const IonSeries& ion_series);
  virtual ~IonSeries();

  IonSeries& operator=(const IonSeries& ion_series);
  const Ion& operator()(const IonType& ion_type, int clvg_site);

  const Peptide& peptide() const { return peptide_; }
  const std::map<IonType, std::vector<Ion> >& ion_matrix() const {
    return ion_matrix_;
  }
  const Ion& ion_matrix(const IonType& ion_type, int clvg_site) {
    return ion_matrix_[ion_type][clvg_site];
  }
  const std::vector<Ion>& sorted_ions() const { return sorted_ions_; };

  int GetNumPredictedIons() const { return sorted_ions_.size(); };
  int GetNumPredictedIons(int lower_site, int upper_site) const;
  int GetNumPredictedIons(double lower_bnd, double upper_bnd) const;

 private:
  void InitIonMatrix(const Parameters& paras, int precursor_charge);
  void InitIonVector(const IonType& ion_type, std::vector<Ion>& ions);
  void InitSortedIons();
  static bool LessIonMZ(const Ion& ion1, const Ion& ion2);

 protected:
  Peptide peptide_;
  std::map<IonType, std::vector<Ion> > ion_matrix_;
  std::vector<Ion> sorted_ions_;
};

} // namespace phos_loc

#endif // IONSERIES_H_
