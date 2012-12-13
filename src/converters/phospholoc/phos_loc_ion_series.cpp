#include "phos_loc_ion_series.h"

namespace phos_loc {

IonSeries::IonSeries() {
}

IonSeries::IonSeries(const Peptide& peptide, const Parameters& paras,
                     int precursor_charge)
    : peptide_(peptide) {
  InitIonMatrix(paras, precursor_charge);
  InitSortedIons();
}

IonSeries::IonSeries(const std::string& pep_seq,
                     const std::vector<LocationMod>& loc_mods,
                     const Parameters& paras, int precursor_charge,
                     char nterm_aa, char cterm_aa)
    : peptide_(pep_seq, loc_mods, paras, nterm_aa, cterm_aa) {
  InitIonMatrix(paras, precursor_charge);
  InitSortedIons();
}

IonSeries::IonSeries(const IonSeries& ion_series)
    : peptide_(ion_series.peptide_),
      ion_matrix_(ion_series.ion_matrix_),
      sorted_ions_(ion_series.sorted_ions_) {
}

IonSeries::~IonSeries() {
}

IonSeries& IonSeries::operator=(const IonSeries& ion_series) {
  peptide_ = ion_series.peptide_;
  ion_matrix_ = ion_series.ion_matrix_;
  sorted_ions_ = ion_series.sorted_ions_;
  return *this;
}

const Ion& IonSeries::operator()(const IonType& ion_type, int clvg_site) {
  return ion_matrix_[ion_type][clvg_site];
}

int IonSeries::GetNumPredictedIons(int lower_site, int upper_site) const {
  if (lower_site > upper_site) // phosphoed upper_site not included
    std::swap(lower_site, upper_site);
  std::map<IonType, std::vector<Ion> >::const_iterator it;
  int num_pred = 0;
  for (it = ion_matrix_.begin(); it != ion_matrix_.end(); ++it) {
    if (it->first.NoCondition()) {
      num_pred += (upper_site - lower_site);
    }
    else {
      int lower_clvgsite =
          it->first.IndexToCleavageSite(lower_site, peptide_.sequence().size());
      int upper_clvgsite =
          it->first.IndexToCleavageSite(upper_site, peptide_.sequence().size());
      if (lower_clvgsite > upper_clvgsite)
        std::swap(lower_clvgsite, upper_clvgsite);
      std::vector<Ion>::const_iterator iter;
      for (iter = it->second.begin() + lower_clvgsite;
           iter != it->second.begin() + upper_clvgsite; ++iter) {
        if (iter->NotValid())
          continue;
        ++num_pred;
      }
    }
  }
  return num_pred;
}

int IonSeries::GetNumPredictedIons(double lower_bnd, double upper_bnd) const {
  if (lower_bnd > upper_bnd)
    std::swap(lower_bnd, upper_bnd);
  std::vector<Ion>::const_iterator it = sorted_ions_.begin();
  int num_pred = 0;
  while (it->m_over_z() < lower_bnd && it != sorted_ions_.end())
    ++it;
  while (it->m_over_z() < upper_bnd && it != sorted_ions_.end()) {
    ++num_pred;
    ++it;
  }
  return num_pred;
}

void IonSeries::InitIonMatrix(const Parameters& paras, int precursor_charge) {
  for (std::vector<IonType>::const_iterator it = paras.specified_ion_types_.begin();
       it != paras.specified_ion_types_.end(); ++it) {
    std::vector<Ion> ions;
    InitIonVector(IonType(*it), ions);
    ion_matrix_.insert(make_pair(*it, ions));
    ions.clear();
    // consider 2+ fragment ions
    if (precursor_charge > 2) {
      IonType ion_t(*it);
      ion_t.set_charge(2);
      InitIonVector(ion_t, ions);
      ion_matrix_.insert(make_pair(ion_t, ions));
      ions.clear();
    }
  }
}

void IonSeries::InitIonVector(const IonType& ion_type, std::vector<Ion>& ions) {
  for (std::string::size_type i = 0; i < peptide_.sequence().size()-1; ++i) {
    double mz = ion_type.GetMoverZ(peptide_, i);
    Ion ion(ion_type, i, mz, 0);
    ions.push_back(ion);
  }
}

void IonSeries::InitSortedIons() {
  std::map<IonType, std::vector<Ion> >::const_iterator it;
  for (it = ion_matrix_.begin(); it != ion_matrix_.end(); ++it) {
    if (it->first.NoCondition()) {
      sorted_ions_.insert(sorted_ions_.end(),
                          it->second.begin(), it->second.end());
    }
    else {
      for (std::vector<Ion>::const_iterator iter = it->second.begin();
           iter != it->second.end(); ++iter) {
        if (iter->NotValid())
          continue;
        sorted_ions_.push_back(*iter);
      }
    }
  }
  std::sort(sorted_ions_.begin(), sorted_ions_.end(), LessIonMZ);
}

bool IonSeries::LessIonMZ(const Ion& ion1, const Ion& ion2) {
  return ion1.m_over_z() < ion2.m_over_z();
}

} // namespace phos_loc
