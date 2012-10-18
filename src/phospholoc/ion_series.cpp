#include "ion_series.h"

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

void IonSeries::InitIonMatrix(const Parameters& paras, int precursor_charge) {
  for (std::vector<IonType>::const_iterator it = paras.specified_ion_types_.begin();
       it != paras.specified_ion_types_.end(); ++it) {
    std::vector<Ion> ions;
    InitIonVector(*it, ions);
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
  if (ion_type.ion_category() == NTERM_IONS)
    for (std::string::size_type i = 0; i < peptide_.sequence().size()-1; ++i) {
      double mz = ion_type.GetMoverZ(peptide_.nterm_ladder(i));
      Ion ion(ion_type, i, mz, 0);
      ions.push_back(ion);
    }
  else if (ion_type.ion_category() == CTERM_IONS)
    for (std::string::size_type i = 0; i < peptide_.sequence().size()-1; ++i) {
      double mz = ion_type.GetMoverZ(peptide_.cterm_ladder(i));
      Ion ion(ion_type, i, mz, 0);
      ions.push_back(ion);
    }
}

void IonSeries::InitSortedIons() {
  std::map<IonType, std::vector<Ion> >::iterator it;
  for (it = ion_matrix_.begin(); it != ion_matrix_.end(); ++it)
    sorted_ions_.insert(sorted_ions_.end(), it->second.begin(), it->second.end());
  std::sort(sorted_ions_.begin(), sorted_ions_.end(), LessIonMZ);
}

bool IonSeries::LessIonMZ(const Ion& ion1, const Ion& ion2) {
  return ion1.m_over_z() < ion2.m_over_z();
}
