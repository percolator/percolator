#include "match.h"

Match::Match() {

}

Match::Match(const Spectrum& spec, const IonSeries& ions,
             const Tolerance& tol)
    : spectrum_(spec),
      ion_series_(ions),
      frag_tol_(tol) {
}

Match::Match(const Match& match)
    : spectrum_(match.spectrum_),
      ion_series_(match.ion_series_),
      frag_tol_(match.frag_tol_),
      peak_to_ions_(match.peak_to_ions_),
      ion_to_peaks_(match.ion_to_peaks_) {
}

Match::~Match() {
}

std::vector<int> Match::PeaksForSingleIon(int start_peak_idx, double ion_mz) {
  std::vector<int> matched_peak_indexes;
  for (std::vector<Peak>::size_type i = start_peak_idx;
       i < spectrum_.peaks().size(); ++i) {
    double peak_mz = spectrum_[i].m_over_z();
    if (IsMatched(peak_mz, ion_mz))
      matched_peak_indexes.push_back(i);
    else if (peak_mz > ion_mz)
      break;
  }
  return matched_peak_indexes;
}

void Match::FindAllPeakIonMatches() {
  int curr_peak_idx = 0;
  std::vector<int> pks_idx;
  InitIonMatchMatrix();
  for (std::vector<Ion>::size_type i = 0;
       i < ion_series_.sorted_ions().size(); ++i) {
    double ion_mz = ion_series_.sorted_ions()[i].m_over_z();
    pks_idx = PeaksForSingleIon(curr_peak_idx, ion_mz);
    if (pks_idx.empty())
      continue;
    curr_peak_idx = pks_idx[0];
    InsertMatchedPeaksForIon(i, pks_idx);
    InsertMatchedIonForPeaks(i, pks_idx);
  }
}

double Match::ProbabilityOfAPeakIonMatch(int peak_depth, double window_width) const {
  return frag_tol_.value * 2 * peak_depth / window_width;
}

double Match::PeptideScoreForAscore(int peak_depth, int lower_site,
                                    int upper_site, double window_width) const {
  double p = ProbabilityOfAPeakIonMatch(peak_depth, window_width);
  int num_pred_ions = GetNumPredictedIons(lower_site, upper_site);
  int num_match_ions = GetNumMatchedIons(peak_depth, lower_site, upper_site);
  if (num_match_ions == 0) return 1.0;
  return utils::CumulativeBinomialProbability(num_pred_ions, num_match_ions, p);
}

double Match::PeptideScoreForAscore(int peak_depth, double window_width) const {
  double p = ProbabilityOfAPeakIonMatch(peak_depth, window_width);
  int num_pred_ions = GetNumPredictedIons();
  int num_match_ions = GetNumMatchedIons(peak_depth);
  if (num_match_ions == 0) return 1.0;
  return utils::CumulativeBinomialProbability(num_pred_ions, num_match_ions, p);
}

void Match::Print() {
  fprintf(stderr, "Matched peaks:\n");
  std::map<int, std::vector<int> >::const_iterator it_p2i;
  for (it_p2i = peak_to_ions_.begin(); it_p2i != peak_to_ions_.end(); ++it_p2i) {
    fprintf(stderr, "%f %f: ", spectrum_[it_p2i->first].m_over_z(),
            spectrum_[it_p2i->first].intensity());
    for (std::vector<int>::const_iterator it_i = it_p2i->second.begin();
         it_i != it_p2i->second.end(); ++it_i) {
      fprintf(stderr, "%s ", ion_series_.sorted_ions()[*it_i].GetIonLabel().c_str());
    }
    fprintf(stderr, "\n");
  }
  fprintf(stderr, "Matched ions:\n");
  std::map<int, std::vector<int> >::const_iterator it_i2p;
  for (it_i2p = ion_to_peaks_.begin(); it_i2p != ion_to_peaks_.end(); ++it_i2p) {
    fprintf(stderr, "%s: ",
            ion_series_.sorted_ions()[it_i2p->first].GetIonLabel().c_str());
    for (std::vector<int>::const_iterator it_p = it_i2p->second.begin();
         it_p != it_i2p->second.end(); ++it_p)
      fprintf(stderr, "%.4f|%.4f ", spectrum_[*it_p].m_over_z(),
              spectrum_[*it_p].intensity());
    fprintf(stderr, "\n");
  }
}

void Match::InitIonMatchMatrix() {
  std::map<IonType, std::vector<Ion> >::const_iterator it;
  std::vector<int> matching_peak_idx; // empty
  int len_match_row = ion_series_.peptide().sequence().size() - 1;
  std::vector<std::vector<int> > match_row(len_match_row, matching_peak_idx);
  for (it = ion_series_.ion_matrix().begin();
       it != ion_series_.ion_matrix().end(); ++it) {
    ion_match_matrix_.insert(make_pair(it->first, match_row));
  }
}

bool Match::IsMatched(double peak_mz, double ion_mz) {
  return fabs(peak_mz + frag_tol_.drift - ion_mz) <= frag_tol_.value;
}

int Match::IndexOfHighestPeak(std::vector<int>& peak_indexes) {
  int idx_highest = 0;
  double max_inten = 0;
  for (std::vector<int>::size_type i = 0; i < peak_indexes.size(); ++i) {
    if (spectrum_[peak_indexes[i]].intensity() > max_inten) {
      max_inten = spectrum_[peak_indexes[i]].intensity();
      idx_highest = i;
    }
  }
  return idx_highest;
}

void Match::InsertMatchedPeaksForIon(int ion_idx, std::vector<int>& peak_indexes,
                                     bool only_keep_highest) {
  // get index of highest matching peak
  std::vector<int> new_peak_indexes;
  if (only_keep_highest) {
    int idx_highest = IndexOfHighestPeak(peak_indexes);
    new_peak_indexes.push_back(peak_indexes[idx_highest]);
  }
  else {
    new_peak_indexes.assign(peak_indexes.begin(),peak_indexes.end());
  }
  // label matched ion
  IonType ion_t = ion_series_.sorted_ions()[ion_idx].ion_type();
  int clvg_site = ion_series_.sorted_ions()[ion_idx].cleavage_site();
  if (ion_to_peaks_.find(ion_idx) == ion_to_peaks_.end()) {
    ion_to_peaks_.insert(make_pair(ion_idx, new_peak_indexes));
    ion_match_matrix_[ion_t][clvg_site] = new_peak_indexes;
  }
  else {
    ion_to_peaks_[ion_idx].insert(ion_to_peaks_[ion_idx].end(),
                                  new_peak_indexes.begin(),
                                  new_peak_indexes.end());
    ion_match_matrix_[ion_t][clvg_site] = new_peak_indexes;
  }
}

void Match::InsertMatchedIonForPeaks(int ion_idx, std::vector<int>& peak_indexes,
                                     bool only_keep_highest) {
  // get index of highest matching peak
  std::vector<int> new_peak_index;
  if (only_keep_highest) {
    int idx_highest = IndexOfHighestPeak(peak_indexes);
    new_peak_index.push_back(peak_indexes[idx_highest]);
  }
  else {
    new_peak_index.assign(peak_indexes.begin(),peak_indexes.end());
  }
  // label matching peak
  for (std::vector<int>::const_iterator it = new_peak_index.begin();
       it != new_peak_index.end(); ++it) {
    if (peak_to_ions_.find(*it) == peak_to_ions_.end()) {
      std::vector<int> ion_indexes(1, ion_idx);
      peak_to_ions_.insert(make_pair(*it, ion_indexes));
    }
    else {
      peak_to_ions_[*it].push_back(ion_idx);
    }
  }
}

int Match::GetNumMatchedIons(int peak_depth, int lower_site, int upper_site) const {
  if (lower_site > upper_site) {
    double swap = lower_site;
    lower_site = upper_site;
    upper_site = swap;
  }
  int num_matches = 0;
  std::map<IonType, std::vector<std::vector<int> > >::const_iterator it;
  for (it = ion_match_matrix_.begin(); it != ion_match_matrix_.begin(); ++it) {
    size_t lower_idx = lower_site;
    size_t upper_idx = upper_site;
    if (it->first.ion_category() == CTERM_IONS) {
      int pep_len = ion_series_.peptide().sequence().size();
      lower_idx = pep_len - upper_site - 2;
      upper_idx = pep_len - lower_site - 2;
    }
    for (std::vector<std::vector<int> >::size_type i = lower_idx;
         i < upper_idx; ++i) {
      for(std::vector<int>::const_iterator it_pk = it->second[i].begin();
          it_pk != it->second[i].end(); ++it_pk) {
        if (spectrum_.peaks()[*it_pk].rank() < peak_depth) {
          ++num_matches;
          break;
        }
      }
    }
  }
  return num_matches;
}

int Match::GetNumMatchedIons(int peak_depth) const {
  int num_matches = 0;
  for (std::map<int, std::vector<int> >::const_iterator it = ion_to_peaks_.begin();
       it != ion_to_peaks_.end(); ++it) {
    for (std::vector<int>::const_iterator it_pk = it->second.begin();
         it_pk != it->second.end(); ++it_pk) {
      if (spectrum_.peaks()[*it_pk].rank() < peak_depth) {
        ++num_matches;
        break;
      }
    }
  }
  return num_matches;
}

int Match::GetNumPredictedIons(int lower_site, int upper_site) const {
  if (lower_site > upper_site) {
    double swap = lower_site;
    lower_site = upper_site;
    upper_site = swap;
  }
  return ion_match_matrix_.size() * (upper_site - lower_site);
}

int Match::GetNumPredictedIons() const {
  return ion_series_.GetNumPredictedIons();
}
