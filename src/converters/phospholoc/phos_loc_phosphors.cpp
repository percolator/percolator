#include "phos_loc_phosphors.h"

namespace phos_loc {

PhosphoRS::PhosphoRS() {
}

PhosphoRS::PhosphoRS(const Spectrum& spec, const std::string& pep_seq,
                     std::vector<std::vector<LocationMod> >& input_var_mod_combs,
                     const Parameters& paras) {
  InitPhosphoRS(spec, pep_seq, input_var_mod_combs, paras);
}

PhosphoRS::~PhosphoRS() {
}

double PhosphoRS::GetPeptideScore(std::vector<int>& phospho_locs) {
  int idx = -1;
  for (std::vector<std::vector<int> >::size_type i = 0;
       i < all_phospho_site_combinations_.size(); ++i) {
    if (AreEqualVectors(phospho_locs, all_phospho_site_combinations_[i])) {
      idx = i;
      break;
    }
  }
  if (idx >= 0)
    return peptide_scores_[idx];
  else
    return 0.0;
}

double PhosphoRS::GetSequenceProbability(std::vector<int>& phospho_locs) {
  int idx = -1;
  for (std::vector<std::vector<int> >::size_type i = 0;
       i < all_phospho_site_combinations_.size(); ++i) {
    if (AreEqualVectors(phospho_locs, all_phospho_site_combinations_[i])) {
      idx = i;
      break;
    }
  }
  if (idx >= 0)
    return sequence_probs_[idx];
  else
    return 0.0;
}

double PhosphoRS::GetPhosphoSiteScore(std::vector<int>& phospho_locs) {
  double site_score = 0.0;
  for (std::vector<int>::const_iterator it = phospho_locs.begin();
       it != phospho_locs.end(); ++it) {
    site_score += phospho_site_probs_[*it];
  }
  return site_score / phospho_locs.size();
}

double PhosphoRS::min_peak_depth_ = 2;
double PhosphoRS::max_peak_depth_ = 8;

void PhosphoRS::InitPhosphoRS(
    const Spectrum& spec,
    const std::string pep_seq,
    std::vector<std::vector<LocationMod> >& input_var_mod_combs,
    const Parameters& paras) {
  // set all phospho-site combinations
  std::vector<int> phospho_locs = GetPhosphoedLocations(input_var_mod_combs[0]);
  int num_phospho_sites = phospho_locs.size();
  std::vector<int> potential_sites = GetPotentialPhosphoSites(pep_seq);
  SetAllPhosphoSiteCombinations(potential_sites, num_phospho_sites);
  // set all phospho-PSMs
  std::vector<LocationMod> loc_mods(input_var_mod_combs[0]);
  for (std::vector<std::vector<int> >::const_iterator it =
           all_phospho_site_combinations_.begin();
       it != all_phospho_site_combinations_.end(); ++it) {
    ResetPhosphoedSites(*it, loc_mods);
    IonSeries ion_series(pep_seq, loc_mods, paras, spec.precursor().charge);
    Match psm(spec, ion_series, paras.frag_tolerance_);
    pep_spec_matches_.push_back(psm);
  }
  // save phospho-site combinations actually occurred in search result
  for (std::vector<std::vector<LocationMod> >::const_iterator it =
           input_var_mod_combs.begin();
       it != input_var_mod_combs.end(); ++it) {
    phospho_locs = GetPhosphoedLocations(input_var_mod_combs[0]);
    input_phospho_site_combinations_.push_back(phospho_locs);
  }
  // get the indexes of input phospho-sites combs in 'all_phospho_site_combinations_'
  MapInputSiteCombinationsInAll();
  // find optimal peak depth for each m/z window separately
  SelectOptimalPeakDepthsForAllWindows();
  // start scoring
  InitPeptideScoresAndSequenceProbs();
  InitPhosphoSiteProbabilities();

}

std::vector<int> PhosphoRS::GetPhosphoedLocations(
    const std::vector<LocationMod>& loc_mods) {
  std::vector<int> phosphoed_locs;
  for (std::vector<LocationMod>::const_iterator it = loc_mods.begin();
       it != loc_mods.end(); ++it) {
    if (it->mod_id == UNIMOD_PHOSPHO_ID)
      phosphoed_locs.push_back(it->aa_idx);
  }
  return phosphoed_locs;
}

std::vector<int> PhosphoRS::GetPotentialPhosphoSites(
    const std::string& pep_seq) {
  std::vector<int> sites;
  for (std::string::size_type i = 0; i < pep_seq.size(); ++i) {
    if (pep_seq[i] == 'S' || pep_seq[i] == 'T' || pep_seq[i] == 'Y')
      sites.push_back(i);
  }
  return sites;
}

void PhosphoRS::SetAllPhosphoSiteCombinations(
    std::vector<int>& potential_phospho_sites, int actual_mod_num) {
  do {
    std::vector<int> comb(potential_phospho_sites.begin(),
                          potential_phospho_sites.begin() + actual_mod_num);
    all_phospho_site_combinations_.push_back(comb);
  } while (phos_loc::NextCombination(potential_phospho_sites.begin(),
                                  potential_phospho_sites.begin() + actual_mod_num,
                                  potential_phospho_sites.end()));
}

void PhosphoRS::ResetPhosphoedSites(const std::vector<int>& one_site_comb,
                                    std::vector<LocationMod>& loc_mods) {
  std::vector<LocationMod> new_loc_mods;
  for (std::vector<LocationMod>::const_iterator it = loc_mods.begin();
       it != loc_mods.end(); ++it) {
    if (it->mod_id != UNIMOD_PHOSPHO_ID)
      new_loc_mods.push_back(*it);
  }
  for (std::vector<int>::const_iterator it = one_site_comb.begin();
       it != one_site_comb.end(); ++it) {
    new_loc_mods.push_back(LocationMod(*it, UNIMOD_PHOSPHO_ID));
  }
  loc_mods.clear();
  loc_mods.assign(new_loc_mods.begin(), new_loc_mods.end());
}

void PhosphoRS::MapInputSiteCombinationsInAll() {
  for (std::vector<std::vector<int> >::iterator it =
           input_phospho_site_combinations_.begin();
       it != input_phospho_site_combinations_.end(); ++it) {
    for (std::vector<std::vector<int> >::size_type i = 0;
         i < all_phospho_site_combinations_.size(); ++i) {
      if (AreEqualVectors(*it, all_phospho_site_combinations_[i])) {
        indexes_in_all_.push_back(i);
        break;
      }
    }
  }
}

void PhosphoRS::InitScoreRanksCurrentWindow() {
  for (std::vector<Match>::size_type i = 0;
       i < pep_spec_matches_.size(); ++i) {
    index_scores_current_window_.push_back(i);
  }
}

void PhosphoRS::SortScoresAtPeakDepth(int peak_depth) {
  std::vector<double> scores_at_peak_depth;
  std::vector<std::vector<double> >::const_iterator it;
  for (it = peptide_scores_current_window_.begin();
       it != peptide_scores_current_window_.end(); ++it) {
    scores_at_peak_depth.push_back(it->at(peak_depth - min_peak_depth_));
  }
  InitScoreRanksCurrentWindow();
  std::sort(index_scores_current_window_.begin(),
            index_scores_current_window_.end(),
            GreaterScore(scores_at_peak_depth));
}

double PhosphoRS::DeltaScoreAtPeakDepthAndRank(int peak_depth, int score_rank) {
  SortScoresAtPeakDepth(peak_depth);
  int top_idx = index_scores_current_window_[0];
  int low_idx = index_scores_current_window_[score_rank];
  double top_score =
      peptide_scores_current_window_[top_idx][peak_depth-min_peak_depth_];
  double low_score =
      peptide_scores_current_window_[low_idx][peak_depth-min_peak_depth_];
  return top_score - low_score;
}

int PhosphoRS::PeakDepthAtMaxScore(std::vector<double>& scores_all_depths) {
  double peak_depth = min_peak_depth_;
  double max_score = 0.0;
  for (std::vector<double>::size_type i = 0;
       i < scores_all_depths.size(); ++i) {
    if (scores_all_depths[i] > max_score) {
      max_score = scores_all_depths[i];
      peak_depth = i + min_peak_depth_;
    }
  }
  return peak_depth;
}

bool PhosphoRS::AreEqualScores(double score1, double score2, double tol) {
 return abs(score1 - score2) < tol;
}

template<typename T>
bool PhosphoRS::AreEqualVectors(std::vector<T>& v1, std::vector<T>& v2) {
  std::sort(v1.begin(), v1.end());
  std::sort(v2.begin(), v2.end());
  return v1 == v2;
}

std::vector<int> PhosphoRS::PeakDepthsAtMaxDeltaScore(
    int rank_to_compare, std::vector<int>& candidate_peak_depths) {
  double max_delta_score = 0.0;
  std::vector<int> updated_candidate_peak_depths;
  for (std::vector<int>::const_iterator it = candidate_peak_depths.begin();
       it != candidate_peak_depths.end(); ++it) {
    double curr_delta = DeltaScoreAtPeakDepthAndRank(*it, rank_to_compare);
    if (AreEqualScores(max_delta_score, curr_delta, 0.001)) {
      updated_candidate_peak_depths.push_back(*it);
    }
    else if (curr_delta > max_delta_score) {
      max_delta_score = curr_delta;
      updated_candidate_peak_depths.clear();
      updated_candidate_peak_depths.push_back(*it);
    }
  }
  return updated_candidate_peak_depths;
}

int PhosphoRS::SelectOptimalPeakDepthInOneWindow() {
  int optimal_peak_depth = min_peak_depth_;
  std::vector<int> candidate_peak_depths;
  for (int pd = min_peak_depth_; pd <= max_peak_depth_; ++pd) {
    candidate_peak_depths.push_back(pd);
  }
  // find peak depth that makes peptide score maximal if only one combination
  if (peptide_scores_current_window_.size() == 1) {
    optimal_peak_depth = PeakDepthAtMaxScore(peptide_scores_current_window_[0]);
    return optimal_peak_depth;
  }
  // find peak depth that make delta peptide score maximal
  int max_rank_to_compare = 3;
  if (peptide_scores_current_window_.size() <= 3)
    max_rank_to_compare = peptide_scores_current_window_.size();
  for (int curr_rank = 1; curr_rank <= max_rank_to_compare; ++curr_rank) {
    candidate_peak_depths =
        PeakDepthsAtMaxDeltaScore(curr_rank, candidate_peak_depths);
    // if optimal depth is uniquely determined
    if (candidate_peak_depths.size() <= 1)
      break;
    SortScoresAtPeakDepth(candidate_peak_depths[0]);
    optimal_peak_depth = PeakDepthAtMaxScore(
        peptide_scores_current_window_[index_scores_current_window_[0]]);
  }

  if (candidate_peak_depths.size() > 1) {
    SortScoresAtPeakDepth(optimal_peak_depth);
    optimal_peak_depth = PeakDepthAtMaxScore(
        peptide_scores_current_window_[index_scores_current_window_[0]]);
  }
  else {
    optimal_peak_depth = candidate_peak_depths[0];
  }
  return optimal_peak_depth;
}

void PhosphoRS::SelectOptimalPeakDepthsForAllWindows() {
  std::vector<double> window_bounds =
      pep_spec_matches_[0].spectrum().window_bounds();
  for (std::vector<double>::size_type i = 1; i < window_bounds.size(); ++i) {
    double lower_bnd = window_bounds[i-1];
    double upper_bnd = window_bounds[i];
    peptide_scores_current_window_.clear();
    for (std::vector<Match>::const_iterator it_m = pep_spec_matches_.begin();
         it_m != pep_spec_matches_.end(); ++it_m) {
      std::vector<double> curr_pep_scores;
      for (int curr_peak_depth = min_peak_depth_;
           curr_peak_depth <= max_peak_depth_; ++curr_peak_depth) {
        double prob = it_m->PeptideScoreForPhosphoRS(curr_peak_depth,
                                                     lower_bnd, upper_bnd);
        curr_pep_scores.push_back(-10*log10(prob));
      }
      peptide_scores_current_window_.push_back(curr_pep_scores);
    }
    int optimal_pd = SelectOptimalPeakDepthInOneWindow();
    optimal_peak_depths_.push_back(optimal_pd);
  }
}

void PhosphoRS::InitPeptideScoresAndSequenceProbs() {
  double sum_seq_probs = 0;
  for (std::vector<Match>::const_iterator it = pep_spec_matches_.begin();
       it != pep_spec_matches_.end(); ++it) {
    double prob = it->PeptideScoreForPhosphoRS(optimal_peak_depths_);
    peptide_scores_.push_back(-10*log10(prob));
    sequence_probs_.push_back(1/prob);
    sum_seq_probs += sequence_probs_.back();
  }
  for (std::vector<double>::iterator it = sequence_probs_.begin();
       it != sequence_probs_.end(); ++it) {
    *it /= sum_seq_probs;
  }
}

void PhosphoRS::InitPhosphoSiteProbabilities() {
  std::string pep_seq = pep_spec_matches_[0].ion_series().peptide().sequence();
  std::vector<int> potential_sites = GetPotentialPhosphoSites(pep_seq);
  for (std::vector<int>::const_iterator it = potential_sites.begin();
       it != potential_sites.end(); ++it) {
    phospho_site_probs_.insert(std::pair<int, double>(*it, 0.0));
  }

  for (size_t i = 0; i < sequence_probs_.size(); ++i) {
    std::vector<int> curr_phos_sites = all_phospho_site_combinations_[i];
    for (std::vector<int>::const_iterator it = curr_phos_sites.begin();
         it != curr_phos_sites.end(); ++it) {
      phospho_site_probs_[*it] += sequence_probs_[i];
    }
  }
}

} // namespace phos_loc
