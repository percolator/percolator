#include "ascore.h"
#include "match.h"

Ascore::Ascore() {
}

Ascore::Ascore(const Spectrum& spec, const std::string& pep_seq,
               std::vector<std::vector<LocationMod> >& input_var_mod_combs,
               const Parameters& paras) {
  InitAscore(spec, pep_seq, input_var_mod_combs, paras);
}

Ascore::~Ascore() {
}

double Ascore::GetPeptideScore(std::vector<LocationMod>& input_var_mod_comb,
                               const Parameters& paras) {
  std::vector<int> phospho_locs = GetPhosphoedLocations(input_var_mod_comb,
                                                        paras.variable_mods_);
  int idx = -1;
  for (std::vector<std::vector<int> >::size_type i = 0;
       i < all_phospho_site_combinations_.size(); ++i) {
    if (IsEqual(phospho_locs, all_phospho_site_combinations_[i]))
      idx = i;
  }
  if (idx != -1)
    return weigthed_average_pep_scores_[idx];
  else
    return 0.0;
}

double Ascore::GetLocalPeptideScore(std::vector<LocationMod>& input_var_mod_comb,
                                    const Parameters& paras){
  std::vector<int> phospho_locs = GetPhosphoedLocations(input_var_mod_comb,
                                                        paras.variable_mods_);
  int idx = -1;
  for (std::vector<std::vector<int> >::size_type i = 0;
       i < input_phospho_site_combinations_.size(); ++i) {
    if (IsEqual(phospho_locs, input_phospho_site_combinations_[i]))
      idx = i;
  }
  if (idx != -1)
    return local_pep_scores_[idx];
  else
    return 0.0;
}

void Ascore::InitAscore(
    const Spectrum& spec,
    const std::string pep_seq,
    std::vector<std::vector<LocationMod> >& input_var_mod_combs,
    const Parameters& paras) {
  // set all phospho-site combinations
  std::vector<int> phospho_locs = GetPhosphoedLocations(input_var_mod_combs[0],
                                                        paras.variable_mods_);
  int num_phospho_sites = phospho_locs.size();
  std::vector<int> potential_sites = GetPotentialPhosphoSites(pep_seq);
  SetAllPhosphoSiteCombinations(potential_sites, num_phospho_sites);
  // set all phospho-PSMs
  std::vector<LocationMod> loc_mods(input_var_mod_combs[0]);
  for (std::vector<std::vector<int> >::const_iterator it =
           all_phospho_site_combinations_.begin();
       it != all_phospho_site_combinations_.end(); ++it) {
    ResetPhosphoedSites(*it, paras.variable_mods_, loc_mods);
    IonSeries ion_series(pep_seq, loc_mods, paras, spec.precursor().charge);
    Match psm(spec, ion_series, paras.frag_tolerance_);
    pep_spec_matches_.push_back(psm);
  }
  // save phospho-site combinations actually occurred in search result
  for (std::vector<std::vector<LocationMod> >::const_iterator it =
           input_var_mod_combs.begin();
       it != input_var_mod_combs.end(); ++it) {
    phospho_locs = GetPhosphoedLocations(input_var_mod_combs[0],
                                         paras.variable_mods_);
    input_phospho_site_combinations_.push_back(phospho_locs);
  }
  // get the indexes of input phospho-sites combs in 'all_phospho_site_combinations_'
  MapInputSiteCombinationsInAll();
  // calculate weighted average of peptide scores
  InitAllPeptideScores(paras);
  CalculateAllWeigthedAveragePeptideScores();
  SortWeightedAveragePeptideScores();
  // calculate local peptide scores for each phospho-peptide in search result
  InitPhosphoSiteProbabilities();
  InitLocalPepScores();
}

std::vector<int> Ascore::GetPhosphoedLocations(
    const std::vector<LocationMod>& loc_mods,
    const std::vector<Modification>& var_mods) {
  std::vector<int> phosphoed_locs;
  for (std::vector<LocationMod>::const_iterator it = loc_mods.begin();
       it != loc_mods.end(); ++it) {
    if (var_mods[it->mod_id].name() == "Phospho")
      phosphoed_locs.push_back(it->aa_idx);
  }
  return phosphoed_locs;
}

std::vector<int> Ascore::GetPotentialPhosphoSites(const std::string& pep_seq) {
  std::vector<int> sites;
  for (std::string::size_type i = 0; i < pep_seq.size(); ++i) {
    if (pep_seq[i] == 'S' || pep_seq[i] == 'T' || pep_seq[i] == 'Y')
      sites.push_back(i);
  }
  return sites;
}

void Ascore::SetAllPhosphoSiteCombinations(
    std::vector<int>& potential_phospho_sites, int actual_mod_num) {
  do {
    std::vector<int> comb(potential_phospho_sites.begin(),
                          potential_phospho_sites.begin() + actual_mod_num);
    all_phospho_site_combinations_.push_back(comb);
  } while (utils::NextCombination(potential_phospho_sites.begin(),
                                  potential_phospho_sites.begin() + actual_mod_num,
                                  potential_phospho_sites.end()));
}

void Ascore::ResetPhosphoedSites(const std::vector<int>& one_site_comb,
                                 const std::vector<Modification>& var_mods,
                                 std::vector<LocationMod>& loc_mods) {
  std::vector<LocationMod> new_loc_mods;
  int phospho_idx;
  for (std::vector<LocationMod>::const_iterator it = loc_mods.begin();
       it != loc_mods.end(); ++it) {
    if (var_mods[it->mod_id].name() != "Phospho")
      new_loc_mods.push_back(*it);
    else
      phospho_idx = it->mod_id;
  }
  for (std::vector<int>::const_iterator it = one_site_comb.begin();
       it != one_site_comb.end(); ++it) {
    new_loc_mods.push_back(LocationMod(*it, phospho_idx));
  }
  loc_mods.clear();
  loc_mods.assign(new_loc_mods.begin(), new_loc_mods.end());
}

void Ascore::MapInputSiteCombinationsInAll() {
  for (std::vector<std::vector<int> >::iterator it =
           input_phospho_site_combinations_.begin();
       it != input_phospho_site_combinations_.end(); ++it) {
    for (std::vector<std::vector<int> >::size_type i = 0;
         i < all_phospho_site_combinations_.size(); ++i) {
      if (IsEqual(*it, all_phospho_site_combinations_[i])) {
        indexes_in_all_.push_back(i);
        break;
      }
    }
  }
}

void Ascore::InitAllPeptideScores(const Parameters& paras) {
  for (std::vector<Match>::const_iterator it = pep_spec_matches_.begin();
       it != pep_spec_matches_.end(); ++it) {
    std::vector<double> logprobs;
    for (int peak_depth = 1;
         peak_depth <= paras.preproc_parameters_[2]; ++peak_depth) {
      double prob = it->PeptideScoreForAscore(peak_depth,
                                              paras.preproc_parameters_[1]);
      logprobs.push_back(-log10(prob));
    }
    pep_scores_.push_back(logprobs);
  }
}

double Ascore::GetWeigthedAverage(const std::vector<double>& values,
                                  const std::vector<double>& weights) {
  double weighted_score = 0;
  double sum_weights = 0;
  for (std::vector<double>::size_type i = 0; i < values.size(); ++i) {
    sum_weights += weights[i];
    weighted_score += weights[i] * values[i];
  }
  return weighted_score / sum_weights;
}

void Ascore::CalculateAllWeigthedAveragePeptideScores() {
  double w[10] = {0.5, 0.75, 1, 1, 1, 1, 0.75, 0.5, 0.25, 0.25};
  std::vector<double> weights(w, w + sizeof(w) / sizeof(double));
  for (std::vector<std::vector<double> >::const_iterator it = pep_scores_.begin();
       it != pep_scores_.end(); ++it) {
    double ws = GetWeigthedAverage(*it, weights);
    weigthed_average_pep_scores_.push_back(ws);
  }
}

// return: a vector of pair<competing site, index of the site combination>
std::vector<std::pair<int, int> > Ascore::GetCompetingPhosphoSites(
    std::vector<int>& curr_site_comb) {
  std::vector<std::pair<int, int> > competing_sites;
  std::set<int> curr_site_set(curr_site_comb.begin(), curr_site_comb.end());
  for (std::vector<int>::const_iterator it_x = curr_site_comb.begin();
       it_x != curr_site_comb.end(); ++it_x) {
    for (std::vector<int>::const_iterator it_y = indexes_.begin();
         it_y != indexes_.end(); ++it_y) {
      std::set<int> ref_site_set(all_phospho_site_combinations_[*it_y].begin(),
                                 all_phospho_site_combinations_[*it_y].end());
      std::set<int> diff;
      std::set_difference(curr_site_set.begin(), curr_site_set.end(),
                          ref_site_set.begin(), ref_site_set.end(),
                          std::inserter(diff, diff.begin()));
      if (diff.size() == 1 && *(diff.begin()) == *it_x) {
        diff.clear();
        std::set_difference(ref_site_set.begin(), ref_site_set.end(),
                            curr_site_set.begin(), curr_site_set.end(),
                            std::inserter(diff, diff.begin()));
        std::pair<int, int> pr(*it_y, *(diff.begin()));
        competing_sites.push_back(pr);
        break;
      }
    }
  }
  return competing_sites;
}

void Ascore::SortWeightedAveragePeptideScores() {
  for (std::vector<double>::size_type i = 0;
       i < weigthed_average_pep_scores_.size(); ++i) {
    indexes_.push_back(i);
  }
  std::sort(indexes_.begin(), indexes_.end(),
            GreaterScore(weigthed_average_pep_scores_));
}

int Ascore::SelectPeakDepth(const std::vector<double>& log_prob1,
                            const std::vector<double>& log_prob2) {
  double max_diff = 0;
  int depth = 1;
  for (std::vector<double>::size_type i = 0; i < log_prob1.size(); ++i) {
    double diff = abs(log_prob1[i] - log_prob2[i]);
    if (diff > max_diff) {
      max_diff = diff;
      depth = i + 1;
    }
  }
  return depth;
}

std::vector<double> Ascore::CalculatePhosphoSiteProbabilities(
    std::vector<int>& curr_phospho_site_comb, int index_in_all) {
  std::vector<std::pair<int, int> > competing_sites;
  competing_sites = GetCompetingPhosphoSites(curr_phospho_site_comb);
  std::vector<double> probs;
  for (std::vector<int>::size_type i = 0;
       i < curr_phospho_site_comb.size(); ++i) {
    int peak_depth = SelectPeakDepth(pep_scores_[index_in_all],
                                     pep_scores_[competing_sites[i].second]);
    double prob = pep_spec_matches_[index_in_all].PeptideScoreForAscore(
        peak_depth, curr_phospho_site_comb[i], competing_sites[i].first);
    probs.push_back(prob);
  }
  return probs;
}

void Ascore::InitPhosphoSiteProbabilities() {
  for (std::vector<std::vector<int> >::size_type i = 0;
       i < input_phospho_site_combinations_.size(); ++i) {
    std::vector<double> site_probs = CalculatePhosphoSiteProbabilities(
        input_phospho_site_combinations_[i], indexes_in_all_[i]);
    phospho_site_probabilities_.push_back(site_probs);
  }
}

void Ascore::InitLocalPepScores() {
  for (std::vector<std::vector<double> >::const_iterator it =
           phospho_site_probabilities_.begin();
       it != phospho_site_probabilities_.end(); ++it) {
    double product = 1;
    for (std::vector<double>::const_iterator it_y = it->begin();
         it_y != it->end(); ++it_y) {
      product *= (*it_y);
    }
    local_pep_scores_.push_back(1 - product);
  }
}

template<typename T>
bool Ascore::IsEqual(std::vector<T>& v1, std::vector<T>& v2) {
  std::sort(v1.begin(), v1.end());
  std::sort(v2.begin(), v2.end());
  return v1 == v2;
}
