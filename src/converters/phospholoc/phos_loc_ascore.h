#ifndef ASCORE_H_
#define ASCORE_H_

#include "phos_loc_match.h"
#include <set>

namespace phos_loc {

class Ascore {
 public:
  Ascore();
  Ascore(const Spectrum& spec, const std::string& pep_seq,
         std::vector<std::vector<LocationMod> >& input_var_mod_combs,
         const Parameters& paras);
  virtual ~Ascore();

  double GetPeptideScore(std::vector<int>& phospho_locs);
  double GetPhosphoSiteScore(std::vector<int>& phospho_locs);

 private:
  void InitAscore(const Spectrum& spec,
                  const std::string pep_seq,
                  std::vector<std::vector<LocationMod> >& input_var_mod_combs,
                  const Parameters& paras);
  std::vector<int> GetPhosphoedLocations(
      const std::vector<LocationMod>& loc_mods);
  std::vector<int> GetPotentialPhosphoSites(const std::string& pep_seq);
  void SetAllPhosphoSiteCombinations(std::vector<int>& potential_phospho_sites,
                                     int actual_mod_num);
  void ResetPhosphoedSites(const std::vector<int>& one_site_comb,
                           std::vector<LocationMod>& loc_mods);
  void MapInputSiteCombinationsInAll();
  void InitAllPeptideScores(const Parameters& paras);
  double GetWeigthedAverage(const std::vector<double>& values,
                            const std::vector<double>& weights);
  void CalculateAllWeigthedAveragePeptideScores();
  std::vector<std::pair<int, int> > GetCompetingPhosphoSites(
      std::vector<int>& curr_site_comb);
  void SortWeightedAveragePeptideScores();
  int SelectPeakDepth(const std::vector<double>& log_prob1,
                      const std::vector<double>& log_prob2);
  std::vector<double> CalculatePhosphoSiteProbabilities(
      std::vector<int>& curr_phospho_site_comb, int index_in_all);
  void InitPhosphoSiteProbabilities();
  void InitPhosphoSiteScores();
  template<typename T>
  bool AreEqualVectors(std::vector<T>& v1, std::vector<T>& v2);
  struct GreaterScore {
    GreaterScore(std::vector<double>& scr)
        : scores(scr) {
    }
    bool operator() (const int& a, const int& b) const {
      return scores[a] > scores[b];
    }
    std::vector<double>& scores;
  };

 private:
  // all phospho-PSMs corresponding to all possible phospho-site combinations
  std::vector<Match> pep_spec_matches_;
  // all possible phospho-site combinations
  std::vector<std::vector<int> > all_phospho_site_combinations_;
  // peptide scores [-log10(cumulative binomial prob)] for all phospho-site
  // combinations and peak depths, one row containing all peak depths for one comb
  std::vector<std::vector<double> > pep_scores_;
  std::vector<double> weigthed_average_pep_scores_;
  // for locking and sorting the above 4 vecotors
  std::vector<int> indexes_;

  // phospho-site combinations actually occurred in search result
  std::vector<std::vector<int> > input_phospho_site_combinations_;
  // index of the above phospho-sites in 'all_phospho_site_combinations'
  std::vector<int> indexes_in_all_;
  // probabilities of phospho-sites, only consider PSMs occurred in search result
  // corresponding to the above 'input_phospho_site_combinations_'
  std::vector<std::vector<double> > phospho_site_probabilities_;
  // local peptide scores derived from 'phospho_site_probabilities_'
  std::vector<double> local_pep_scores_;

  DISALLOW_COPY_AND_ASSIGN(Ascore);
};

} // namespace phos_loc

#endif // ASCORE_H_
