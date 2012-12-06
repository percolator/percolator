#ifndef PHOSPHORS_H_
#define PHOSPHORS_H_

#include "match.h"
#include <cmath>

class PhosphoRS {
 public:
  PhosphoRS();
  PhosphoRS(const Spectrum& spec, const std::string& pep_seq,
            std::vector<std::vector<LocationMod> >& input_var_mod_combs,
            const Parameters& paras);
  virtual ~PhosphoRS();

  double GetPeptideScore(std::vector<int>& phospho_locs);
  double GetSequenceProbability(std::vector<int>& phospho_locs);
  double GetPhosphoSiteScore(std::vector<int>& phospho_locs);

 private:
  void InitPhosphoRS(const Spectrum& spec,
                     const std::string pep_seq,
                     std::vector<std::vector<LocationMod> >& input_var_mod_combs,
                     const Parameters& paras);
  // begin: this part should be refractored into a new class
  std::vector<int> GetPhosphoedLocations(
      const std::vector<LocationMod>& loc_mods);
  std::vector<int> GetPotentialPhosphoSites(const std::string& pep_seq);
  void SetAllPhosphoSiteCombinations(std::vector<int>& potential_phospho_sites,
                                     int actual_mod_num);
  void ResetPhosphoedSites(const std::vector<int>& one_site_comb,
                           std::vector<LocationMod>& loc_mods);
  void MapInputSiteCombinationsInAll();
  // end
  void InitScoreRanksCurrentWindow();
  void SortScoresAtPeakDepth(int peak_depth);
  double DeltaScoreAtPeakDepthAndRank(int peak_depth, int score_rank);
  int PeakDepthAtMaxScore(std::vector<double>& scores_all_depths);
  bool AreEqualScores(double score1, double score2, double tol);
  template<typename T>
  bool AreEqualVectors(std::vector<T>& v1, std::vector<T>& v2);
  std::vector<int> PeakDepthsAtMaxDeltaScore(
      int rank_to_compare, std::vector<int>& candidate_peak_depths);
  int SelectOptimalPeakDepthInOneWindow();
  void SelectOptimalPeakDepthsForAllWindows();
  void InitPeptideScoresAndSequenceProbs();
  void InitPhosphoSiteProbabilities();

  struct GreaterScore {
    GreaterScore(std::vector<double>& scr)
        : scores(scr) {
    }
    bool operator() (const int& a, const int& b) const {
      return scores[a] > scores[b];
    }
    std::vector<double>& scores;
  };


 protected:
 private:
  // variables with '<--' should be refactored into a new class
  // all phospho-PSMs corresponding to all possible phospho-site combinations
  std::vector<Match> pep_spec_matches_; // <--
  // all possible phospho-site combinations
  std::vector<std::vector<int> > all_phospho_site_combinations_; // <--
  // optimal peak depth for each m/z window
  std::vector<int> optimal_peak_depths_; // size=#windows
  // 2D array for current window containing local peptide scores for each combination
  // #rows=#site_combinations, #columns=#peak_depths
  std::vector<std::vector<double> > peptide_scores_current_window_;
  // index of a score ranked i-th in peptide_scores_current_window_
  // is index_scores_current_window_[i]
  std::vector<int> index_scores_current_window_;
  // cumulative binomial probs -10*log(Prob), Prob=SUM_[k..n](C(n,k)*p^k*(1-p)^(n-k))
  // p=N_peaks*d/w, d is specified ion tolerance, w is full mass range of a spectrum
  std::vector<double> peptide_scores_; // size=#site_combinations
  // (1/Prob) normalized by sum of all (1/Prob)
  std::vector<double> sequence_probs_;
  // for a specified site, sum over all sequence_probs_ that containing this site
  std::map<int, double> phospho_site_probs_; // <site_index, site_prob>

  // phospho-site combinations actually occurred in search result
  std::vector<std::vector<int> > input_phospho_site_combinations_; // <--
  // index of the above phospho-sites in 'all_phospho_site_combinations'
  std::vector<int> indexes_in_all_; // <--

  static double min_peak_depth_;
  static double max_peak_depth_;

  DISALLOW_COPY_AND_ASSIGN(PhosphoRS);
};

#endif // PHOSPHORS_H_
