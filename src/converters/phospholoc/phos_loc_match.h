#ifndef MATCH_H_
#define MATCH_H_

#include "phos_loc_peptide.h"
#include "phos_loc_ion_series.h"
#include "phos_loc_spectrum.h"
#include "phos_loc_parameters.h"
#include "phos_loc_utils.h"
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>

namespace phos_loc {

class Match {
 public:
  Match();
  Match(const Spectrum& spec, const IonSeries& ions,
        const Tolerance& tol = Tolerance());
  Match(const Match& match);
  virtual ~Match();

  const Spectrum& spectrum() const { return spectrum_; }
  const IonSeries& ion_series() const { return ion_series_; }
  const std::map<int, std::vector<int> >& peak_to_ions() const {
    return peak_to_ions_;
  }
  const std::map<int, std::vector<int> >& ion_to_peaks() const {
    return ion_to_peaks_;
  }
  const std::map<IonType, std::vector<std::vector<int> > > ion_match_matrix() {
    return ion_match_matrix_;
  }

  std::vector<int> PeaksForSingleIon(int start_peak_idx, double ion_mz);
  void FindAllPeakIonMatches();
  int NumAllPeaks(std::vector<int>& peak_depths) const;
  double ProbabilityOfAPeakIonMatch(int peak_depth, double window_width) const;
  double ProbabilityOfAPeakIonMatch(std::vector<int>& peak_depths) const;
  double PeptideScoreForAscore(int peak_depth, int lower_site,
                               int upper_site, double window_width = 100) const;
  double PeptideScoreForAscore(int peak_depth, double window_width = 100) const;
  double PeptideScoreForPhosphoRS(int peak_depth, double lower_bnd,
                                  double upper_bnd) const;
  double PeptideScoreForPhosphoRS(std::vector<int>& optimal_peak_depths) const;
  void Print();

 private:
  void InitIonMatchMatrix();
  bool IsMatched(double peak_mz, double ion_mz);
  int IndexOfHighestPeak(std::vector<int>& peak_indexes);
  void InsertMatchedPeaksForIon(int ion_idx, std::vector<int>& peak_indexes,
                                bool only_keep_highest = true);
  void InsertMatchedIonForPeaks(int ion_idx, std::vector<int>& peak_indexes,
                                bool only_keep_highest = true);
  int GetNumMatchedIons(int peak_depth, int lower_site, int upper_site) const;
  int GetNumMatchedIons(int peak_depth) const;
  int GetNumMatchedIons(int peak_depth, double lower_bnd, double upper_bnd) const;
  int GetNumPredictedIons(int lower_site, int upper_site) const;
  int GetNumPredictedIons() const;
  int GetNumPredictedIons(double lower_bnd, double upper_bnd) const;

  Spectrum spectrum_;
  IonSeries ion_series_;
  Tolerance frag_tol_;
  // each vector<int> corresponds to the indexes of matched ions
  // in ion_series.sorted_ions_, first int is the index of peaks
  std::map<int, std::vector<int> > peak_to_ions_;
  // each vector<int> corresponds to the indexes of matched peaks
  // in spectrum.peaks_, first int is the index of sorted_ions
  std::map<int, std::vector<int> > ion_to_peaks_;
  // corresponding to IonSeries.ion_matrix_
  std::map<IonType, std::vector<std::vector<int> > > ion_match_matrix_;
};

} // namespace phos_loc

#endif // MATCH_H_
