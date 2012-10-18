#ifndef SPECTRUM_H_
#define SPECTRUM_H_

#include "misc.h"
#include "peak.h"
#include <string>
#include <vector>
#include <cstdio>
#include <cmath>
#include <algorithm>


class Spectrum {
 public:
  Spectrum();
  Spectrum(const std::string& dta_file_path, MassType precursor_mass_t = MONO);
  Spectrum(const Spectrum& spec);
  virtual ~Spectrum();

  Spectrum& operator=(const Spectrum& spec);
  Peak& operator[](const int& idx);

  Precursor precursor() const { return precursor_; }
  int scan_num(bool second = false) const {
    return second ? scan_num2_ : scan_num1_;
  }
  char ms_level() const { return ms_level_; }
  double retention_time() const { return retention_time_; }
  double retention_time_apex() const { return retention_time_apex_; }
  const std::vector<Peak>& peaks() const { return peaks_; }
  Peak base_peak() const { return base_peak_; }
  double total_ion_current() const { return total_ion_current_; }
  ActivationType activation_type() const { return activation_type_; }
  SpectraFormat format() const { return format_; }
  std::string file_path() const { return file_path_; }

  void set_precursor(const Precursor& prec) { precursor_ = prec; }
  void set_scan_num(int scan_num, bool second) {
    if (second) scan_num2_ = scan_num;
    else scan_num1_ = scan_num;
  }
  void set_ms_level(char ms_level) { ms_level_ = ms_level; }
  void set_retention_time(double rt) { retention_time_ = rt; }
  void set_retention_time_apex(double rt) { retention_time_apex_ = rt; }
  void set_peaks(const std::vector<Peak>& peaks) {
    peaks_.assign(peaks.begin(), peaks.end());
  }
  void set_base_peak(const Peak& peak) { base_peak_ = peak; }
  void set_total_ion_current(double tic) { total_ion_current_ = tic; }

  void Print();
  void Preprocess(PreprocessMethod prep_method, double* prep_paras,
                  Tolerance tolerance);

 private:
  void ReadDtaFile(const std::string& dta_file_path,
                   MassType precursor_mass_t = AVERAGE);
  bool GetLine(char* line, FILE* file);
  static bool LessPeakMZ(const Peak& peak1, const Peak& peak2);
  static bool HigherPeakInternsity(const Peak& peak1, const Peak& peak2);
  double CalcLossMZ(double neutral_loss);
  void RemovePhosphoLossPeaks(Tolerance tolerance);
  std::vector<int> PeaksForSingleIon(int start_peak_idx, double ion_mz,
                                     Tolerance tol);
  bool IsMatched(double peak_mz, double ion_mz, Tolerance tol);
  void SplitPeaksIntoWindows(double start_mz, double win_width,
                             int max_num_peaks_per_win);
  void RankAndFilterPeaksInSingleWindow(std::vector<Peak>& win_peaks,
                                        int win_id, int max_num_peaks_per_win);

 protected:
  Precursor precursor_;
  int scan_num1_;
  int scan_num2_;
  char ms_level_;
  double retention_time_;  // in second
  double retention_time_apex_;
  std::vector<Peak> peaks_;
  Peak base_peak_;
  double total_ion_current_;
  ActivationType activation_type_;
  SpectraFormat format_;
  std::string file_path_;
};

#endif // SPECTRUM_H_
