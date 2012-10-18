#include "spectrum.h"

Spectrum::Spectrum()
    : scan_num1_(0),
      scan_num2_(0),
      ms_level_(2),
      retention_time_(0),
      retention_time_apex_(0),
      total_ion_current_(0),
      activation_type_(CID),
      format_(DTA),
      file_path_("") {
}

Spectrum::Spectrum(const std::string& dta_file_path, MassType precursor_mass_t)
    : file_path_(dta_file_path) {
  ReadDtaFile(dta_file_path, precursor_mass_t);
}

Spectrum::Spectrum(const Spectrum& spec)
    : precursor_(spec.precursor_),
      scan_num1_(spec.scan_num1_),
      scan_num2_(spec.scan_num2_),
      ms_level_(spec.ms_level_),
      retention_time_(spec.retention_time_),
      retention_time_apex_(spec.retention_time_apex_),
      peaks_(spec.peaks_),
      base_peak_(spec.base_peak_),
      total_ion_current_(spec.total_ion_current_),
      activation_type_(spec.activation_type_),
      format_(spec.format_),
      file_path_(spec.file_path_) {
}

Spectrum::~Spectrum() {
}

Spectrum& Spectrum::operator=(const Spectrum& spec) {
  if (this != &spec) {
    precursor_ = spec.precursor_;
    scan_num1_ = spec.scan_num1_;
    scan_num2_ = spec.scan_num2_;
    ms_level_ = spec.ms_level_;
    retention_time_ = spec.retention_time_;
    retention_time_apex_ = spec.retention_time_apex_;
    peaks_ = spec.peaks_;
    base_peak_ = spec.base_peak_;
    total_ion_current_ = spec.total_ion_current_;
    activation_type_ = spec.activation_type_;
    format_ = spec.format_;
    file_path_ = spec.file_path_;
  }
  return *this;
}

Peak& Spectrum::operator[](const int& idx) {
  return peaks_[idx];
}

void Spectrum::Print() {
  fprintf(stdout, ">Peaks:\n");
  fprintf(stdout, "Precursor chg:%d mz:%f inten:%f\n", precursor_.charge,
          precursor_.m_over_z, precursor_.intensity);
  for (std::vector<Peak>::const_iterator it = peaks_.begin();
       it != peaks_.end(); ++it) {
    fprintf(stdout, "%f %f %d %d\n", it->m_over_z(), it->intensity(),
            it->rank(), it->window_id());
  }
}

void Spectrum::Preprocess(PreprocessMethod prep_method, double* prep_paras,
                          Tolerance tolerance) {
  if (TOPN_WINDOW == prep_method) {
    RemovePhosphoLossPeaks(tolerance);
    SplitPeaksIntoWindows(prep_paras[0], prep_paras[1], (int) prep_paras[2]);
  }
}

void Spectrum::ReadDtaFile(const std::string& dta_file_path,
                           MassType precursor_mass_t) {
  FILE* dta_file;
  dta_file = fopen(dta_file_path.c_str(), "r");
  if (NULL == dta_file) {
    fprintf(stderr, "- Cannot open file: %s.\n", dta_file_path.c_str());
    return;
  }

  char line[MAX_LENGTH_SPECTRUM_LINE];
  int charge;
  double mhplus, mz, inten;

  // read in precursor
  if (!GetLine(line, dta_file)) {
    fprintf(stderr, "- Error occurred or empty file: %s.\n", dta_file_path.c_str());
    return;
  }
  if (2 != sscanf(line, "%lf %d", &mhplus, &charge)) {
    fprintf(stderr, "- Error in data format: %s.\n", dta_file_path.c_str());
    return;
  }
  if (MONO == precursor_mass_t)
    mz = (mhplus + (charge - 1) * PROTON) / charge;
  else
    mz = (mhplus + (charge - 1) * HYDROGEN[AVERAGE]) / charge;
  precursor_ = Precursor(charge, mz, precursor_mass_t, 0);

  // read in peaks
  while (!feof(dta_file)) {
    if (!GetLine(line, dta_file)) {
      return;
    }
    if (2 != sscanf(line, "%lf %lf", &mz, &inten)) {
      fprintf(stderr, "- Error in data format: %s.\n", dta_file_path.c_str());
      return;
    }
    peaks_.push_back(Peak(mz, inten));
  }

  fclose(dta_file);
}

bool Spectrum::GetLine(char* line, FILE* file) {
  while (!feof(file)) {
    if (NULL == fgets(line, MAX_LENGTH_SPECTRUM_LINE, file))
      return false;
    if (strspn(line, " \t\n\r") == strlen(line))
      continue;
    return true;
  }
  return false;
}

bool Spectrum::LessPeakMZ(const Peak& peak1, const Peak& peak2) {
  return peak1.m_over_z() < peak2.m_over_z();
}

bool Spectrum::HigherPeakInternsity(const Peak& peak1, const Peak& peak2) {
  return peak1.intensity() > peak2.intensity();
}

double Spectrum::CalcLossMZ(double neutral_loss) {
  return (precursor_.m_over_z - (neutral_loss / precursor_.charge));
}

void Spectrum::RemovePhosphoLossPeaks(Tolerance tol) {
  std::vector<double> neutral_peak_mz;
  MassType mt = tol.mass_type;
  neutral_peak_mz.push_back(CalcLossMZ(PHOSPHOLOSS[mt] + H2O[mt]));
  neutral_peak_mz.push_back(CalcLossMZ(PHOSPHOLOSS[mt] + H2O[mt] - C13ISOTOPE));
  neutral_peak_mz.push_back(CalcLossMZ(PHOSPHOLOSS[mt] + H2O[mt] - 2*C13ISOTOPE));
  neutral_peak_mz.push_back(CalcLossMZ(PHOSPHOLOSS[mt]));
  neutral_peak_mz.push_back(CalcLossMZ(PHOSPHOLOSS[mt] - C13ISOTOPE));
  neutral_peak_mz.push_back(CalcLossMZ(PHOSPHOLOSS[mt] - 2*C13ISOTOPE));
  neutral_peak_mz.push_back(CalcLossMZ(H2O[mt] + H2O[mt]));
  neutral_peak_mz.push_back(CalcLossMZ(H2O[mt]));

  int curr_peak_idx = 0;
  std::vector<bool> remove_pos(peaks_.size(), false);
  for (std::vector<double>::iterator it = neutral_peak_mz.begin();
       it != neutral_peak_mz.end(); ++it) {
    std::vector<int> matched_peak_idx = PeaksForSingleIon(curr_peak_idx,
                                                          *it, tol);
    for (std::vector<int>::iterator it_idx = matched_peak_idx.begin();
         it_idx != matched_peak_idx.end(); ++it_idx)
      remove_pos[*it_idx] = true;
    if (!matched_peak_idx.empty())
      curr_peak_idx = matched_peak_idx[0];
  }

  std::vector<Peak> new_peaks;
  for (std::vector<int>::size_type i = 0; i < peaks_.size(); ++i) {
    if (false == remove_pos[i])
      new_peaks.push_back(peaks_[i]);
    else
      fprintf(stderr, "DEL: %f %f\n", peaks_[i].m_over_z(),
              peaks_[i].intensity());
  }
  peaks_.clear();
  peaks_.assign(new_peaks.begin(), new_peaks.end());
}

std::vector<int> Spectrum::PeaksForSingleIon(int start_peak_idx,
                                             double ion_mz,
                                             Tolerance tol) {
  std::vector<int> matched_peak_indexes;
  for (std::vector<Peak>::size_type i = start_peak_idx;
       i < peaks_.size(); ++i) {
    double peak_mz = peaks_[i].m_over_z();
    if (IsMatched(peak_mz, ion_mz, tol))
      matched_peak_indexes.push_back(i);
    else if (peak_mz > ion_mz)
      break;
  }
  return matched_peak_indexes;
}

bool Spectrum::IsMatched(double peak_mz, double ion_mz, Tolerance tol) {
  return fabs(peak_mz - ion_mz) <= tol.value;
}

void Spectrum::SplitPeaksIntoWindows(double start_mz, double win_width,
                                     int max_num_peaks_per_win) {
  if (-1 == start_mz)
    start_mz = peaks_.begin()->m_over_z();

  int start_peak_idx = 0;
  while (peaks_[start_peak_idx].m_over_z() < start_mz)
    ++start_peak_idx;

  int curr_win_id = 0;
  std::vector<Peak> new_peaks;
  std::vector<Peak> win_peaks;
  for (std::vector<Peak>::const_iterator it = peaks_.begin() + start_peak_idx;
       it != peaks_.end(); ++it) {
    if (it->m_over_z() < (curr_win_id + 1) * win_width + start_mz) {
      win_peaks.push_back(*it);
    }
    else {
      RankAndFilterPeaksInSingleWindow(win_peaks, curr_win_id, max_num_peaks_per_win);
      new_peaks.insert(new_peaks.end(), win_peaks.begin(), win_peaks.end());
      win_peaks.clear();
      ++curr_win_id;
      win_peaks.push_back(*it);
    }
  }
  RankAndFilterPeaksInSingleWindow(win_peaks, curr_win_id, max_num_peaks_per_win);
  new_peaks.insert(new_peaks.end(), win_peaks.begin(), win_peaks.end());
  win_peaks.clear();
  peaks_.clear();
  peaks_.assign(new_peaks.begin(), new_peaks.end());
}

void Spectrum::RankAndFilterPeaksInSingleWindow(std::vector<Peak>& win_peaks,
                                                int win_id,
                                                int max_num_peaks_per_win) {
  if (win_peaks.empty()) return;
  std::sort(win_peaks.begin(), win_peaks.end(), HigherPeakInternsity);
  int max_num = std::min((int)win_peaks.size(), max_num_peaks_per_win);
  for (int i = 0; i < max_num; ++i) {
    win_peaks[i].set_rank(i);
    win_peaks[i].set_window_id(win_id);
  }
  win_peaks.erase(win_peaks.begin()+max_num, win_peaks.end());
  std::sort(win_peaks.begin(), win_peaks.end(), LessPeakMZ);
}
