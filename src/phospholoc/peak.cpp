#include "peak.h"


Peak::Peak()
    : m_over_z_(0),
      intensity_(0),
      rank_(0),
      window_id_(0) {
}

Peak::Peak(double m_over_z, double intensity, int rank, int window_id)
    : m_over_z_(m_over_z),
      intensity_(intensity),
      rank_(rank),
      window_id_(window_id) {
}

Peak::Peak(const Peak& peak)
    : m_over_z_(peak.m_over_z_),
      intensity_(peak.intensity_),
      rank_(peak.rank_),
      window_id_(peak.window_id_) {
}

Peak::~Peak() {
}

Peak& Peak::operator=(const Peak& peak) {
  if (this != &peak) {
    m_over_z_ = peak.m_over_z_;
    intensity_ = peak.intensity_;
    rank_ = peak.rank_;
    window_id_ = peak.window_id_;
  }
  return *this;
}
