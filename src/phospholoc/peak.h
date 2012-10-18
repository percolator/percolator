#ifndef PEAK_H_
#define PEAK_H_

#include "misc.h"
#include <vector>

class Ion;

class Peak {
 public:
  Peak();
  Peak(double m_over_z, double intensity, int rank = 0, int window_id = 0);
  Peak(const Peak& peak);
  virtual ~Peak();

  Peak& operator=(const Peak &peak);

  double m_over_z() const { return m_over_z_; }
  double intensity() const { return intensity_; }
  int rank() const { return rank_; }
  int window_id() const { return window_id_; }

  void set_m_over_z(double m_over_z) { m_over_z_ = m_over_z; }
  void set_intensity(double intensity) { intensity_ = intensity; }
  void set_rank(int rank) { rank_ = rank; }
  void set_window_id(int window_id) { window_id_ = window_id; }

 protected:
  double m_over_z_;
  double intensity_;
  // for preprocess method: window top-N
  int rank_;  // 'rank' of this peak intensity in 'window_id'
  int window_id_;
};

#endif // PEAK_H_
