#ifndef ION_H_
#define ION_H_

#include "ion_type.h"

class Ion {
 public:
  Ion();
  Ion(const IonType& ion_type, int cleavage_site,
      double m_over_z, double intensity = 0);
  Ion(const IonType& ion_type, int cleavage_site, int cleavage_site2,
      double m_over_z, double intensity = 0);
  Ion(const Ion& ion);
  virtual ~Ion();

  Ion& operator=(const Ion& ion);

  IonType ion_type() const { return ion_type_; }
  char cleavage_site(bool second = false) const {
    return second ? cleavage_site2_ : cleavage_site_;
  }
  double m_over_z() const { return m_over_z_; }
  double intensity() const { return intensity_; }

  void set_ion_type(const IonType& ion_type) { ion_type_ = ion_type; }
  void set_cleavage_site(int clvg_site, bool second = false) {
    if (second) cleavage_site2_ = clvg_site;
    else cleavage_site_ = clvg_site;
  }
  void set_m_over_z(double mz) { m_over_z_ = mz; }
  void set_intensity(double intensity) { intensity_ = intensity; }

  std::string GetIonLabel() const;

 protected:
  IonType ion_type_;
  char cleavage_site_;
  char cleavage_site2_;
  double m_over_z_;
  double intensity_;
};

#endif // ION_H_
