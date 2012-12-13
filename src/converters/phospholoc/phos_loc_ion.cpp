#include "phos_loc_ion.h"

namespace phos_loc {

Ion::Ion()
    : cleavage_site_(0),
      cleavage_site2_(0),
      m_over_z_(0),
      intensity_(0) {
}

Ion::Ion(const IonType& ion_type, int cleavage_site,
         double m_over_z, double intensity)
    : ion_type_(ion_type),
      cleavage_site_(cleavage_site),
      cleavage_site2_(0),
      m_over_z_(m_over_z),
      intensity_(intensity) {
}

Ion::Ion(const IonType& ion_type, int cleavage_site,int cleavage_site2,
         double m_over_z, double intensity)
    : ion_type_(ion_type),
      cleavage_site_(cleavage_site),
      cleavage_site2_(cleavage_site2),
      m_over_z_(m_over_z),
      intensity_(intensity) {
}

Ion::Ion(const Ion& ion)
    : ion_type_(ion.ion_type_),
      cleavage_site_(ion.cleavage_site_),
      cleavage_site2_(ion.cleavage_site2_),
      m_over_z_(ion.m_over_z_),
      intensity_(ion.intensity_) {
}

Ion::~Ion() {
}

Ion& Ion::operator=(const Ion& ion) {
  if (this != &ion) {
    ion_type_ = ion.ion_type_;
    cleavage_site_ = ion.cleavage_site_;
    cleavage_site2_ = ion.cleavage_site2_;
    m_over_z_ = ion.m_over_z_;
  }
  return *this;
}

std::string Ion::GetIonLabel() const {
  char s[256];
  sprintf(s, "%s_%d_%d=%f", ion_type_.ion_name().c_str(),
          cleavage_site_, ion_type_.charge(), m_over_z_);
  std::string label(s);
  return label;
}

} // namespace phos_loc
