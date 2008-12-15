#include "PSMDescription.h"
#include "DescriptionOfCorrect.h"

PSMDescription::PSMDescription()
{
}

PSMDescription::~PSMDescription()
{
}

void PSMDescription::calcRegressionFeature() {
  string::size_type pos1 = peptide.find('.');
  string::size_type pos2 = peptide.find('.',++pos1);
  string pep = peptide.substr(pos1,pos2-pos1);
  pI = DescriptionOfCorrect::isoElectricPoint(pep);
  if (retentionFeatures) {
    DescriptionOfCorrect::fillFeaturesAllIndex(pep, retentionFeatures);
  }
//  cout <<  peptide << " " << pep << " " << retentionFeatures[0] << endl;
}

