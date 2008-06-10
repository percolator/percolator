#ifndef DESCRIPTIONOFCORRECT_H_
#define DESCRIPTIONOFCORRECT_H_
#include<string>
using namespace std;

class DescriptionOfCorrect
{
public:
  DescriptionOfCorrect();
  virtual ~DescriptionOfCorrect();
protected:
  double isoElectricPoint(string peptide);
  double kyteDolittle(string peptide);
  static float krokhin_index['Z'-'A'+1],hessa_index['Z'-'A'+1],kytedoolittle_index['Z'-'A'+1];
};

#endif /*DESCRIPTIONOFCORRECT_H_*/
