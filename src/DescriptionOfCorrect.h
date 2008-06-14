#ifndef DESCRIPTIONOFCORRECT_H_
#define DESCRIPTIONOFCORRECT_H_
#include<string>
using namespace std;

class DescriptionOfCorrect
{
public:
  DescriptionOfCorrect();
  virtual ~DescriptionOfCorrect();
  void fillFeaturesAllIndex(const string& peptide, double *features);
  double isoElectricPoint(const string& peptide);
protected:
  inline double indexSum(const float *index, const string& peptide);
  inline double indexN(const float *index, const string& peptide);
  inline double indexC(const float *index, const string& peptide);
  inline double indexNC(const float *index, const string& peptide);
  double* fillAAFeatures(const string& pep, double *feat);
  double* fillFeaturesIndex(const string& peptide, const float *index, double *features);
  double kyteDolittle(string peptide);
  static float krokhin_index['Z'-'A'+1],hessa_index['Z'-'A'+1],kytedoolittle_index['Z'-'A'+1];
  static string aaAlphabet,isoAlphabet;
  static float pKiso[7]; 
  static float pKN,pKC;
};

#endif /*DESCRIPTIONOFCORRECT_H_*/
