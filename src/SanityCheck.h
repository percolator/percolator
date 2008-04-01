#ifndef SANITYCHECK_H_
#define SANITYCHECK_H_
class Scores;
class Normalizer;

class SanityCheck
{
public:
  SanityCheck();
  virtual ~SanityCheck();
    
  void readWeights(istream & weightStream, double * w);
  int getInitDirection(Scores * testset,Scores * trainset, Normalizer * pNorm,double *w,double test_fdr);
  virtual bool validateDirection(double* w);
  void resetDirection(double* w);

  static void setInitWeightFN(string fn) {initWeightFN=fn;}
  static bool overRule;
protected:
  virtual void getDefaultDirection(double *w);
  int initPositives;
  double fdr;
  static string initWeightFN;
  Scores *pTestset, *pTrainset;
};

#endif /*SANITYCHECK_H_*/
