#ifndef POSTERIORESTIMATOR_H_
#define POSTERIORESTIMATOR_H_
#include<vector>
#include<utility>
using namespace std;
#include "LogisticRegression.h"

class PosteriorEstimator
{
public:
  PosteriorEstimator();
  virtual ~PosteriorEstimator();
  bool parseOptions(int argc, char **argv);   
  void run();
  void estimate( vector<pair<double,bool> >& combined, LogisticRegression& lr, double pi0);
  void getPValues(const vector<pair<double,bool> >& combined, vector<double>& p);
  void getQValues(double pi0,const vector<pair<double,bool> >& combined, vector<double>& q);
  double estimatePi0(vector<pair<double,bool> >& combined, const unsigned int numBoot=100);
protected:
  void finishStandalone(vector<pair<double,bool> >& combined, LogisticRegression& lr, double pi0);
  void binData(const vector<pair<double,bool> >& combined, vector<double>& medians, 
               vector<unsigned int>& negatives, vector<unsigned int>& sizes);
  string targetFile,decoyFile;

};

#endif /*POSTERIORESTIMATOR_H_*/
