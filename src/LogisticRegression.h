#ifndef LOGISTICREGRESSION_H_
#define LOGISTICREGRESSION_H_
#include<vector>
using namespace std;
#include "BaseSpline.h"

class LogisticRegression : public BaseSpline
{
public:
  LogisticRegression();
  virtual ~LogisticRegression();
  virtual double predict(const double x);
  void setData(const vector<double>& xx, const vector<unsigned int>& yy, const vector<unsigned int>& mm) {BaseSpline::setData(xx);y=yy;m=mm;}
  void setCutOff(double pi0);
protected:
  virtual void calcPZW();
  virtual void initg();
  vector<unsigned int> y,m;
  double cutOff,pi0;
  Vec p;
};

#endif /*LOGISTICREGRESSION_H_*/
