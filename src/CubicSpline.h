#ifndef CUBICSPLINE_H_
#define CUBICSPLINE_H_
#include <math.h>

class CubicSpline
{
 public:
  CubicSpline();
  CubicSpline(vector<double> &x, vector<double> &y, bool l=false);
  virtual ~CubicSpline();
  void setData(vector<double> &xx, vector<double> &yy, bool l=false);
  void removeDuplicates();
  void calcDeriv();
  double interpolate(double xx);
  double linearInterpolate(double xx);
  static double logit(double p) {return log(p/(1-p));}
  static double invlogit(double l) {double el=exp(l);return el/(1+el);}
 protected:
  vector<double> y,x,d2y;
  double yp0,ypn; 
  bool derivCalculated,natural0,naturaln,logged; 
};

#endif /*CUBICSPLINE_H_*/
