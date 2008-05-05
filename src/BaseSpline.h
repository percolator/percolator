#ifndef BASESPLINE_H_
#define BASESPLINE_H_

#include "Transform.h"
#include "ArrayLibrary.h"

class BaseSpline
{
public:
//  BaseSpline() : pTransf(NULL) {;}
//  virtual ~BaseSpline() {if (pTransf) delete pTransf;}
  BaseSpline() {;}
  virtual ~BaseSpline() {;}
  double splineEval(double xx);
  static double convergeDelta;
  void iterativeReweightedLeastSquares();
  void predict(const vector<double>& x, vector<double>& predict);
  void setData(const vector<double>& x);
  virtual double predict(double xx) {return splineEval(xx);}
protected:
  virtual void calcPZW() {;}
  virtual void initg() {int n=x.size();g.resize(n);gnew.resize(n);w.resize(n);z.resize(n,0.5);gamma.resize(n-2);}
  void initiateQR();
  double crossValidation(double alpha);
  pair<double,double> alphaLinearSearch(double min_p,double max_p, double p1, double p2, double cv1, double cv2);

  Transform transf;

  PackedMatrix Q,Qt,R;
  Vec gnew,w,z,dx;
  Vec g,gamma;
  vector<double> x;
};

#endif /*BASESPLINE_H_*/
