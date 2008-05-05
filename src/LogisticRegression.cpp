#include "ArrayLibrary.h"
#include "LogisticRegression.h"

LogisticRegression::LogisticRegression()
{
}

LogisticRegression::~LogisticRegression()
{
}

void invlogit(double & out,double in) {
  double e = exp(in);
  out = e/(1+e);
}

double invlogit(double in) {
  double e = exp(in);
  return e/(1+e);
}

double logit(double p) {
  return log(p/(1-p));
}

double LogisticRegression::predict(const double xx) {
  if (xx<cutOff)
    return 1.0;
  return min(1.0,pi0*exp(splineEval(xx)));  
}


void LogisticRegression::calcPZW() {
  for (int ix=z.size();ix--;) {
    double e = exp(g[ix]);
    p[ix] = min(max(e/(1+e),Numerical::epsilon),1-Numerical::epsilon);
    w[ix] = m[ix]*p[ix]*(1-p[ix]);
    z[ix] = g[ix] + (((double)y[ix])-p[ix]*((double)m[ix]))/w[ix];
  }
}

void LogisticRegression::initg() {
  BaseSpline::initg();
  int n=x.size();
  p.resize(n);  
  for (int ix=g.size();ix--;) {
    double p = (y[ix]+0.1)/(m[ix]+0.1);
    gnew[ix] = log(p/(1-p));
  }
}

void LogisticRegression::setCutOff(double pi_0) {
  pi0 = pi_0;
  int ix=g.size();
  for (;ix--;) {
    if (pi0*exp(g[ix])>1)
      break; 
  }
  if (ix<0) ix=0;
  cutOff = x[ix];
}


