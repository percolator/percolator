#ifndef TRANSFORM_H_
#define TRANSFORM_H_

#include <math.h>

class Transform
{
public:
  Transform(double delt=0.0,bool dLt=false,bool dL=false) : delta(delt),doLogit(dLt),doLog(dL) {;}
  ~Transform() {;}
  double operator() (double xx) {
    if ((!doLogit) && (!doLog)) return xx;
    if (delta>0) {
      if (doLogit) xx *= (1-2*delta);
      xx += delta;
    }
    if (doLogit) return log(xx/(1.-xx));
    return log(xx);
  }
private:
  double delta;
  bool doLogit,doLog;
};

#endif /*TRANSFORM_H_*/
