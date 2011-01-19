/*******************************************************************************
 Copyright 2006-2009 Lukas Käll <lukas.kall@cbr.su.se>

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

 *******************************************************************************/

#ifndef BASESPLINE_H_
#define BASESPLINE_H_

#include <assert.h>
#include "Transform.h"
#include "Numerical.h"
#include "Vector.h"
#include "Matrix.h"

class BaseSpline {
  public:
    //  BaseSpline() : pTransf(NULL) {;}
    //  virtual ~BaseSpline() {if (pTransf) delete pTransf;}
    BaseSpline(){};
    virtual ~BaseSpline(){};
    double splineEval(double xx);
    static double convergeEpsilon;
    static double stepEpsilon;
    void iterativeReweightedLeastSquares();
    void predict(const vector<double>& x, vector<double>& predict);
    void setData(const vector<double>& x);
    double predict(double xx) {
      return splineEval(xx);
    }
    static void solveInPlace(Matrix& mat, Vector& res);
  protected:
    virtual void calcPZW() {
      ;
    }
    virtual void initg() {
      int n = x.size();
      g = Vector(n,0);
      gnew = Vector(n,0);
      w = Vector(n,0);
      z = Vector(n,0.5);
      gamma = Vector(n-2,0);
    }
    virtual void limitg() {
      ;
    }
    virtual void limitgamma() {
      ;
    }
    void initiateQR();
    double crossValidation(double alpha);
    pair<double, double> alphaLinearSearch(double min_p, double max_p,
                                           double p1, double p2,
                                           double cv1, double cv2);
    void testPerformance();
    Transform transf;

    Matrix Q, Qt, R;
    Vector gnew, w, z, dx;
    Vector g, gamma;
    vector<double> x;
};

#endif /*BASESPLINE_H_*/
