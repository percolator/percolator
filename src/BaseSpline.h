/*******************************************************************************
 Copyright (c) 2008 Lukas Käll

 Permission is hereby granted, free of charge, to any person
 obtaining a copy of this software and associated documentation
 files (the "Software"), to deal in the Software without
 restriction, including without limitation the rights to use,
 copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the
 Software is furnished to do so, subject to the following
 conditions:

 The above copyright notice and this permission notice shall be
 included in all copies or substantial portions of the Software. 
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 OTHER DEALINGS IN THE SOFTWARE.
 
 $Id: BaseSpline.h,v 1.5 2008/05/20 00:24:43 lukall Exp $
 
 *******************************************************************************/

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
  static double convergeEpsilon;
  static double stepEpsilon;
  void iterativeReweightedLeastSquares();
  void predict(const vector<double>& x, vector<double>& predict);
  void setData(const vector<double>& x);
  double predict(double xx) {return splineEval(xx);}
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
