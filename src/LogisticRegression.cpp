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
 
 $Id: LogisticRegression.cpp,v 1.3 2008/05/20 00:24:43 lukall Exp $
 
 *******************************************************************************/

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



