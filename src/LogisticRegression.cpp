/*******************************************************************************
 Copyright (c) 2008-9 Lukas Käll

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

 $Id: LogisticRegression.cpp,v 1.7 2009/02/12 14:31:47 lukall Exp $

 *******************************************************************************/
#include<iterator>
#include<vector>
#include<algorithm>
#include<numeric>
#include<functional>
#include<fstream>
#include<sstream>
using namespace std;

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
    assert(isfinite(g[ix]));
    assert(isfinite(e));
    assert(isfinite(p[ix]));
    assert(isfinite(w[ix]));
    assert(isfinite(z[ix]));
  }
}

void LogisticRegression::initg() {
  BaseSpline::initg();
  int n=x.size();
  p.resize(n);
  for (int ix=g.size();ix--;) {
    double p = (y[ix]+0.05)/(m[ix]+0.1);
    gnew[ix] = log(p/(1-p));
    assert(isfinite(p));
    assert(isfinite(g[ix]));
  }
//#define OUTPUT_DEBUG_FILES
#undef OUTPUT_DEBUG_FILES
#ifdef OUTPUT_DEBUG_FILES
  ofstream drFile("decoyRate.bins",ios::out),xvalFile("xvals.bins",ios::out);
  ostream_iterator<double> xvalIt(xvalFile,"\n");

  copy(x.begin(),x.end(),xvalIt);

  for(size_t yix=0;yix<y.size();++yix) {
    drFile << y[yix]/(double) m[yix] << endl;
  }
  drFile.close();
#endif
}



