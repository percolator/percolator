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
 
 $Id: Transform.h,v 1.2 2008/05/07 21:25:08 lukall Exp $
 
 *******************************************************************************/
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
