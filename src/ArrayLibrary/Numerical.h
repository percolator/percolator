/*******************************************************************************
 Copyright (c) 2008 Oliver Serang

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
 
 $Id: Numerical.h,v 1.2 2008/05/07 21:25:08 lukall Exp $
 
 *******************************************************************************/

#ifndef _Numerical_H
#define _Numerical_H

#include <iostream>
#include <math.h>
#include <limits>

using namespace std;

const double Pi = 3.14159;

class Numerical
{
 public:
  static bool isPos(double d);
  static bool isNonpos(double d);
  static bool isNeg(double d);
  static bool isNonneg(double d);
  static bool isZero(double d);
  static bool isNonzero(double d);
  static bool isEqual(double a, double b);
  static bool isInequal(double a, double b);

  static double epsilon;
  static double infinity;
};

#endif

