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
 
 $Id: Numerical.cpp,v 1.2 2008/05/07 21:25:08 lukall Exp $
 
 *******************************************************************************/
#include "Numerical.h"

double Numerical::epsilon = 1e-5;
double Numerical::infinity = numeric_limits<double>::infinity();

bool Numerical::isPos(double d)
{
  return d > epsilon;
}

bool Numerical::isNonpos(double d)
{
  return d <= epsilon;
}

bool Numerical::isNeg(double d)
{
  return d < -epsilon;
}

bool Numerical::isNonneg(double d)
{
  return d >= -epsilon;
}

bool Numerical::isZero(double d)
{
  return fabs(d) <= epsilon;
}

bool Numerical::isNonzero(double d)
{
  return fabs(d) > epsilon;
}

bool Numerical::isEqual(double a, double b)
{
  return isZero( b - a );
}

bool Numerical::isInequal(double a, double b)
{
  return isNonzero( b - a );
}
