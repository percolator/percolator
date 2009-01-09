/*******************************************************************************
 Copyright (c) 2008-9 Oliver Serang

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
 
 $Id: Random.cpp,v 1.4 2009/01/09 14:41:00 lukall Exp $
 
 *******************************************************************************/
#include <cstdlib>
#include "Random.h"

double Random::uniform(double a, double b)
{
  if ( a > b )
    {
      if ( Numerical::isPos(a-b) )
	throw OrderException();
      else
	{
	  // they are negligibly close
	  // just return one of the limits
	  return b;
	}
    }

  double u = rand() % 1000000 / 999999.0;
  return u * (b-a) + a;
}

double Random::standardNormal()
{
  double u1 = uniform(0.0, 1.0);
  double u2 = uniform(0.0, 1.0);

  double z = sqrt( -2.0 * log(u1) ) * cos(2*Pi*u2);

  return z;
}

double Random::normal(double mean, double var)
{
  return standardNormal() * sqrt(var) + mean;
}
