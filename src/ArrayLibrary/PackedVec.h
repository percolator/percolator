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
 
 $Id: PackedVec.h,v 1.2 2008/05/07 21:25:07 lukall Exp $
 
 *******************************************************************************/

#include "ArrayLibrary.h"
#include "Numerical.h"

#ifndef _PackedVec_H
#define _PackedVec_H

#include <iostream>

using namespace std;


class PackedVec
{
public:
  void push_back(int index, double d)
  {
    nonNull.push_back(index);
    data.push_back(d);
  }

  const double & operator[] (int i) const
  {
    return data[i];
  }

  double & operator[] (int i)
  {
    return data[i];
  }
  
  int index(int i) const
  {
    return nonNull[i];
  }
  
  int packedSize() const
  {
    return data.size();
  }

  void clear()
  {
    data.clear();
    nonNull.clear();
  }

  // exception classes
  class PackedVecFormatException {};
private:
  Array<double> data;
  PackedSet nonNull;
};

double norm(const PackedVec & rhs);

double pDot(const PackedVec & lhs, const Vec & rhs);
double pDot(const Vec & lhs, const PackedVec & rhs);
double pDot(const PackedVec & lhs, const PackedVec & rhs);

const Vec & operator += (Vec & lhs, const PackedVec & rhs);
const Vec & operator -= (Vec & lhs, const PackedVec & rhs);
const Vec & operator *= (Vec & lhs, const PackedVec & rhs);
const Vec & operator /= (Vec & lhs, const PackedVec & rhs);
const PackedVec & operator /= (PackedVec & lhs, const double & rhs);

PackedVec operator + (const PackedVec & lhs, const PackedVec & rhs);
PackedVec operator - (const PackedVec & lhs, const PackedVec & rhs);

const PackedVec & operator += (PackedVec & lhs, const PackedVec & rhs);
const PackedVec & operator -= (PackedVec & lhs, const PackedVec & rhs);

PackedVec operator *(double d, const PackedVec & rhs);

istream & operator >>(istream & is, PackedVec & rhs);
ostream & operator <<(ostream & os, const PackedVec & rhs);

PackedVec projectSingle(const PackedVec & basis, const PackedVec & v);
PackedVec projectSingle(const PackedVec & basis, const Vec & v);

#endif

