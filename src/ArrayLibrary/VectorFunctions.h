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
 
 $Id: VectorFunctions.h,v 1.2 2008/05/07 21:25:07 lukall Exp $
 
 *******************************************************************************/
#ifndef _VectorFunctions_H
#define _VectorFunctions_H

#include <iostream>
#include "ArrayLibrary.h"

using namespace std;

// sequence utilities
OrderedArray<int> seq(int lowest, int highest);
Vec seq(double lowest, double highest, double step);

// operators and functions for Vec, PackedSet, PackedVec
double dot(const Vec & lhs, const Vec & rhs);
double norm(const Vec & lhs);

double min(const Vec & lhs);
OrderedArray<int> argmin(const Vec & lhs);

double max(const Vec & lhs);
OrderedArray<int> argmax(const Vec & lhs);

const Vec & operator +=(Vec & lhs, const Vec & rhs);
const Vec & operator +=(Vec & lhs, double rhs);

const Vec & operator -=(Vec & lhs, const Vec & rhs);
const Vec & operator -=(Vec & lhs, double rhs);

const Vec & operator *=(Vec & lhs, const Vec & rhs);
const Vec & operator *=(Vec & lhs, double rhs);

const Vec & operator /=(Vec & lhs, const Vec & rhs);
const Vec & operator /=(Vec & lhs, double rhs);

Vec operator +(const Vec & lhs, const Vec & rhs);
Vec operator +(const Vec & lhs, double rhs);

Vec operator +(double lhs, const Vec & rhs);

Vec operator -(const Vec & lhs, const Vec & rhs);
Vec operator -(const Vec & lhs, double rhs);

Vec operator -(double lhs, const Vec & rhs);

Vec operator *(const Vec & lhs, const Vec & rhs);
Vec operator *(const Vec & lhs, double rhs);

Vec operator *(double lhs, const Vec & rhs);

Vec operator /(const Vec & lhs, const Vec & rhs);
Vec operator /(const Vec & lhs, double rhs);

Vec operator /(double lhs, const Vec & rhs);

// unary - operator
Vec operator -(const Vec & rhs);

// projection functions
Vec projectSingle(const Vec & basis, const Vec & v);

#endif

