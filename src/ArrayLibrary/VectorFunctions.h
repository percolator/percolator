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

