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

