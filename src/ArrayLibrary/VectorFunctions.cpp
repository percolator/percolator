#include "VectorFunctions.h"

OrderedArray<int> seq(int lowest, int highest)
{
  OrderedArray<int> result( highest - lowest + 1 );
  
  int k;
  for (k=0; k<result.size(); k++)
    {
      result[k] = lowest + k;
    }

  return result;
}

Vec seq(double lowest, double highest, double step)
{
  Array<double> result( int( (highest - lowest) / step ) + 1 );

  int k;
  for (k=0; k<result.size(); k++)
    {
      result[k] = lowest + step * k;
    }

  return result;
}

double dot(const Vec & lhs, const Vec & rhs)
{
  sizeCheck(lhs, rhs);

  int k;
  
  double sum = 0;
  for (k=0; k<lhs.size(); k++)
    {
      sum += lhs[k] * rhs[k];
    }
  return sum;
}

double norm(const Vec & lhs)
{
  return sqrt( dot(lhs, lhs) );
}

double min(const Vec & lhs)
{
  double lowest = Numerical::infinity;
  
  modify(Arithmetic::minEq<double>, lowest, lhs);

  return lowest;
}

OrderedArray<int> argmin(const Vec & lhs)
{
  double lowest = min(lhs);
  return suchThat(Numerical::isZero, map(Arithmetic::sub<double>, lhs, lowest) );
}

double max(const Vec & lhs)
{
  double highest = -Numerical::infinity;
  
  modify(Arithmetic::maxEq<double>, highest, lhs);

  return highest;
}

OrderedArray<int> argmax(const Vec & lhs)
{
  double highest = max(lhs);
  return suchThat(Numerical::isZero, map(Arithmetic::sub<double>, lhs, highest) );
}

const Vec & operator +=( Vec & lhs, const Vec & rhs)
{
  modify(Arithmetic::addEq<double>, lhs, rhs);
  return lhs;
}

const Vec & operator +=( Vec & lhs, double rhs)
{
  modify(Arithmetic::addEq<double>, lhs, rhs);
  return lhs;
}

const Vec & operator -=( Vec & lhs, const Vec & rhs)
{
  modify(Arithmetic::subEq<double>, lhs, rhs);
  return lhs;
}

const Vec & operator -=( Vec & lhs, double rhs)
{
  modify(Arithmetic::subEq<double>, lhs, rhs);
  return lhs;
}

const Vec & operator *=( Vec & lhs, const Vec & rhs)
{
  modify(Arithmetic::multEq<double>, lhs, rhs);
  return lhs;
}

const Vec & operator *=( Vec & lhs, double rhs)
{
  modify(Arithmetic::multEq<double>, lhs, rhs);
  return lhs;
}

const Vec & operator /=( Vec & lhs, const Vec & rhs)
{
  modify(Arithmetic::divEq<double>, lhs, rhs);
  return lhs;
}

const Vec & operator /=( Vec & lhs, double rhs)
{
  modify(Arithmetic::divEq<double>, lhs, rhs);
  return lhs;
}

Vec operator +(const Vec & lhs, const Vec & rhs)
{
  return map(Arithmetic::add<double>, lhs, rhs);
}

Vec operator +(const Vec & lhs, double rhs)
{
  return map(Arithmetic::add<double>, lhs, rhs);
}

Vec operator +(double lhs, const Vec & rhs)
{
  return map(Arithmetic::add<double>, lhs, rhs);
}

Vec operator -(const Vec & lhs, const Vec & rhs)
{
  return map(Arithmetic::sub<double>, lhs, rhs);
}

Vec operator -(const Vec & lhs, double rhs)
{
  return map(Arithmetic::sub<double>, lhs, rhs);
}

Vec operator -(double lhs, const Vec & rhs)
{
  return map(Arithmetic::sub<double>, lhs, rhs);
}

Vec operator *(const Vec & lhs, const Vec & rhs)
{
  return map(Arithmetic::mult<double>, lhs, rhs);
}


Vec operator *(const Vec & lhs, double rhs)
{
  return map(Arithmetic::mult<double>, lhs, rhs);
}

Vec operator *(double lhs, const Vec & rhs)
{
  return map(Arithmetic::mult<double>, lhs, rhs);
}

Vec operator /(const Vec & lhs, const Vec & rhs)
{
  return map(Arithmetic::div<double>, lhs, rhs);
}

Vec operator /(const Vec & lhs, double rhs)
{
  return map(Arithmetic::div<double>, lhs, rhs);
}

Vec operator /(double lhs, const Vec & rhs)
{
  return map(Arithmetic::div<double>, lhs, rhs);
}

// unary
Vec operator -(const Vec & rhs)
{
  return map(Arithmetic::sub<double>, 0.0, rhs);
}

// projection functions
Vec projectSingle(const Vec & basis, const Vec & v)
{
  return ( dot(basis, v)/dot(basis, basis) ) * basis;
}

