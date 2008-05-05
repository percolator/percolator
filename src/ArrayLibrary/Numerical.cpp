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
