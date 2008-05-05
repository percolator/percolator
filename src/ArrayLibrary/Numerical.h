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

