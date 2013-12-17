// Written by Oliver Serang 2009
// see license for more information

#ifndef _Random_H
#define _Random_H

#include <iostream>
#include <cmath>

#include "Array.h"
#include "Numerical.h"

using namespace std;

class Random
{
 public:
  static double uniform(double a, double b);
  static int inRange(int a, int b);

  static double standardNormal();
  static double normal(double mean, double var);

  static void fillRandomUniform(Array<double> & lhs, double low, double high);

  class SamplingException {};
private:
  static Numerical samplingChecker;
};

#endif

