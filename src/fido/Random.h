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
  
  inline static void setSeed(unsigned long s) { seed_ = s; }
  static unsigned long lcg_rand();
private:
  static Numerical samplingChecker;
  static unsigned long seed_;
};

#endif

