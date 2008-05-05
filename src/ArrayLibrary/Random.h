#ifndef _Random_H
#define _Random_H

#include <iostream>
#include <math.h>

#include "Numerical.h"

using namespace std;

class Random
{
public:
  static double uniform(double a, double b);

  static double standardNormal();
  static double normal(double mean, double var);
  
  class OrderException {};
private:

};

#endif

