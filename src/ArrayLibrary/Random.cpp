#include "Random.h"

double Random::uniform(double a, double b)
{
  if ( a > b )
    {
      if ( Numerical::isPos(a-b) )
	throw OrderException();
      else
	{
	  // they are negligibly close
	  // just return one of the limits
	  return b;
	}
    }

  double u = rand() % 1000000 / 999999.0;
  return u * (b-a) + a;
}

double Random::standardNormal()
{
  double u1 = uniform(0.0, 1.0);
  double u2 = uniform(0.0, 1.0);

  double z = sqrt( -2.0 * log(u1) ) * cos(2*Pi*u2);

  return z;
}

double Random::normal(double mean, double var)
{
  return standardNormal() * sqrt(var) + mean;
}
