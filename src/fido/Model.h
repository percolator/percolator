// Written by Oliver Serang 2009
// see license for more information

#ifndef _Model_H
#define _Model_H

using namespace std;

#include <math.h>
#include "Combinatorics.h"

#include <iostream>

class Model
{
 protected:
 public:
  double alpha, beta, gamma;
  Model()
    {
      alpha = beta = gamma = -1;
    }
  Model(double a, double b, double g)
    {
      alpha = a;
      beta = b;
      gamma = g;
    }

  friend bool operator ==(const Model & lhs, const Model & rhs)
  {
    return lhs.alpha == rhs.alpha && lhs.beta == rhs.beta && lhs.gamma == rhs.gamma;
  }

  void setAlphaBeta(double a, double b)
  {
    alpha = a;
    beta = b;
  }

  double associatedEmission() const
  {
    return alpha;
  }

  double spontaneousEmission() const
  {
    return beta;
  }

  double probabilityNoEmissionFrom(int numActiveProts) const
  {
    //    return (1-spontaneousEmission())*pow( (1-associatedEmission()), numActiveProts);

    // using log for greater precision
    return pow(2.0, log2( 1-spontaneousEmission() )+numActiveProts * log2(1-associatedEmission()) );
  }

  double probabilityProteins( int totalProts, int activeProts ) const
  {
    //    return Combinatorics::binomial(totalProts, activeProts) * pow(gamma, activeProts) * pow(1-gamma, totalProts-activeProts);

    // using log for greater precision
    return pow(2.0, Combinatorics::logBinomial(totalProts, activeProts) + activeProts*log2(gamma) + (totalProts-activeProts) * log2(1-gamma) );
  }

  friend ostream & operator <<(ostream & os, const Model & m)
  {
    os << "alpha = " << m.alpha << ", \t beta = " << m.beta << ", \t gamma = " << m.gamma << endl;
    return os;
  }
};

class RealRange
{
 protected:
  double value, min, max, resolution;
  int count;
 public:
  RealRange()
    {
      min = 1;
      max = 0;
      resolution = -1;
    }
  RealRange(double small, double res, double large)
    {
      min = small;
      resolution = res;
      max = large;

      start();
    }

  // note: this does not include state!
  friend bool operator ==(const RealRange & lhs, const RealRange & rhs)
  {
    return lhs.min == rhs.min && lhs.max == rhs.max && lhs.resolution == rhs.resolution;
  }

  double getValue() const
  {
    return value;
  }

  void setValue(double v)
  {
    value = v;
  }
  
  void setCount(int c)
  {
    count = c;
  }

  void start()
  {
    value = min;
    count = 0;
  }

  void advance()
  {
    value += resolution;
    count++;
  }

  bool inRange() const
  {
    return value < (max - 1e-5);
  }

  int getCount() const
  {
    return count;
  }
  
  int maxCount() const
  {
    return int((max-min)/resolution + 1);
  }
};

class GridModel : public Model
{
 protected:
  RealRange alphaRange, betaRange;
  int count;
 public:
  GridModel()
    {}

 GridModel(RealRange aR, RealRange bR, double g):
  alphaRange(aR), betaRange(bR)
    {
      gamma = g;
    }

  // note: this does not include state!
  // to compare state only, you need to cast as a model
  friend bool operator ==(const GridModel & lhs, const GridModel & rhs)
  {
    return lhs.alphaRange == rhs.alphaRange && lhs.betaRange == rhs.betaRange && lhs.gamma == rhs.gamma;
  }

  void start()
  {
    alphaRange.start();
    betaRange.start();

    setAlphaBeta( alphaRange.getValue(), betaRange.getValue() );

    count = 0;
  }

  int maxRowIndex() const
  {
    return betaRange.maxCount();
  }

  int maxColIndex() const
  {
    return alphaRange.maxCount();
  }

  int getRowIndex() const
  {
    return betaRange.getCount();
  }

  int getColIndex() const
  {
    return alphaRange.getCount();
  }

  void advance()
  {
   // rectangular prior 
    alphaRange.advance();

    if ( ! alphaRange.inRange() )
      {
	betaRange.advance();
	alphaRange.start();
      }

    setAlphaBeta( alphaRange.getValue(), betaRange.getValue() );

    count++;
  }

  /***
  void advance()
  {
    // triangular prior (alpha >= beta)
    alphaRange.advance();

    if ( ! alphaRange.inRange() )
      {
	betaRange.advance();
	alphaRange.start();
	
	alphaRange.setValue( betaRange.getValue() );

	// note: assumes same bounding region
	alphaRange.setCount( betaRange.getCount() );
      }

    setAlphaBeta( alphaRange.getValue(), betaRange.getValue() );

    count++;
  }
  ***/

  int getCount() const
  {
    return count;
  }

  bool inRange() const
  {
    return alphaRange.inRange() && betaRange.inRange();
  }
};

#endif

