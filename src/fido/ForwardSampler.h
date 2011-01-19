// Written by Oliver Serang 2009
// see license for more information

#ifndef _ForwardSampler_H
#define _ForwardSampler_H

#include "LPSolver.h"

class ForwardSampler : public LPSolver
{
public:
  ForwardSampler();
  ForwardSampler(char * fname);
  Vector solve();

  //  void loadArray(istream & fin);
  //  void load(istream & fin);

  void saveArray(ostream & fout);

  void checkFeasibility() const;
protected:
  void descentStep();
  
  void clear();

  // returns false if optimal
  bool randomStep();

  int advance(const Vector & p);
  void addConstraint(int index);
  void addConstraintLinearProjection(int index);
  void addConstraintGramSchmidt(int index);

  //  void normalizeLP();
  //  void addPositivityConstraints();

  Vector getProjection() const;

  Matrix getSpanningVectors() const;
  Matrix getOldSpanning();

  Matrix S;
  Matrix SSTransInverse;

  Numerical hitChecker;

  Set currentGliders;

  GramSchmidt gs;
};

#endif

