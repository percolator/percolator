// Written by Oliver Serang 2009
// see license for more information

#ifndef _Matrix_H
#define _Matrix_H

#include "Vector.h"

class Matrix
{
public:
  Matrix()
    {
      cols = 0;
    }
  Matrix(int r, int c)
    {
      rows = Array<Vector>(r, Vector(c) );
      cols = c;
    }
  Matrix(int c)
    {
      cols = c;
    }
  Matrix(const Array<Array<double> > & rhs)
    {
      *this = rhs;
    }
  Matrix(const Array<Vector> & rhs)
    {
      cols = 0;

      if ( rhs.size() > 0 )
	cols = rhs[0].size();

      rows = rhs;
    }
  explicit Matrix(const Vector & rhs)
    {
      cols = rhs.size();
      rows = Array<Vector>(1, rhs);
    }

  const Matrix & operator =(const Array<Array<double> > & rhs);

  virtual ~Matrix() {}

  virtual void add(const Vector & rhs);

  Matrix transpose() const;

  void reduce();

  Vector & operator [](int k) 
  {
    return rows[k];
  }

  const Vector & operator [](int k) const
  {
    return rows[k];
  }

  Matrix operator [](const Set & rhs) const
  {
    return Matrix(rows[rhs]);
  }

  int numRows() const
  {
    return rows.size();
  }
  int numCols() const
  {
    return cols;
  }

  Array<Set> transposeIndices() const;

  Set prec(double val) const;
  Set succ(double val) const;
  Set precEq(double val) const;
  Set succEq(double val) const;

  void embedInIdentity();

  void displayMatrix() const;

  void verifyMatrixIntegrity() const;

  void invert();
  Matrix inverse() const;

  void eliminate(int k, Matrix & id);
  void eliminateRRE(int k, Matrix & id);

  Matrix RRE();

  const Matrix & operator +=(const Matrix & rhs);
  const Matrix & operator -=(const Matrix & rhs);
  Matrix operator -() const;

  Matrix operator *(double scale) const;

  double maxMagnitude() const;

  double frobeniusNorm() const;

  static Matrix identityMatrix(int n);
  static Matrix diagonalMatrix(const Vector & vec);

  Vector column(int k) const;

  Array<Array<double> > unpack() const;

  friend ostream & operator <<(ostream & os, const Matrix & rhs);
  friend istream & operator >>(istream & is, Matrix & rhs);

  class NonSquareMatrixException {};
  class SizeException {};
  class MatrixFormatException {};
  class MatrixIntegrityException {};
  class SingularMatrixException {};

  friend class GramSchmidt;

protected:
  int cols;
  Array<Vector> rows;
};

Matrix ShermanMorrison(const Matrix & mat, const Vector & u, const Vector & v);
Vector fastSolve(const Matrix & A, const Vector & b);
Vector solve(const Matrix & A, const Vector & b);

Matrix operator *(const Matrix & lhs, const Matrix & rhs);
Vector operator *(const Matrix & lhs, const Vector & rhs);
Vector operator *(const Vector & lhs, const Matrix & rhs);

Matrix operator +(const Matrix & lhs, const Matrix & rhs);
Matrix operator -(const Matrix & lhs, const Matrix & rhs);

Array<Array<double> > operator -(const Array<Array<double> > & rhs);


#endif

