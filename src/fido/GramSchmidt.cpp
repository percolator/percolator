// Written by Oliver Serang 2009
// see license for more information

#include "GramSchmidt.h"

Numerical GramSchmidt::zeroVectorChecker(1e-8);

void GramSchmidt::add(const Vector & rhs)
{
  Vector ortho = projectOrthogonal(rhs);

  if ( zeroVectorChecker.isZero( ortho.norm() ) )
    {
      cerr << "GramSchmidt projectiong norm is: " << ortho.norm() << endl;
      cerr << "\ttable had " << numRows() << " vectors and this was added: " << endl;
      cerr << rhs << endl << endl;

      throw ZeroVectorException();
    }

  Matrix::add( (1/ortho.norm()) * ortho);
  original.add(rhs);
}

Vector GramSchmidt::projectOrthogonal(const Vector & rhs) const
{
  Vector ortho = rhs;
  //  cout << "Proj in: " << rhs << endl;

  int k;
  for (k=0; k<numRows(); k++)
    {
      //            ortho = ortho - ( ((*this)[k] * ortho)/( (*this)[k] * (*this)[k] ) ) * (*this)[k];

      // they are unit vectors-- don't worry
      ortho.addEqScaled( -((*this)[k] * ortho) , (*this)[k] );
    }

  //  cout << "Proj out: " << ortho << endl;
  return ortho;
}

void GramSchmidt::swapToLast(int a)
{
  for (int k=a; k<numRows() - 1; k++)
    {
      swapWithNext(k);
    }
}

void GramSchmidt::swapWithNext(int k)
{
  // something is wrong with this code...

  Vector first = rows[k];
  Vector second = rows[k+1];

  second += projectSingle(rows[k], original[k+1]);
  first -= projectSingle(second, original[k]);

  rows[k] = second;
  rows[k+1] = first;
  swap( original.rows[k+1], original.rows[k] );
}

Vector GramSchmidt::projectSingle(const Vector & basis, const Vector & v) const
{
  return ( basis * v / (basis * basis) ) * basis;
}
