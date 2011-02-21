// Written by Oliver Serang 2009
// see license for more information

#include "Matrix.h"

const Matrix & Matrix::operator =(const Array<Array<double> > & rhs)
{
  // note: this will not accept arrays with no rows
  Matrix temp(rhs.size(), rhs[0].size());
  for (int k=0; k<rhs.size(); k++)
    {
      temp.rows[k] = Vector(rhs[k]);
    }

  *this = temp;
  return *this;
}

void Matrix::add(const Vector & rhs)
{
  if ( numRows() > 0 && rhs.size() != numCols() )
    {
      cerr << "Cannot add vector with size " << rhs.size() <<
	" to a matrix with " << numCols() << " columns" << endl;
      throw MatrixIntegrityException();
    }

  if ( numRows() == 0 )
    {
      cols = rhs.size();
    }

  rows.add(rhs);
}

Vector operator *(const Matrix & lhs, const Vector & rhs)
{
  //  if ( lhs.numCols() == 0 )
  //    throw Matrix::SizeException();

  Array<double> result( lhs.numRows() );

  int k;
  for (k=0; k<lhs.numRows(); k++)
    {
      result[k] = lhs[k] * rhs;
    }

  return Vector(result);
}

Vector operator *(const Vector & lhs, const Matrix & rhs)
{
  // performance note: may want to switch to transposed indices later
  
  Matrix trans = rhs.transpose();
  return trans * lhs;
}

void Matrix::displayMatrix() const
{
  //  cout << endl;
  int k;
  cout << "corner\t";
  for ( k=0; k<numCols(); k++)
    {
      cout << "col" << k;
      if ( k != numCols() - 1 )
	cout << "\t";
    }
  cout << endl;
  for ( k=0; k<numRows(); k++)
    {
      cout << "row" << k << "\t";
      rows[k].displayVector();
    }
}

void Matrix::verifyMatrixIntegrity() const
{
  for (int k=0; k<numRows(); k++)
    {
      if ( rows[k].size() != numCols() )
	throw MatrixIntegrityException();
    }
}

ostream & operator <<(ostream & os, const Matrix & rhs)
{
  os << "( " << rhs.numCols() << " : " << rhs.rows << " )";
  return os;
}

istream & operator >>(istream & is, Matrix & rhs)
{
  char c;

  is >> c;
  if ( c != '(' )
    throw Matrix::MatrixFormatException();

  is >> rhs.cols;

  is >> c;
  if ( c != ':' )
    throw Matrix::MatrixFormatException();

  is >> rhs.rows;

  is >> c;
  if ( c != ')' )
    throw Matrix::MatrixFormatException();

  rhs.verifyMatrixIntegrity();

  return is;
}

Matrix Matrix::operator -() const
{
  Matrix result = *this;

  for (int k=0; k<result.numRows(); k++)
    {
      result.rows[k] = -result[k];
    }
  
  return result;
}


Matrix Matrix::transpose() const
{
  Array<Array<double> > resultArray(numCols(), Array<double>(numRows(), 0.0) );

  int k,j;
  for (k=0; k<numRows(); k++)
    {
      for (j=0; j<numCols(); j++)
	{
	  resultArray[j][k] = rows[k][j];
	}
    }

  return Matrix(resultArray);
}

Matrix Matrix::operator *(double scale) const
{
  Matrix result = *this;
  for (int k=0; k<result.numRows(); k++)
    {
      result[k] *= scale;
    }

  return result;
}

Matrix operator *(const Matrix & lhs, const Matrix & rhs)
{
  if ( lhs.numCols() != rhs.numRows() )
    {
      throw Matrix::SizeException();
    }

  Matrix result( lhs.numRows(), rhs.numCols() );

  // performance note: you might want to rewrite this without using
  // transpose
  Matrix rhsTrans = rhs.transpose();

  int k,j;
  for (k=0; k<lhs.numRows(); k++)
    {
      Array<double> nextRow( rhs.numCols() );
      for (j=0; j<rhsTrans.numRows(); j++)
	{
	  //	  nextRow[j] = lhs[k].unpack() * rhsTrans[j].unpack();
	  nextRow[j] = lhs[k] * rhsTrans[j];
	}
      result[k] = Vector(nextRow);
    }

  return result;
}

const Matrix & Matrix::operator +=(const Matrix & rhs)
{
  if ( numRows() != rhs.numRows() && numCols() != rhs.numCols() )
    {
      throw SizeException();
    }
  
  int k;
  for (k=0; k<numRows(); k++)
    {
      rows[k] += rhs[k];
    }

  return *this;
}

const Matrix & Matrix::operator -=(const Matrix & rhs)
{
  if ( numRows() != rhs.numRows() && numCols() != rhs.numCols() )
    {
      throw SizeException();
    }
  
  int k;
  for (k=0; k<numRows(); k++)
    {
      rows[k] -= rhs[k];
    }

  return *this;
}

Matrix ShermanMorrison(const Matrix & mat, const Vector & u, const Vector & v)
{
  Matrix result = mat;

  double coef =  1/(1+v*mat*u);

  if ( Vector::sparseChecker.isZero(coef) )
    {
      cerr << "Sherman morrison has a division by a very small amount; still LI?" << endl;
      exit(1);
    }

  // performance note: the outer product should be more efficient,
  // especially for sparse vectors
  Matrix lhs = Matrix(coef * (mat*u) ).transpose();
  Matrix rhs = Matrix(v*mat);

  //  cout << "\tTerm to add in ShermanMorrison is: " << endl << lhs*rhs << endl << endl;
  result -= lhs * rhs;

  return result;
}

double Matrix::frobeniusNorm() const
{
  double tot = 0;

  for (int k=0; k<numRows(); k++)
    {
      tot += pow(rows[k].norm(), 2.0);
    }

  return sqrt(tot);;
}

Matrix operator +(const Matrix & lhs, const Matrix & rhs)
{
  Matrix result = lhs;
  result += rhs;
  return result;
}

Matrix operator -(const Matrix & lhs, const Matrix & rhs)
{
  Matrix result = lhs;
  result -= rhs;
  return result;
}

Set Matrix::prec(double val) const
{
  Set result;

  for (int k=0; k<numRows(); k++)
    {
      if ( rows[k].prec(val) )
	result.add( k );
    }

  return result;
}

Set Matrix::succ(double val) const
{
  Set result;

  for (int k=0; k<numRows(); k++)
    {
      if ( rows[k].succ(val) )
	result.add( k );
    }

  return result;
}

Set Matrix::precEq(double val) const
{
  Set result;

  for (int k=0; k<numRows(); k++)
    {
      if ( rows[k].precEq(val) )
	result.add( k );
    }

  return result;
}

Set Matrix::succEq(double val) const
{
  Set result;

  for (int k=0; k<numRows(); k++)
    {
      if ( rows[k].succEq(val) )
	result.add( k );
    }

  return result;
}

void Matrix::embedInIdentity()
{
  if ( numRows() != numCols() )
    throw NonSquareMatrixException();

  for (int k=0; k<numRows(); k++)
    {
      rows[k].resize( rows[k].size() + 1 );
    }
  
  Vector idRow( numCols() + 1 );
  idRow.addElement( numCols(), 1.0 );

  cols++;

  add( idRow );
}

Matrix Matrix::diagonalMatrix( const Vector & v )
{
  Matrix result( v.size() , v.size() );
  
  for (int k=0; k<v.size(); k++)
    {
      result.rows[k].addElement( k, v[k] );
    }

  return result;
}

void Matrix::invert()
{
  if ( numRows() != numCols() )
    throw NonSquareMatrixException();

  Matrix id = identityMatrix( numRows() );

  //  cout << "Pre inverse: " << endl;
  //  displayMatrix();
  //  cout << endl;

  // eliminate
  for (int k=0; k<numRows(); k++)
    {
      eliminate(k, id);

      /***
      cout << "After row " << k << " mat and ID are: " << endl;
      displayMatrix();
      cout << endl;
      id.displayMatrix();
      cout << endl << endl;
      ***/
    }

  /***
  cout << "Hopefully looks like ID Matrix: " << endl;
  cout << *this << endl;
  displayMatrix();
  ***/

  *this = id;
}

void Matrix::eliminate(int k, Matrix & id)
{
  // first find the row with the largest element k
  for (int i=k+1; i<numRows(); i++)
    {
      if ( fabs( rows[i][k] ) > fabs( rows[k][k] ) )
	{
	  swap(rows[i], rows[k]);
	  swap(id.rows[i], id.rows[k]);
	}
    }

  // first normalize rows[k] so that it's first element is 1
  if ( Vector::sparseChecker.isZero(rows[k][k] ) )
    {
      cerr << "Error: largest value in column is " << rows[k][k] << endl;
      throw SingularMatrixException();
    }

  double valueK = rows[k][k];

  id.rows[k] /= valueK;
  rows[k] /= valueK;

  //  id.rows[k] *= 1/rows[k].leadingElement();
  //  rows[k] *= 1/rows[k].leadingElement();

  for (int j=0; j<numRows() && j < id.numRows(); j++)
    {
      if ( j == k )
	continue;

      //      int i = rows[k].leadingIndex();

      double mult = rows[j][k];

      //      if ( mult != 0.0 )
      if ( Vector::sparseChecker.isNonzero(mult) )
	{
	  id.rows[j].addEqScaled( -mult, id.rows[k]);
	  rows[j].addEqScaled( -mult, rows[k] );
	}
    }
}

void Matrix::eliminateRRE(int k, Matrix & id)
{
  // first find the row with the largest element k
  for (int i=k+1; i<numRows(); i++)
    {
      if ( fabs( rows[i][k] ) > fabs( rows[k][k] ) )
	{
	  swap(rows[i], rows[k]);
	  swap(id.rows[i], id.rows[k]);
	}
    }

  // first normalize rows[k] so that it's first element is 1
  if ( Vector::sparseChecker.isZero(rows[k][k] ) )
    {
      cerr << "Error: largest value in column is " << rows[k][k] << endl;
      throw SingularMatrixException();
    }

  double valueK = rows[k][k];

  id.rows[k] /= valueK;
  rows[k] /= valueK;

  //  id.rows[k] *= 1/rows[k].leadingElement();
  //  rows[k] *= 1/rows[k].leadingElement();

  for (int j=0; j<numRows() && j < id.numRows(); j++)
    {
      if ( j == k )
	continue;

      //      int i = rows[k].leadingIndex();

      double mult = rows[j][k];

      // note: this prevents unnecessary operations, but could be done
      // more cleanly later

      //      if ( mult != 0.0 )
      if ( Vector::sparseChecker.isNonzero(mult) )
	{
	  id.rows[j].addEqScaled( -mult, id.rows[k] );
	  rows[j].addEqScaled( -mult, rows[k] );
	}
    }
}

Array<Array<double> > Matrix::unpack() const
{
  Array<Array<double> > result;

  for (int k=0; k<numRows(); k++)
    {
      result.add( rows[k].unpack() );
    }

  return result;
}

Matrix Matrix::inverse() const
{
  Matrix temp = *this;
  temp.invert();
  return temp;
}

Array<Set> Matrix::transposeIndices() const
{
  Array<Set> transposedIndices(numCols());

  for (int k=0; k<numRows(); k++)
    {
      for (Set::Iterator iter = rows[k].beginNonzero(); iter != rows[k].endNonzero(); iter++)
	{
	  transposedIndices[ *iter ].add(k);
	}
    }

  return transposedIndices;
}

Matrix Matrix::identityMatrix(int n)
{
  Matrix result(n);
  for (int k=0; k<n; k++)
    {
      /***
      Array<double> idRow(n, 0.0);
      idRow[k] = 1.0;

      result.add( Vector(idRow) );
      ***/

      Vector idRowK(n);
      idRowK.addElement(k, 1.0);

      result.add(idRowK);
    }

  return result;
}

Matrix Matrix::RRE()
{
  Matrix subID = identityMatrix(numRows());
  for (int k=0; k<numCols(); k++)
    {
      eliminateRRE(k, subID);

    }

  return subID[ Set::FullSet(0, numCols()-1) ];
}

 /***
Matrix Matrix::transpose() const
{
  Matrix result(numCols(), numRows());

  for (int k=0; k<numRows(); k++)
    {
      for (int j=0; j<numCols(); j++)
	{
	  result.rows[j][k] = (*this)[k][j];
	}
    }

  return result;
}
 ***/

Vector Matrix::column(int k) const
{
  Vector result( numRows() );
  
  for (int i=0; i<numRows(); i++)
    {
      result.addElement(i, rows[i][k] );
    }

  return result;
}

Array<Array<double> > operator -(const Array<Array<double> > & rhs)
{
  Array<Array<double> > result = rhs;

  for (int k=0; k<result.size(); k++)
    {
      result[k] = -result[k];
    }
  
  return result;
}

double Matrix::maxMagnitude() const
{
  // valid if this is not an empty matrix ( it is fine if there are no
  // nonzeros )
  double maxVal = 0.0;

  for (int k=0; k<numRows(); k++)
    {
      const Vector & row = rows[k];

      for (Set::Iterator iter = row.beginNonzero(); iter != row.endNonzero(); iter++)
	{
	  maxVal = ::max( maxVal, fabs(row[ *iter ]) );
	}
    }

  return maxVal;
}

Vector solve(const Matrix & A, const Vector & b)
{
  // check if square

  return A.inverse()*b;
}

Vector fastSolve(const Matrix & A, const Vector & b)
{
  // check if square

  Vector cumulative = b;
  Vector lastLTerm = b;

  double beta = A.maxMagnitude()+1;
  double gamma = beta*2*A.numRows();

  Matrix L = Matrix::identityMatrix( A.numRows() ) - A*(1/gamma);

  // expansion: I+L+L^2+L^3+...
  // *b = b + Lb +LLb+...

  // note: untidy: fix later
  //  for (int k=1; k<Matrix::SeriesTerms; k++)
  for (int k=1; k<100; k++)
    {
      lastLTerm = L * lastLTerm;
      //      cumulative.addEqScaled( 1, b );
      //      cumulative.addEqScaled( -1, lastLTerm);
      cumulative.addEqScaled( 1, lastLTerm);
    }

  return cumulative / gamma;
}
