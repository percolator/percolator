/*******************************************************************************
 Copyright (c) 2008 Oliver Serang

 Permission is hereby granted, free of charge, to any person
 obtaining a copy of this software and associated documentation
 files (the "Software"), to deal in the Software without
 restriction, including without limitation the rights to use,
 copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the
 Software is furnished to do so, subject to the following
 conditions:

 The above copyright notice and this permission notice shall be
 included in all copies or substantial portions of the Software. 
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 OTHER DEALINGS IN THE SOFTWARE.
 
 $Id: MatrixFunctions.cpp,v 1.2 2008/05/07 21:25:08 lukall Exp $
 
 *******************************************************************************/
#include "MatrixFunctions.h"

Matrix identityMatrix(int n)
{
  Matrix result(n);

  int k;
  for (k=0; k<n; k++)
    {
      result[k] = Vec(n, 0.0);
      result[k][k] = 1.0;
    }

  return result;
}

Matrix makeMatrix(int a, int b)
{
  Matrix res(a);
  int k;
  for (k=0; k<a; k++)
    {
      res[k] = Vec(b, 0.0);
    }

  return res;
}

Matrix matrixMult(const Matrix & lhs, const Matrix & rhs)
{
  Matrix res = makeMatrix(lhs.size(), rhs[0].size());

  int row, col;
  for (row = 0; row < res.size(); row++)
    {
      for (col = 0; col < res[0].size(); col++)
	{
	  //	  res[row][col] = dot(lhs[row], getCol(rhs, col) );

	  double tot = 0.0;
	  for (int k=0; k<lhs[row].size(); k++)
	    {
	      tot += lhs[row][k] * rhs[k][col];
	    }
	  res[row][col] = tot;
	}
    }

  return res;
}

Matrix operator *(const Matrix & lhs, const Matrix & rhs)
{
  return matrixMult(lhs, rhs);
}

Matrix matrixInverse(Matrix inv)
{
  int row;
  Matrix id = identityMatrix(inv.size());
  for (row = 0; row<inv.size(); row++)
    {
      id[row] /= inv[row][row];
      inv[row] /= inv[row][row];
      
      int otherRow;
      for (otherRow = 0; otherRow < inv.size(); otherRow++)
	{
	  if ( row != otherRow )
	    {
	      id[otherRow] -= inv[otherRow][row]*id[row];
	      inv[otherRow] -= inv[otherRow][row]*inv[row];
	    }
	}
    }

  return id;
}

Matrix transpose(const Matrix & mat)
{
  Matrix res( mat[0].size() );
  
  int k;
  for (k=0; k<res.size(); k++)
    {
      res[k] = getCol(mat, k);
    }

  return res;
}

Matrix diagonal(const Vec & v)
{
  Matrix res( v.size() );
  int k;
  for (k=0; k<res.size(); k++)
    {
      res[k] = Vec(v.size(), 0.0);
      res[k][k] = v[k];
    }

  return res;
}
