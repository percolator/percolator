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
