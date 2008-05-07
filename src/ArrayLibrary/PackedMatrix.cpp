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
 
 $Id: PackedMatrix.cpp,v 1.2 2008/05/07 21:25:08 lukall Exp $
 
 *******************************************************************************/
#include<cstdlib>
using namespace std;
#include "PackedMatrix.h"

PackedMatrix identityPackedMatrix(int n)
{
  PackedMatrix result(n);

  int k;
  for (k=0; k<n; k++)
    {
      result[k] = PackedVec();
      result[k].push_back(k,1.0);
    }

  return result;
}

PackedMatrix makePackedMatrix(int a, int b)
{
  return makePackedMatrix(a);
}

PackedMatrix makePackedMatrix(int a)
{
  PackedMatrix res(a);
  int k;
  for (k=0; k<a; k++)
    {
      res[k] = PackedVec();
    }

  return res;
}

ostream & operator <<(ostream & os, const PackedMatrix & rhs)
{
  if (rhs.size() == 0)
  {
    os << "0";
    return os;
  }
  os << "{";
  for (int row=0; row<rhs.size(); row++)
  {
    os << rhs[row];
    if ( row != rhs.size() - 1 )
      os << "," << endl;
  }
  os << "}" << endl; 
  return os;
}



PackedMatrix matrixMult(const PackedMatrix & lhs, const PackedMatrix & rhs)
{
  PackedMatrix res = makePackedMatrix(lhs.size());
  PackedMatrix trhs = transpose(rhs);
  
  int row, col;
  for (row = 0; row < lhs.size(); row++)
  {
    for (col = 0; col < trhs.size(); col++)
	{
      double prod = pDot(lhs[row],trhs[col]);
      if (Numerical::isNonzero(prod)) {
        res[row].push_back(col,prod);
      }
	}
  }
  return res;
}

PackedMatrix matrixAdd(const PackedMatrix & lhs, const PackedMatrix & rhs)
{
  PackedMatrix res = makePackedMatrix(lhs.size());
  
  int row;
  for (row = 0; row < lhs.size(); row++)
  {
    res[row] = lhs[row] +rhs[row];
  }
  return res;
}

PackedMatrix operator *(const double & lhs, const PackedMatrix & rhs)
{
  PackedMatrix res = makePackedMatrix(rhs.size());
  
  int row;
  for (row = 0; row < rhs.size(); row++)
  {
    res[row] = lhs*rhs[row];
  }
  return res;
}


PackedMatrix operator *(const PackedMatrix & lhs, const PackedMatrix & rhs)
{
  return matrixMult(lhs, rhs);
}
Vec operator *(const PackedMatrix & lhs, const Vec & rhs)
{
//  PackedMatrix tlhs = transpose(lhs);
  return map<double,PackedVec,Vec>(pDot,lhs,rhs);
}

PackedMatrix operator +(const PackedMatrix & lhs, const PackedMatrix & rhs)
{
  return matrixAdd(lhs, rhs);
}

PackedMatrix matrixInverse(PackedMatrix mat)
{
  PackedMatrix inv = identityPackedMatrix(mat.size());
  solveEquation<PackedVec>(mat,inv);
  return inv;
}

size_t numCol(const PackedMatrix & mat)
{
  int col=-1;
  for (int i=0; i<mat.size();++i) {
    int last=mat[i].packedSize()-1;
    if (last>=0)
      col=max(col,mat[i].index(last));
  }

  return (col+1);
}


PackedMatrix transpose(const PackedMatrix & mat)
{
  PackedMatrix res =makePackedMatrix( numCol(mat) );
  
  int i,j;
  for (i=0; i<mat.size();++i) {
    for (j=0;j<mat[i].packedSize();++j) {
      res[mat[i].index(j)].push_back(i,mat[i][j]);
    }
  }

  return res;
}

PackedMatrix diagonalPacked(const Vec & v)
{
  PackedMatrix res( v.size() );
  int k;
  for (k=0; k<res.size(); k++)
    {
      res[k] = PackedVec();
      res[k].push_back(k,v[k]);
    }

  return res;
}


