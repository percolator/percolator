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
 
 $Id: PackedMatrix.h,v 1.4 2008/05/15 21:31:51 lukall Exp $
 
 *******************************************************************************/
#ifndef _PackedMatrix_H
#define _PackedMatrix_H

#include <assert.h>
#include "ArrayLibrary.h"

using namespace std;

typedef Array<PackedVec> PackedMatrix;

PackedMatrix identityPackedMatrix(int n);
PackedMatrix packedMatrixInverse(PackedMatrix inv);

PackedMatrix diagonalPacked(const Vec & v);
PackedMatrix transpose(const PackedMatrix & mat);
PackedMatrix matrixInverse(PackedMatrix mat);
//template<class T> void solveEquation(PackedMatrix& mat,Array<T>& res);

PackedMatrix matrixMult(const PackedMatrix & lhs, const PackedMatrix & rhs);
PackedMatrix operator *(const PackedMatrix & lhs, const PackedMatrix & rhs);
PackedMatrix operator *(const double & lhs, const PackedMatrix & rhs);
Vec          operator *(const PackedMatrix & lhs, const Vec & rhs);
PackedMatrix matrixAdd(const PackedMatrix & lhs, const PackedMatrix & rhs);
PackedMatrix operator +(const PackedMatrix & lhs, const PackedMatrix & rhs);
PackedMatrix makePackedMatrix(int a, int b);
PackedMatrix makePackedMatrix(int a);
size_t numCol(const PackedMatrix & mat);

ostream & operator <<(ostream & os, const PackedMatrix & rhs);

template<typename T> void solveEquation(PackedMatrix& mat,Array<T>& res) {
  assert(mat.size()==res.size());
  // Current implementation requires a quadratic mat
  int nCol = mat.size();

  PackedVec nonEmpty;
  int col,row,rowPos;
  for (col= 0; col<nCol; col++) {
    // find the non-null elements in this column (below row "col")
    nonEmpty.clear();
    int pivotPos(-1);
    for(row=col;row<mat.size();row++) {
      for(rowPos=mat[row].packedSize();rowPos--;) {
        if (mat[row].index(rowPos)==col) {
          nonEmpty.push_back(row,mat[row][rowPos]);
          if (row==col)
            pivotPos=nonEmpty.packedSize()-1;
        }
      }
    }
    //find most significant row
    double maxVal(0.0);
    int maxRow(-1),maxRowPos(-1);
    for (rowPos=nonEmpty.packedSize();rowPos--;) {
      double val = nonEmpty[rowPos];
      if (fabs(val)>fabs(maxVal)) {
        maxVal = val;
        maxRow = nonEmpty.index(rowPos);
        maxRowPos = rowPos;
      }    
    }
    // Put the most significant row at row "col"
    if (maxRow!=col) {
      swap(mat[col],mat[maxRow]);
      swap(res[col],res[maxRow]);
      if (pivotPos>=0)
        swap(nonEmpty[maxRowPos],nonEmpty[pivotPos]); 
    }
    // Divide the row with maxVal
    mat[col] /= maxVal;
    res[col] /= maxVal;
    // subtract the row from other rows
    for (rowPos=nonEmpty.packedSize();rowPos--;) {
      row = nonEmpty.index(rowPos);
      if (row == col) continue;
      // If the pivotRow was empty (prior to swap) at col=row do not process this row
      if (pivotPos<0 && row == maxRow) continue;
      double val = nonEmpty[rowPos];
      mat[row] -= val*mat[col];
      res[row] -= val*res[col];
    }
  }
  // Go bottom up and clear upper halfmatrix
  
  for (col=mat.size(); col--;) {
    nonEmpty.clear();
    for(row=0;row<col;row++) {
      for(rowPos=mat[row].packedSize();rowPos--;) {
        if (mat[row].index(rowPos)==col) {
          nonEmpty.push_back(row,mat[row][rowPos]);
        }
      }
    }
    // subtract the row from other rows
    for (rowPos=nonEmpty.packedSize();rowPos--;) {
      row = nonEmpty.index(rowPos);
      double val = nonEmpty[rowPos];
      mat[row] -= val*mat[col];
      res[row] -= val*res[col];
    }
  }
  return;
}

#endif

