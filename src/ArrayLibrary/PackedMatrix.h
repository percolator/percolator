#ifndef _PackedMatrix_H
#define _PackedMatrix_H

#include "PackedVec.h"

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
    for(row=col;row--;) {
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

