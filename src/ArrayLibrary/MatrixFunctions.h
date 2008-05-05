#ifndef _MatrixFunctions_H
#define _MatrixFunctions_H

#include <iostream>
#include "ArrayLibrary.h"

using namespace std;

typedef Array<Vec> Matrix;

class MatrixDimensionException {};

Matrix identityMatrix(int n);
Matrix matrixInverse(Matrix inv);

Matrix diagonal(const Vec & v);
Matrix transpose(const Matrix & mat);

Matrix matrixMult(const Matrix & lhs, const Matrix & rhs);
Matrix operator *(const Matrix & lhs, const Matrix & rhs);
Matrix makeMatrix(int a, int b);

#endif

