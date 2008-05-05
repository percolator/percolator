#ifndef _ArrayLibrary_H
#define _ArrayLibrary_H

#include <iostream>
#include <math.h>
#include <vector>

//#define SAFE_ARRAYS

using namespace std;

// forward declarations
template <typename T> class Array;
template <typename T> class OrderedArray;
template <typename T> class PackedArray;

#include "Array.h"
#include "OrderedArray.h"
#include "ArrayTypes.h"
#include "Numerical.h"
#include "Arithmetic.h"
#include "VectorFunctions.h"
#include "MatrixFunctions.h"
#include "Random.h"
#include "PackedVec.h"
#include "PackedSetFunctions.h"
#include "PackedMatrix.h"

#include "Array.cpp"
#include "OrderedArray.cpp"
#include "Arithmetic.cpp"

#endif

