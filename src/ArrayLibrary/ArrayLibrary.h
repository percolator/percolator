/*******************************************************************************
 Copyright (c) 2008-9 Oliver Serang

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

 $Id: ArrayLibrary.h,v 1.6 2009/03/31 09:26:57 lukall Exp $

 *******************************************************************************/
#ifndef _ArrayLibrary_H
#define _ArrayLibrary_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <iostream>
#include <math.h>
#include <vector>

#ifdef DEBUG
#define SAFE_ARRAYS
#endif

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

