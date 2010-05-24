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

 $Id: OrderedArray.h,v 1.4 2009/03/30 03:13:31 cegrant Exp $

 *******************************************************************************/

// this class behaves like Array, except it
// requires a total order
// a < operator is required for T
//
// objects must be manually inserted in order

#ifndef _OrderedArray_H
#define _OrderedArray_H

template <typename T> class Array;

template <class T>
class OrderedArray : public Array<T> {
     public:
          // constructors
          OrderedArray() {}
          OrderedArray(int n);

          virtual void push_back(const T& element);

          class OutOfOrderException {};
     protected:

};

#endif

