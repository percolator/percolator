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
 
 $Id: OrderedArray.cpp,v 1.3 2008/05/20 00:24:43 lukall Exp $
 
 *******************************************************************************/

#include "OrderedArray.h"

template <typename T>
OrderedArray<T>::OrderedArray(int n) :
  Array<T>(n)
{ 
}

template <typename T>
void OrderedArray<T>::push_back(const T & element)
{
#ifdef SAFE_ARRAYS
  if ( Array<T>::size() != 0 && !( Array<T>::back() < element ) ) {
    cerr << "Trying to insert " << element << " after " << Array<T>::back() << ", that is not in order." << endl;
    cerr.flush();
    throw OutOfOrderException();
  }
#endif

  Array<T>::push_back(element);
}

