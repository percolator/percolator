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
 
 $Id: Arithmetic.h,v 1.2 2008/05/07 21:25:07 lukall Exp $
 
 *******************************************************************************/

#ifndef _Arithmetic_H
#define _Arithmetic_H

#include <iostream>

using namespace std;


class Arithmetic
{
public:
  template <typename T>
    static T add(T lhs, T rhs);
  template <typename T>
    static void addEq(T & lhs, T rhs);

  template <typename T>
    static T sub(T lhs, T rhs);
  template <typename T>
    static void subEq(T & lhs, T rhs);

  template <typename T>
    static T mult(T lhs, T rhs);
  template <typename T>
    static void multEq(T & lhs, T rhs);

  template <typename T>
    static T div(T lhs, T rhs);
  template <typename T>
    static void divEq(T & lhs, T rhs);

  template <typename T>
    static T min(T lhs, T rhs);
  template <typename T>
    static void minEq(T & lhs, T rhs);

  template <typename T>
    static T max(T lhs, T rhs);
  template <typename T>
    static void maxEq(T & lhs, T rhs);

private:

};

#endif

