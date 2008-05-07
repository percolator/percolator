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
 
 $Id: Super.h,v 1.2 2008/05/07 21:25:08 lukall Exp $
 
 *******************************************************************************/
#include <vector>
#include <iostream>

using namespace std;

// forward declarations
template <typename T> class Array;
template <typename T> class OrderedArray;

/***
// specific types of Arrays
typedef OrderedArray<int> PackedSet;
typedef Array<double> Vec;
typedef PackedArray<double> PackedVec;
***/

// operators for Vec, PackedSet, PackedVec

template <typename T>
class Array
{
public:
  // constructors
  Array() {}

  Array(int n);
  Array(int n, const T & defaultValue);

  // important vector functions
  virtual const T & operator [] (int k) const;
  virtual T & operator [] (int k);
  virtual void push_back(const T & element);

  virtual T & back();
  virtual const T & back() const;

  // accessors
  int size() const;

  // novel functions

  // exception classes
  class OutOfBoundsException {};
  class SizeException {};
  
protected:
  void boundsCheck(int i) const;
  vector<T> data;
};



template <typename T>
class OrderedArray : public Array<T>
{
public:
  virtual void push_back(const T & element);

  class OutOfOrderException {};
protected:
  
};


Array<int> seq(int lowest, int highest);
Array<double> seq(double lowest, double highest, double step);

Array<int> seq(int lowest, int highest)
{
  Array<int> result( highest - lowest + 1, 0 );
  
  int k;
  for (k=0; k<result.size(); k++)
    {
      result[k] = lowest + k;
    }

  return result;
}

Array<double> seq(double lowest, double highest, double step)
{
  Array<double> result( int( ( highest - lowest) / step ) + 1 , 0);

  int k;
  for (k=0; k<result.size(); k++)
    {
      result[k] = lowest + step * k;
    }

  return result;
}



// this class behaves like vector, except it has
// added functionality to map functions over the
// Array

template <typename R, typename T>
  Array<R> map(const Array<T> & lhs, const Array<T> & rhs, R (*func)(const T & lhs, const T & rhs) );

template <typename R, typename T>
  Array<R> map(const Array<T> & lhs, const T & rhs, R (*func)(const T & lhs, const T & rhs) );

template <typename T, typename R>
  void modify(Array<T> & lhs, const Array<T> & rhs, void (*func)(T & lhs, const T & rhs) );

template <typename T>
OrderedArray<int> suchThat(Array<T> & lhs, const Array<T> & rhs, void (*func)(T & lhs, const T & rhs) );

template <typename T>
ostream & operator << (ostream & os, const Array<T> & rhs);

template <typename T>
istream & operator >> (istream & is, Array<T> & rhs);

// member functions

template <typename T>
Array<T>::Array(int n) : data(n)
{
  
}

template <typename T>
Array<T>::Array(int n, const T & filler) : data(n, filler)
{
  
}

template <typename T>
const T & Array<T>::operator [] (int k) const
{
  boundsCheck(k);
  return data[k];
}

template <typename T>
T & Array<T>::operator [] (int k)
{
  boundsCheck(k);
  return data[k];
}

template <typename T>
void Array<T>::push_back(const T & element)
{
  data.push_back( element );
}

template <typename T>
T & Array<T>::back()
{
  boundsCheck(0);
  return data.back();
}

template <typename T>
const T & Array<T>::back() const
{
  boundsCheck(0);
  return data.back();
}

template <typename T>
int Array<T>::size() const
{
  return data.size();
}

template <typename T>
void Array<T>::boundsCheck(int i) const
{
  if ( i < 0 || i >= size() )
    throw OutOfBoundsException();
}

// nonmember functions

template <typename T, typename R>
Array<R> map(const Array<T> & lhs, const Array<T> & rhs, R (*func)(const T & parL, const T & parR) )
{
  // precondition: lhs and rhs have the same size
  // postcondition: returns an array filled with
  // { ..., func(lhs[i], rhs[i]), ... }

  if ( lhs.size() != rhs.size() )
    {
      throw Array<T>::SizeException();
    }

  Array<R> result;
  
  int k;
  for (k=0; k<lhs.size(); k++)
    {
      result.push_back( func(lhs[k], rhs[k] ) );
    }

  return result;
}

template <typename T, typename R>
Array<R> map(const Array<T> & lhs, const T & rhs, R (*func)(const T & parL, const T & parR) )
{
  // precondition: 
  // postcondition: returns an array filled with
  // { ..., func(lhs[i], rhs), ... }

  Array<R> result;
  
  int k;
  for (k=0; k<lhs.size(); k++)
    {
      result.push_back( func(lhs[k], rhs ) );
    }

  return result;
}

template <typename T, typename R>
void modify(Array<T> & lhs, const Array<T> rhs, void (*func)(T & parL, const T & parL) )
{
  // precondition: lhs and rhs have the same size
  // postcondition: modifies each element of lhs:
  // func( lhs[i], rhs[i] );
}

template <typename T>
OrderedArray<int> suchThat(Array<T> & lhs, bool (*verify)(const T & par) )
{
  // precondition:
  // postcondition: returns a PackedSet containing
  // the indices of lhs for which
  // func(lhs[i]) returns true

  OrderedArray<int> result;

  int k;
  for (k=0; k<lhs.size(); k++)
    {
      if ( verify( lhs[k] ) )
	result.push_back(k);
    }

  return result;
}

template <typename T>
ostream & operator << (ostream & os, const Array<T> & rhs)
{
  os << "{";

  int k;
  for (k=0; k<rhs.size(); k++)
    {
      os << rhs[k];
      if ( k != rhs.size() - 1 )
	cout << ", ";
    }

  os << "}";
  return os;
}

template <typename T>
istream & operator >> (istream & is, Array<T> & rhs);


// this class behaves like Array, except it
// requires a total order
// a < operator is required for T
// 
// objects must be manually inserted in order

template <typename T>
void OrderedArray<T>::push_back(const T & element)
{
  if ( this->size() != 0 && !( this->back() < element ) )
    throw OutOfOrderException();

  this->data.push_back(element);
}

