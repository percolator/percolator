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
 
 $Id: Array.cpp,v 1.2 2008/05/07 21:25:07 lukall Exp $
 
 *******************************************************************************/

template <typename T> class OrderedArray;

#include "OrderedArray.h"

template <typename T>
Array<T>::Array(int n) :
  data(n)
{
}

template <typename T>
Array<T>::Array(int n, const T & element) :
  data(n, element)
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
void Array<T>::clear()
{
  data.clear();
}

template <typename T>
void Array<T>::resize(int n)
{
  if ( n < 0 )
    throw ResizeException();

  data.resize(n);
}

template <typename T>
void Array<T>::resize(int n, const T & element)
{
  if ( n < 0 )
    throw ResizeException();

  data.resize(n,element);
}


template <typename T>
int Array<T>::size() const
{
  return data.size();
}

template <typename T>
bool Array<T>::inBounds(int i) const
{
  return i >= 0 && i < size();
}

template <typename T>
void Array<T>::boundsCheck(int i) const
{
#ifdef SAFE_ARRAYS
  if ( ! inBounds(i) )
    throw OutOfBoundsException();
#endif
}

// nonmember functions

// error checking
template <typename T, typename R>
void sizeCheck(const Array<T> & lhs, const Array<R> & rhs)
{
#ifdef SAFE_ARRAYS
  if ( lhs.size() != rhs.size() )
    {
      throw typename Array<T>::SizeException();
    }
#endif
}

// concatonate arrays together lhs to rhs
template <typename T>
Array<T> concatonate(const Array<T> & lhs, const Array<T> & rhs)
{
  Array<T> result(lhs.size() + rhs.size());
  int k;
  for (k=0; k<lhs.size(); k++)
    {
      result[k] = lhs[k];
    }
  
  for (k=0; k<rhs.size(); k++)
    {
      result[ k + lhs.size() ] = rhs[k];
    }
  return result;
}

// get all values at the column
template <typename T>
Array<T> getCol(const Array<Array<T> > & lhs, const int & col)
{
  Array<T> result( lhs.size() );

  int k;
  for (k=0; k<lhs.size(); k++)
    {
      result[k] = lhs[k][col];
    }

  return result;
}

// index function
template <typename T>
T index(const Array<T> & lhs, const int & rhs)
{
  return lhs[rhs];
}

// indexOf function
template <typename T>
int indexOf(const Array<T> & lhs, const T & rhs)
{
  int k;
  for (k=0; k<lhs.size(); k++)
    {
      if ( lhs[k] == rhs )
	return k;
    }
  return -1;
}

// unary mapping functions

// passed by const reference
template <typename R, typename T>
Array<R> map( R (*func)(const T & lhs), const Array<T> & lhs)
{
  Array<R> result(lhs.size());

  int k;
  for (k=0; k<lhs.size(); k++)
    {
      result[k] = func(lhs[k]);
    }

  return result;
}

// passed by value
template <typename R, typename T>
Array<R> map( R (*func)(T lhs), const Array<T> & lhs)
{
  Array<R> result(lhs.size());

  int k;
  for (k=0; k<lhs.size(); k++)
    {
      result[k] = func(lhs[k]);
    }

  return result;
}


// binary mapping functions

// passed by const reference
template <typename R, typename T, typename S>
Array<R> map(R (*func)(const T & parL, const S & parR), const Array<T> & lhs, const Array<S> & rhs)
{
  // precondition: lhs and rhs have the same size
  // postcondition: returns an array filled with
  // { ..., func(lhs[i], rhs[i]), ... }

  sizeCheck(lhs, rhs);

  Array<R> result(lhs.size());
  
  int k;
  for (k=0; k<lhs.size(); k++)
    {
      result[k] =  func(lhs[k], rhs[k] );
    }

  return result;
}

template <typename R, typename T, typename S>
Array<R> map(R (*func)(const T & parL, const S & parR), const Array<T> & lhs, const S & rhs)
{
  // precondition: 
  // postcondition: returns an array filled with
  // { ..., func(lhs[i], rhs), ... }

  Array<R> result(lhs.size());
  
  int k;
  for (k=0; k<lhs.size(); k++)
    {
      result[k] = func(lhs[k], rhs );
    }

  return result;
}

template <typename R, typename T, typename S>
Array<R> map(R (*func)(const T & parL, const S & parR), const T & lhs, const Array<S> & rhs)
{
  // precondition: 
  // postcondition: returns an array filled with
  // { ..., func(lhs[i], rhs), ... }

  Array<R> result(rhs.size());
  
  int k;
  for (k=0; k<rhs.size(); k++)
    {
      result[k] = func(lhs, rhs[k] );
    }

  return result;
}

// passed by value
template <typename R, typename T, typename S>
Array<R> map(R (*func)(T parL, S parR), const Array<T> & lhs, const Array<S> & rhs)
{
  // precondition: lhs and rhs have the same size
  // postcondition: returns an array filled with
  // { ..., func(lhs[i], rhs[i]), ... }

  sizeCheck(lhs, rhs);

  Array<R> result(lhs.size());
  
  int k;
  for (k=0; k<lhs.size(); k++)
    {
      result[k] = func(lhs[k], rhs[k] );
    }

  return result;
}

template <typename R, typename T, typename S>
Array<R> map(R (*func)(T parL, S parR), const Array<T> & lhs, S rhs)
{
  // precondition: 
  // postcondition: returns an array filled with
  // { ..., func(lhs[i], rhs), ... }

  Array<R> result(lhs.size());
  
  int k;
  for (k=0; k<lhs.size(); k++)
    {
      result[k] = func(lhs[k], rhs );
    }

  return result;
}

template <typename R, typename T, typename S>
Array<R> map(R (*func)(T parL, S parR), T lhs, const Array<S> & rhs)
{
  // precondition: 
  // postcondition: returns an array filled with
  // { ..., func(lhs[i], rhs), ... }

  Array<R> result(rhs.size());
  
  int k;
  for (k=0; k<rhs.size(); k++)
    {
      result[k] = func(lhs, rhs[k]) ;
    }

  return result;
}

// modify functions

// passed by const reference
template <typename T, typename S>
void modify(void (*func)(T & parL, const S & parL), Array<T> & lhs, const Array<S> & rhs)
{
  // precondition: lhs and rhs have the same size
  // postcondition: modifies each element of lhs:
  // func( lhs[i], rhs[i] );

  sizeCheck(lhs, rhs);

  int k;
  for (k=0; k<lhs.size(); k++)
    {
      func( lhs[k], rhs[k] );
    }
}

template <typename T, typename S>
void modify(void (*func)(T & parL, const S & parL), Array<T> & lhs, const S & rhs)
{
  // precondition: lhs and rhs have the same size
  // postcondition: modifies each element of lhs:
  // func( lhs[i], rhs[i] );

  int k;
  for (k=0; k<lhs.size(); k++)
    {
      func( lhs[k], rhs );
    }
}

template <typename T, typename S>
void modify(void (*func)(T & parL, const S & parL), T & lhs, const Array<S> & rhs)
{
  // precondition: lhs and rhs have the same size
  // postcondition: modifies each element of lhs:
  // func( lhs[i], rhs[i] );

  int k;
  for (k=0; k<rhs.size(); k++)
    {
      func( lhs, rhs[k] );
    }
}

// passed by value 
template <typename T, typename S>
void modify(void (*func)(T & parL, S parL), Array<T> & lhs, const Array<S> & rhs)
{
  // precondition: lhs and rhs have the same size
  // postcondition: modifies each element of lhs:
  // func( lhs[i], rhs[i] );

  sizeCheck(lhs, rhs);

  int k;
  for (k=0; k<lhs.size(); k++)
    {
      func( lhs[k], rhs[k] );
    }
}

template <typename T, typename S>
void modify(void (*func)(T & parL, S parL), Array<T> & lhs, S rhs)
{
  // precondition: lhs and rhs have the same size
  // postcondition: modifies each element of lhs:
  // func( lhs[i], rhs[i] );

  int k;
  for (k=0; k<lhs.size(); k++)
    {
      func( lhs[k], rhs );
    }
}

template <typename T, typename S>
void modify(void (*func)(T & parL, S parL), T & lhs, const Array<S> & rhs)
{
  // precondition: lhs and rhs have the same size
  // postcondition: modifies each element of lhs:
  // func( lhs[i], rhs[i] );

  int k;
  for (k=0; k<rhs.size(); k++)
    {
      func( lhs, rhs[k] );
    }
}


// such that functions
template <typename T>
OrderedArray<int> suchThat(bool (*verify)(const T & par), const Array<T> & lhs )
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
OrderedArray<int> suchThat(bool (*verify)(T par), const Array<T> & lhs )
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
  if ( rhs.size() == 0 )
    {
      os << '0';
      return os;
    }

  int w = os.width();
  os << "{";

  int k;
  for (k=0; k<rhs.size(); k++)
    {
      os.width(10);
      os << rhs[k];
      if ( k != rhs.size() - 1 )
	os << ", ";
    }

  os << "}";

  os.width(w);
  return os;
}

template <typename T>
istream & operator >> (istream & is, Array<T> & rhs)
{
  char delim;
  is >> delim;

  rhs.clear();

  // the character for an empty array
  if ( delim == '0')
    return is;

  if ( delim != '{' )
    {
      throw typename Array<T>::ArrayFormatException();
    }
  T element;
  
  // now that a { character has been read
  // it is guaranteed that at least one
  // element will occur

  do
    {
      is >> element;
      rhs.push_back(element);
      is >> delim;
      
      if ( delim != ',' && delim != '}')
	throw typename Array<T>::ArrayFormatException();

    } while ( delim != '}' );
  return is;
}
