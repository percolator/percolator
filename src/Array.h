// Written by Oliver Serang 2009
// see license for more information

// this class behaves like vector, except it has
// added functionality to map functions over the
// Array
#ifndef _Array_H
#define _Array_H


#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <sstream>
#include <algorithm>
#include <functional>

//class Set;
using namespace std;



template <typename T>
class Array 
{
public:
  // constructors
  Array(int n) : data(static_cast<std::size_t>(n)){};
  Array(std::size_t n) : data(n) {};
  Array(int n, const T & element) : data(std::size_t(n), element){};
  Array(const vector<T> & newvector) : data(newvector){};
  Array():data() {};

  // explicit Array(int n);
  // explicit Array(std::size_t n);
  // Array(int n, const T & defaultValue);
  // Array(const vector<T> & newvector);
  // destructor
  virtual ~Array() {};

  // important vector functions
  virtual const T & operator [] (int k) const;
  virtual T & operator [] (int k);
  virtual const T & operator [] (std::size_t k) const;
  virtual T & operator [] (std::size_t k);
  virtual Array<T> operator [](const Array<int> & rhs) const;

  //Set operator ==(const T & rhs) const;
  bool operator ==(const Array<T> & rhs) const;

  Array<int> sort();
  Array<int> sortA();
  
  vector<T> getVector();
  vector<T> getVector() const;

  virtual void add(const T & element);
  virtual void append(const Array<T> & elements);

  virtual T & back();
  virtual const T & back() const;

  virtual void clear();
  virtual void resize(int n);
  
  void remove(int k);

  // accessors
  virtual std::size_t size() const;

  // bound checking function
  bool inBounds(int i) const;

  class Iterator {
      //  protected:
   public:
      const Array<T>*array;
      int location;
   public:
      Iterator() : array(nullptr), location(-1) {}

      Iterator(const Array<T>* a, int loc) : array(a), location(loc) {}

      int getLocation() const {
        return location;
      }

      const Iterator & operator++(int) {
        location++;
        return *this;
      }

      bool operator!=(const Iterator & rhs) const {
        return !(*this == rhs);
      }

      bool operator==(const Iterator & rhs) const {
        return array == rhs.array && location == rhs.location;
      }
      bool operator <(const Iterator & rhs) {
        return array == rhs.array && location < rhs.location;
      }
      const T & operator *() const {
        return (*array)[location];
      }
  };

  Iterator begin() const {
    return Iterator(this, 0);
  }

  Iterator end() const {
    return Iterator(this, static_cast<int>(size()));
  }

  // exception classes
  class OutOfBoundsException {};
  class SizeException {};
  class FormatException {};
  class ResizeException {};
  
protected:
  // error checking functions
  void boundsCheck(int i) const;

  vector<T> data;
};

// non-member error checking functions
template <typename T, typename R>
  void sizeCheck(const Array<T> & lhs, const Array<R> & rhs);

// concatenate
template <typename T>
Array<T> concatonate(const Array<T> & lhs, const Array<T> & rhs);

template <typename T>
ostream & operator <<(ostream & os, const Array<T> & rhs);

template <typename T>
ostream & operator >>(ostream & os, Array<T> & rhs);

//#include "Set.h"
//#include "Array.cpp"


// member functions
class Set;

template <typename T>
const T & Array<T>::operator [] (int k) const
{
  boundsCheck(k);
  return data[static_cast<std::size_t>(k)];
}

template <typename T>
T & Array<T>::operator [] (int k)
{
  boundsCheck(k);
  return data[static_cast<std::size_t>(k)];
}

template <typename T>
const T & Array<T>::operator [] (std::size_t k) const
{
  boundsCheck(static_cast<int>(k));
  return data[k];
}

template <typename T>
T & Array<T>::operator [] (std::size_t k)
{
  boundsCheck(static_cast<int>(k));
  return data[k];
}

template <typename T>
void Array<T>::add(const T & element)
{
  data.push_back( element );
}

template <typename T>
void Array<T>::append(const Array<T> & elements)
{
  for (std::size_t k=0; k<elements.size(); k++)
    {
      add(elements[k]);
    }
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
void Array<T>::remove(int k)
{
  boundsCheck(k);
  
  int j;
  for (j=k; j<size()-1; j++)
    {
      data[j] = data[j+1];
    }
  
  resize( size() - 1 );
}


template <typename T>
void Array<T>::resize(int n)
{
  if ( n < 0 )
    throw ResizeException();

  data.resize(static_cast<std::size_t>(n));
}

template <typename T>
size_t Array<T>::size() const
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

template <typename T>
Array<T> Array<T>::operator [](const Array<int> & rhs) const
{
  Array<T> result(static_cast<int>(rhs.size()));

  int counter = 0;
  for (Array<int>::Iterator iter = rhs.begin(); iter != rhs.end(); iter++, counter++)
    {
      result[ counter ] = (*this)[ *iter ];
    }

  return result;
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

template <typename T>
ostream & operator <<(ostream & os, const Array<T> & rhs)
{
  if ( rhs.size() == 0 )
    {
      os << '0';
      return os;
    }

  int w = static_cast<int>(os.width());
  os << "{ ";

  std::size_t k;
  for (k=0; k<rhs.size(); k++)
    {
      os.width(4);
      os << rhs[k];
      if ( k != rhs.size() - 1 )
	os << " , ";
    }

  os << " }";

  os.width(w);
  return os;
}

template <typename T>
istream & operator >>(istream & is, Array<T> & rhs)
{
  char delim;
  is >> delim;

  rhs.clear();

  // the character for an empty array
  if ( delim == '0')
    return is;

  if ( delim != '{' )
    {
      throw typename Array<T>::FormatException();
    }
  T element;
  
  // now that a { character has been read
  // it is guaranteed that at least one
  // element will occur

  do
    {
      is >> element;
      rhs.add(element);
      is >> delim;
      
      if ( delim != ',' && delim != '}')
	throw typename Array<T>::FormatException();

    } while ( delim != '}' );
  return is;
}

// template <typename T>
// Set Array<T>::operator ==(const T & rhs) const
// {
//   Set result;
// 
//   for (int k=0; k<size(); k++)
//     {
//       if ( (*this)[k] == rhs )
// 	{
// 	  result.add(k);
// 	}
//     }
// 
//   return result;
// }

template <typename T>
bool Array<T>::operator ==(const Array<T> & rhs) const {
  if ( size() != rhs.size() )
    return false;

  for (int k=0; k<size(); k++) {
    if ( (*this)[k] != rhs[k] ) {
  	  return false;
	  }
  }
  return true;
}

template <typename T>
Array<int> Array<T>::sort() {
  vector<pair<T, int> > sortie(size());
  for (std::size_t k = 0; k < size(); k++) {
    sortie[k] = pair<T, int>( (*this)[k], k);
  }

  ::sort( sortie.begin(), sortie.end() , std::greater<std::pair<T, int> >() );
  
  Array<int> result(static_cast<int>(size()));
  for (std::size_t k = 0; k < size(); k++) {
    (*this)[k] = sortie[k].first;
    result[k] = sortie[k].second;
  }
  return result;
}

template <typename T>
Array<int> Array<T>::sortA() {
  vector<pair<T, int> > sortie(size());

  for (int k = 0; k < size(); k++) {
    sortie[k] = pair<T, int>( (*this)[k], k);
  }

  ::sort(sortie.begin(), sortie.end());
  
  Array<int> result(size());
  for (std::size_t k = 0; k < size(); k++) {
    (*this)[k] = sortie[k].first;
    result[k] = sortie[k].second;
  }

  return result;
}

template <typename T>
vector<T> Array<T>::getVector() {
  return data;
}

template <typename T>
vector<T> Array<T>::getVector() const {
  return data;
}


#endif

