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
#include "Streamable.h"

//class Set;
using namespace std;



template <typename T>
class Array 
{
public:
  // constructors
 Array():data() {}

  explicit Array(int n);
  Array(int n, const T & defaultValue);
  Array(const vector<T> & newvector);
  // destructor
  virtual ~Array() {};

  // important vector functions
  virtual const T & operator [] (int k) const;
  virtual T & operator [] (int k);
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
  virtual int size() const;

  // bound checking function
  bool inBounds(int i) const;

  class Iterator
  {
    //  protected:
  public:
    const Array<T>*array;
    int location;
  public:
    Iterator()
      {
	array = NULL;
	location = -1;
      }
    Iterator(const Array<T>*a, int loc)
      {
	array = a;
	location = loc;
      }
    int getLocation() const
    {
      return location;
    }
    const Iterator & operator ++(int)
      {
	location++;
	return *this;
      }
    bool operator !=(const Iterator & rhs)
    {
      return ! ((*this) == rhs);
    }
    bool operator ==(const Iterator & rhs)
    {
      return array == rhs.array && location == rhs.location;
    }
    bool operator <(const Iterator & rhs)
    {
      return array == rhs.array && location < rhs.location;
    }
    const T & operator *() const
    {
      return (*array)[location];
    }
  };

  Iterator begin() const
  {
    return Iterator(this, 0);
  }

  Iterator end() const
  {
    return Iterator(this, size());
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

// concatonate
template <typename T>
Array<T> concatonate(const Array<T> & lhs, const Array<T> & rhs);

template <typename T>
ostream & operator <<(ostream & os, const Array<T> & rhs);

template <typename T>
ostream & operator >>(ostream & os, Array<T> & rhs);

//#include "Set.h"
#include "Array.cpp"

#endif

