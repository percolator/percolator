// this class behaves like vector, except it has
// added functionality to map functions over the
// Array

#ifndef _Array_H
#define _Array_H

template <typename T>
class Array
{
public:
  // constructors
  Array():data() {}
  virtual ~Array() {}

  Array(int n);
  Array(int n, const T & defaultValue);

  // important vector functions
  virtual const T & operator [] (int k) const;
  virtual T & operator [] (int k);
  virtual void push_back(const T & element);

  virtual T & back();
  virtual const T & back() const;

  virtual void clear();
  virtual void resize(int n);
  virtual void resize(int n, const T & defaultValue);

  // accessors
  virtual int size() const;

  // novel functions
  bool inBounds(int i) const;

  // exception classes
  class OutOfBoundsException {};
  class SizeException {};
  class ArrayFormatException {};
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

// get column
template <typename T>
Array<T> getCol(const Array<Array<T> > & lhs, const int & col);

// binary index function for the map function
template <typename T>
T index(const Array<T> & lhs, const int & rhs);

// indexOf function to get the location of an item
template <typename T>
int indexOf(const Array<T> & lhs, const T & rhs);

// unary mapping functions
// distribute a unary function call over an array
template <typename R, typename T>
  Array<R> map( R (*func)(const T & lhs), const Array<T> & lhs);

template <typename R, typename T>
  Array<R> map( R (*func)(T lhs), const Array<T> & lhs);


// binary mapping functions
// distribute a binary function call over two arrays
// or over an array and an argument

// passed by const reference
template <typename R, typename T, typename S>
  Array<R> map( R (*func)(const T & lhs, const S & rhs), const Array<T> & lhs, const Array<S> & rhs);

template <typename R, typename T, typename S>
  Array<R> map( R (*func)(const T & lhs, const S & rhs), const Array<T> & lhs, const S & rhs);

template <typename R, typename T, typename S>
  Array<R> map(R (*func)(const T & parL, const S & parR), const T & lhs, const Array<S> & rhs);

// passed by value
template <typename R, typename T, typename S>
  Array<R> map( R (*func)(T parL, S parR), const Array<T> & lhs, const Array<S> & rhs);

template <typename R, typename T, typename S>
  Array<R> map(R (*func)(T parL, S parR), const Array<T> & lhs, S rhs);

template <typename R, typename T, typename S>
  Array<R> map(R (*func)(T parL, S parR), T lhs, const Array<S> & rhs);

// modify functions
// these are the same as mapping functions, but they
// modify the contents of the arguments

// passed by const reference
template <typename T, typename S>
  void modify( void (*func)(T & lhs, const S & rhs), Array<T> & lhs, const Array<S> & rhs);

template <typename T, typename S>
  void modify( void (*func)(T & lhs, const S & rhs), Array<T> & lhs, const S & rhs);

template <typename T, typename S>
  void modify( void (*func)(T & lhs, const S & rhs), T & lhs, const Array<S> & rhs);

// passed by value
template <typename T, typename S>
  void modify( void (*func)(T & lhs, S rhs), const Array<T> & lhs, const Array<S> & rhs);

template <typename T, typename S>
  void modify( void (*func)(T & lhs, S rhs), const Array<T> & lhs, S rhs);

template <typename T, typename S>
  void modify( void (*func)(T & lhs, S rhs), T & lhs, const Array<S> & rhs);

// boolean such that
// returns an ordered array of indices
// where evaluating the funciton on that
// element returns true

// passed by reference
template <typename T>
OrderedArray<int> suchThat( bool (*func)(const T & rhs), const Array<T> & rhs);

// passed by value
template <typename T>
OrderedArray<int> suchThat( bool (*func)(T rhs), const Array<T> & rhs);

// stream output and input functions
template <typename T>
ostream & operator << (ostream & os, const Array<T> & rhs);

template <typename T>
istream & operator >> (istream & is, Array<T> & rhs);

#endif

