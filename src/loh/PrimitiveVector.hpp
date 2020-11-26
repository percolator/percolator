#ifndef _PRIMITIVEVECTOR_HPP
#define _PRIMITIVEVECTOR_HPP

#include <utility>
#include <stdlib.h>
#include <assert.h>

template <typename T, unsigned long START_SIZE=4>
class PrimitiveVector {
private:
  unsigned long _size;
  unsigned long _capacity;
  T*__restrict _data;

public:
  PrimitiveVector():
    _size(0),
    _capacity(START_SIZE),
    _data( (T*) calloc(_capacity, sizeof(T)) )
  {
    #ifdef SAFE
    assert(_data != NULL);
    #endif
  }
  
  PrimitiveVector(unsigned long sz):
    _size(sz),
    _capacity(sz),
    _data( (T*) calloc(_capacity, sizeof(T)) )
  {
    #ifdef SAFE
    assert(_data != NULL);
    #endif
  }

  PrimitiveVector(unsigned long sz, const T & fill):
    _size(sz),
    _capacity(sz),
    _data( (T*) calloc(_capacity, sizeof(T)) )
  {
    #ifdef SAFE
    assert(_data != NULL);
    #endif
    for (unsigned long i=0; i<_size; ++i)
      _data[i] = fill;
  }

  PrimitiveVector(const PrimitiveVector & rhs):
    _size(rhs._size),
    _capacity(rhs._capacity),
    _data( (T*) calloc(_capacity, sizeof(T)) )
  {
    #ifdef SAFE
    assert(_data != NULL);
    #endif
    for (unsigned long i=0; i<_size; ++i)
      _data[i] = rhs._data[i];
  }

  PrimitiveVector(PrimitiveVector && rhs):
    _size(0),
    _capacity(0),
    _data(NULL)
  {
    std::swap(_size, rhs._size);
    std::swap(_capacity, rhs._capacity);
    std::swap(_data, rhs._data);
  }
  
  ~PrimitiveVector() {
    free( _data );
  }

  std::pair<T*, unsigned long> liberate_pointer() {
    T*data = _data;
    unsigned long size = _size;

    // Reconstruct:
    _size = 0;
    _capacity = START_SIZE;
    _data = (T*) calloc(_capacity, sizeof(T));
    #ifdef SAFE
    assert(_data != NULL);
    #endif
    
    return {data, size};
  }

  const PrimitiveVector<T> & operator =(const PrimitiveVector & rhs) {
    _size = rhs._size;
    _capacity = rhs._capacity;

    _data = (T*)realloc(_data, _capacity*sizeof(T));
    #ifdef SAFE
    assert(_data != NULL);
    #endif
    for (unsigned long i=0; i<_size; ++i)
      _data[i] = rhs._data[i];

    return *this;
  }

  const PrimitiveVector<T> & operator =(PrimitiveVector && rhs) {
    std::swap(_size, rhs._size);
    std::swap(_capacity, rhs._capacity);
    std::swap(_data, rhs._data);
    return *this;
  }

  void push_back(const T & element) {
    if (_capacity == _size) {
      // Need to resize:
      _capacity += (_capacity>>1) + 1;
      _data = (T*)realloc(_data, _capacity*sizeof(T));
      #ifdef SAFE
      assert(_data != NULL);
      #endif
    }
    
    _data[_size] = element;
    ++_size;
  }

  void resize(const unsigned long new_size) {
    // If it needs to expand, expand
    // else do nothing
    _size = new_size;
    if (_capacity < new_size) {
      _capacity = new_size + 1;
      _data = (T*)realloc(_data, _capacity*sizeof(T));
      #ifdef SAFE
      assert(_data != NULL);
      #endif
    }
  }

  void reserve(const unsigned long new_size) {
    // If it needs to expand, expand
    if (_capacity < new_size) {
      _capacity = new_size + 1;
      _data = (T*)realloc(_data, _capacity*sizeof(T));
      #ifdef SAFE
      assert(_data != NULL);
      #endif
    }
  }

  void reserve_more(const unsigned long growth) {
    reserve(size() + growth);
  }

  const T & back() const {
    return _data[_size-1];
  }
  
  const T & pop_back() {
    --_size;
    return _data[_size];
  }

  const T*begin() const {
    return _data;
  }
  const T*end() const {
    return _data+size();
  }
  T*begin() {
    return _data;
  }
  T*end() {
    return _data+size();
  }
  
  const T & operator [](unsigned long i) const {
    return _data[i];
  }

  T & operator [](unsigned long i) {
    return _data[i];
  }

  void clear() {
    _size = 0;
    _capacity = START_SIZE;
    _data = (T*)realloc(_data, _capacity*sizeof(T));
  }

  unsigned long size() const {
    return _size;
  }
};

#endif
