#ifndef _SparseArray_H
#define _SparseArray_H

#include <iostream>

using namespace std;

#include "ArrayLibrary.h"

template <typename T>
class SparseArray : public Array<T>
{
public:

  virtual void push_back(int index, const T & element);
  virtual void resize(int n);

  int packedSize() const;

  const SparseArray<T> & operator +=(const SparseArray<T> & rhs);

private:
  OrderedArray<int> nonNull;
};



#endif

