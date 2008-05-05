#ifndef _PackedArray_H
#define _PackedArray_H

#include "OrderedArray.h"

using namespace std;

template <typename T>
struct PackedPair
{
  PackedPair() {}
  PackedPair(int i, const T & v)
  {
    index = i;
    value = v;
  }
  int index;
  T value;
};

operator < (const PackedPair & lhs, const PackedPair & rhs)
{
  return lhs.index < rhs.index;
}


template <typename T>
class PackedArray : public OrderedArray< PackedPair<T> >
{
public:
  virtual const T & operator [] (int k) const;
  virtual T & operator [] (int k);
  virtual void push_back(const PackedPair<T> & element);
protected:

};

#include "PackedArray.cpp"

#endif

