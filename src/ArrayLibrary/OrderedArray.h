// this class behaves like Array, except it
// requires a total order
// a < operator is required for T
// 
// objects must be manually inserted in order

#ifndef _OrderedArray_H
#define _OrderedArray_H

using namespace std;

template <class T>
class OrderedArray : public Array<T>
{
public:
  // constructors
  OrderedArray() {}
  OrderedArray(int n);

  virtual void push_back(const T & element);

  class OutOfOrderException {};
protected:
  
};

#endif

