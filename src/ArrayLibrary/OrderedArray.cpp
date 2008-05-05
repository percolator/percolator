#include "OrderedArray.h"

template <typename T>
OrderedArray<T>::OrderedArray(int n) :
  Array<T>(n)
{ 
}

template <typename T>
void OrderedArray<T>::push_back(const T & element)
{
#ifdef SAFE_ARRAYS
  if ( Array<T>::size() != 0 && !( Array<T>::back() < element ) )
    throw OutOfOrderException();
#endif

  Array<T>::push_back(element);
}

