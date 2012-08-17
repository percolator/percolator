// Written by Oliver Serang 2009
// see license for more information

#ifndef _Set_H
#define _Set_H

#include "Array.h"

class Set : public Array<int>
{
 public:
  Set() {}

  void add(int el);
  const Set & operator &=(const Set & rhs);
  const Set & operator |=(const Set & rhs);
  bool isEmpty() const
  {
    return size() == 0;
  }

  int find(int x) const;
  int findHelper(int low, int high, int x) const;

  Set reindexTo(const Set & rhs);
  Set reindexToFind(const Set & rhs);

  class UnsortedOrderException {};
  class InvalidBaseException {};


  Array<int> operator [](const Set & rhs) const
  {
    return Array<int>::operator [](rhs);
  }

  const int & operator [](int k) const
  {
    return Array<int>::operator [](k);
  }

  // for hashing sets
  static unsigned int sumSetElements(const Set & s)
  {
    unsigned int sum = 0;
    for (int k=0; k<s.size(); k++)
      sum += s[k];

    return sum;
  }

  using Array<int>::size;
  using Array<int>::Iterator;
  using Array<int>::begin;
  using Array<int>::end;

  Set without( const Set & rhs ) const;

  static Set FullSet(int low, int high);
  static Set SingletonSet(int value)
  {
    return FullSet(value, value);
  }
  
  int randomElement() const
  {
    return (*this)[ rand() % size() ];
  }

 private:

  bool verify() const;
};

Set operator &(Set lhs, const Set & rhs);
Set operator |(Set lhs, const Set & rhs);


#endif

