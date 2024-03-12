// Written by Oliver Serang 2009
// see license for more information

//#ifndef _HASHTABLE_H
//#define _HASHTABLE_H

#include "Array.h"
#include <list>
#include <iostream>

using namespace std;

#include <iterator>
#include <functional>

#define default_size 1001

template <typename D>
class HashTable
{
 public:

 HashTable(unsigned int (*hashFunction) (const D & key), int size = default_size):
  table(size)
  {
    defined_hash = hashFunction;
    numberOfElements = 0;
  }
  virtual ~HashTable()
    {}

  // returns true if the element has been added, false if it is already there
  bool add(const D & data);
  int lookup(const D & d) const;
  const D & getItem(int num);

  int numElements() const
  {
    return numberOfElements;
  }

  const D & operator [] (int k) const
  {
    return itemsByNumber[k];
  }

  const Array<D> & getItemsByNumber() const
  {
    return itemsByNumber;
  }

 private:

  struct Node
  {
    Node(const D & d, int uniq)
    {
      data = d;
      code = uniq;
    }
    bool operator ==(const Node & rhs)
    {
      return data == rhs.data;
    }
    D data;
    int code;
  };

  int searchList(const list<Node> & li, const D & data) const;

  virtual unsigned int hash(const D & key) const;
  unsigned int (*defined_hash)(const D & d);

  int numberOfElements;

  Array<list<Node> > table;
  Array<D> itemsByNumber;
};

#include "HashTable.cpp"

//#endif
