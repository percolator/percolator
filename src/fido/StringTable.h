// Written by Oliver Serang 2009
// see license for more information

#ifndef _StringTable_H
#define _StringTable_H

#include <iostream>
#include <string>
#include "HashTable.h"
#include "StringFunctions.h"

using namespace std;

class StringTable : public HashTable<string>
{
public:
 StringTable() :
  HashTable<string>(StringFunctions::simpleStringHash)
  {
  }
 StringTable(int size) :
  HashTable<string>(StringFunctions::simpleStringHash, size)
  {
  }

  static StringTable AddElements(const Array<string> & elements)
  {
    StringTable result;
    for (int k=0; k<elements.size(); k++)
      {
	result.add(elements[k]);
      }

    return result;
  }


private:
};

#endif

