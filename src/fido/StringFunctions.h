// Written by Oliver Serang 2009
// see license for more information

#ifndef _StringFunctions_H
#define _StringFunctions_H

#include <iostream>

using namespace std;


class StringFunctions
{
 public:
  static unsigned int simpleStringHash(const string & str)
  {
    unsigned long int i = 0;
    int n = str.length();
    for (int k=0; k < n; k++)
      {
	// appropriate if the characters in this string are late...                                                       
	i += int(str[k]) - int('0');
	i *= 255;
      }
    return int (i << 2);
  }

  static unsigned int nucleotideStringHash(const string & str);
 private:
  
};

#endif

