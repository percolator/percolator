// Written by Oliver Serang 2009
// see license for more information

#ifndef _StringFunctions_H
#define _StringFunctions_H

#include <iostream>

using namespace std;

class StringFunctions
{
public:
  static unsigned int simpleStringHash(const string &str)
  {
    unsigned long int i = 0;
    size_t n = str.length();
    for (size_t k = 0; k < n; k++)
    {
      // appropriate if the characters in this string are late...
      i += static_cast<unsigned long int>(int(str[k]) - int('0'));
      i *= 255;
    }
    unsigned int result = static_cast<unsigned int>(i << 2);
    return result;
  }

  static unsigned int nucleotideStringHash(const string &str);

private:
};

#endif
