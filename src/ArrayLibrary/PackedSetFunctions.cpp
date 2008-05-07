/*******************************************************************************
 Copyright (c) 2008 Oliver Serang

 Permission is hereby granted, free of charge, to any person
 obtaining a copy of this software and associated documentation
 files (the "Software"), to deal in the Software without
 restriction, including without limitation the rights to use,
 copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the
 Software is furnished to do so, subject to the following
 conditions:

 The above copyright notice and this permission notice shall be
 included in all copies or substantial portions of the Software. 
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 OTHER DEALINGS IN THE SOFTWARE.
 
 $Id: PackedSetFunctions.cpp,v 1.2 2008/05/07 21:25:07 lukall Exp $
 
 *******************************************************************************/
#include "PackedSetFunctions.h"

PackedSet intersection(const PackedSet & lhs, const PackedSet & rhs)
{
  PackedSet result;

  int kL, kR;
  for (kL=0, kR=0; kL < lhs.size() && kR < rhs.size(); )
    {
      int iL = lhs[kL];
      int iR = rhs[kR];

      if ( iL < iR )
	{
	  result.push_back(iL);
	  kL++;
	}
      else if ( iL > iR )
	{
	  result.push_back(iR);
	  kR++;
	}
      else
	{
	  // equality case
	  result.push_back(iL);
	  kL++;
	  kR++;
	}
    }

  // add anything left over
  if ( kL < lhs.size() )
    {
      for ( ; kL<lhs.size(); kL++ )
	{
	  result.push_back( lhs[kL] );
	}
    }
  if ( kR < rhs.size() )
    {
      for ( ; kR<rhs.size(); kR++ )
	{
	  result.push_back( rhs[kR] );
	}
    }

  return result;
}

Array<int> intersection(const Array<int> & lhs, const Array<int> & rhs)
{
  cout << "In inter..." << endl;

  Array<int> result;

  if ( lhs.size() == 0 || rhs.size() == 0 )
    return result;

  int big = max(lhs.back(), rhs.back());

  Array<int> set(big, 1);

  int k;
  cout << "\t\tinsert lhs" << endl;
  for (k=0; k<lhs.size(); k++)
    {
      set[ lhs[k] ] = 0;
    }

  cout << "\t\tinsert rhs" << endl;
  for (k=0; k<rhs.size(); k++)
    {
      set[ rhs[k] ] = 0;
    }

  cout << "\t\textract" << endl;
  for (k=0; k<set.size(); k++)
    {
      if ( ! set[k] ) 
	result.push_back( k );
    }

  cout << "Done inter" << endl;
  return result;
}

