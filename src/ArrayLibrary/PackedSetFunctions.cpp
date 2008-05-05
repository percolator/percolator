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

