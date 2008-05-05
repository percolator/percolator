#include "PackedVec.h"

double norm(const PackedVec & rhs)
{
  return sqrt(pDot(rhs, rhs));
}

double pDot(const PackedVec & lhs, const Vec & rhs)
{
  double tot = 0;

  int k;
  for (k=0; k<lhs.packedSize(); k++)
    {
      tot += lhs[k] * rhs[ lhs.index(k) ];
    }

  return tot;
}

double pDot(const Vec & lhs, const PackedVec & rhs)
{
  return pDot(rhs, lhs);
}

double pDot(const PackedVec & lhs, const PackedVec & rhs)
{
  double tot = 0;
  
  int kL, kR;
  for (kL=0, kR=0; kL < lhs.packedSize() && kR < rhs.packedSize(); )
    {
      int iL = lhs.index(kL);
      int iR = rhs.index(kR);

      if ( iL < iR )
	{
	  kL++;
	}
      else if ( iL > iR )
	{
	  kR++;
	}
      else
	{
	  // equality case
	  tot += lhs[ kL ] * rhs[ kR ];
	  kL++;
	  kR++;
	}
    }
  return tot;
}

const Vec & operator += (Vec & lhs, const PackedVec & rhs)
{
  int k;
  for (k=0; k<rhs.packedSize(); k++)
    {
      lhs[ rhs.index(k) ] += rhs[k];
    }

  return lhs;
}

const Vec & operator -= (Vec & lhs, const PackedVec & rhs)
{
  int k;
  for (k=0; k<rhs.packedSize(); k++)
    {
      lhs[ rhs.index(k) ] -= rhs[k];
    }

  return lhs;
}

const Vec & operator *= (Vec & lhs, const PackedVec & rhs)
{
  int k;
  for (k=0; k<rhs.packedSize(); k++)
    {
      lhs[ rhs.index(k) ] *= rhs[k];
    }

  return lhs;
}

const Vec & operator /= (Vec & lhs, const PackedVec & rhs)
{
  int k;
  for (k=0; k<rhs.packedSize(); k++)
    {
      lhs[ rhs.index(k) ] /= rhs[k];
    }

  return lhs;
}

const PackedVec & operator /= (PackedVec & lhs, const double & rhs)
{
  int k;
  for (k=0; k<lhs.packedSize(); k++)
  {
    lhs[ k ] /= rhs;
  }

  return lhs;
}



PackedVec operator + (const PackedVec & lhs, const PackedVec & rhs)
{
  PackedVec result;

  int kL, kR;
  for (kL=0, kR=0; kL < lhs.packedSize() && kR < rhs.packedSize(); )
    {
      int iL = lhs.index(kL);
      int iR = rhs.index(kR);

      if ( iL < iR )
	{
	  result.push_back(iL, lhs[kL]);
	  kL++;
	}
      else if ( iL > iR )
	{
	  result.push_back(iR, rhs[kR]);
	  kR++;
	}
      else
	{
	  // equality case
	  double res = lhs[kL]+rhs[kR];
	  if ( Numerical::isNonzero(res) )
	    result.push_back(iL, res);
	  kL++;
	  kR++;
	}
    }
  // add anything left over
  if ( kL < lhs.packedSize() )
    {
      for ( ; kL<lhs.packedSize(); kL++ )
	{
	  result.push_back( lhs.index(kL), lhs[kL] );
	}
    }
  if ( kR < rhs.packedSize() )
    {
      for ( ; kR<rhs.packedSize(); kR++ )
	{
	  result.push_back( rhs.index(kR), rhs[kR] );
	}
    }

  return result;
}

PackedVec operator - (const PackedVec & lhs, const PackedVec & rhs)
{
  PackedVec result;

  int kL, kR;
  for (kL=0, kR=0; kL < lhs.packedSize() && kR < rhs.packedSize(); )
    {
      int iL = lhs.index(kL);
      int iR = rhs.index(kR);

      if ( iL < iR )
	{
	  result.push_back(iL, lhs[kL]);
	  kL++;
	}
      else if ( iL > iR )
	{
//      result.push_back(iR, rhs[kR]); LK
      result.push_back(iR, -rhs[kR]);
	  kR++;
	}
      else
	{
	  // equality case
	  double res = lhs[kL] - rhs[kR];
	  if ( Numerical::isNonzero(res) )
	    result.push_back(iL, res);
	  kL++;
	  kR++;
	}
    }
  // add anything left over
  if ( kL < lhs.packedSize() )
    {
      for ( ; kL<lhs.packedSize(); kL++ )
	{
	  result.push_back( lhs.index(kL), lhs[kL] );
	}
    }
  if ( kR < rhs.packedSize() )
    {
      for ( ; kR<rhs.packedSize(); kR++ )
	{
	  result.push_back( rhs.index(kR), -rhs[kR] );
	}
    }

  return result;
}

const PackedVec & operator += (PackedVec & lhs, const PackedVec & rhs)
{
  lhs = lhs + rhs;
  return lhs;
}

const PackedVec & operator -= (PackedVec & lhs, const PackedVec & rhs)
{
  lhs = lhs - rhs;
  return lhs;
}

PackedVec operator *(double d, const PackedVec & rhs)
{
  PackedVec result = rhs;

  for (int k=0; k<rhs.packedSize(); k++)
    {
      result[k] *= d;
    }
  
  return result;
}

istream & operator >>(istream & is, PackedVec & rhs)
{
  char delim;

  is >> delim;

  rhs.clear();
  if ( delim == '0')
    return is;

  if ( delim != '{' )
    throw PackedVec::PackedVecFormatException();
  
  int i;
  double d;
  for ( ; ; )
    {
      is >> i >> delim >> d;

      rhs.push_back(i, d);
      if ( delim != ':' )
	throw PackedVec::PackedVecFormatException();

      is >> delim;

      if ( delim == '}' )
	break;

      if ( delim != ',' )
	throw PackedVec::PackedVecFormatException();
    }
  return is;
}

ostream & operator <<(ostream & os, const PackedVec & rhs)
{
  if (rhs.packedSize() == 0)
    {
      os << "0";
      return os;
    }

  os << "{ ";
  for (int k=0; k<rhs.packedSize(); k++)
    {
      os << rhs.index(k) << ":" << rhs[k];
      if ( k != rhs.packedSize() - 1 )
	os << ", ";
    }
  os << " }";
  
  return os;
}

// projection functions
PackedVec projectSingle(const PackedVec & basis, const PackedVec & v)
{
  return ( pDot(basis, v)/pDot(basis, basis) ) * basis;
}

PackedVec projectSingle(const PackedVec & basis, const Vec & v)
{
  return ( pDot(basis, v)/pDot(basis, basis) ) * basis;
}

