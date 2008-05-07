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
 
 $Id: SparseArray.cpp,v 1.2 2008/05/07 21:25:07 lukall Exp $
 
 *******************************************************************************/
template <typename T>
void SparseArray<T>::push_back(int index, const T & element)
{
  nonNull.push_back(index);

  if ( size() <= index )
    resize(index+1);

  (*this)[index] = element;
}

template <typename T>
int SparseArray<T>::packedSize() const
{
  return nonNull.size();
}

template <typename T>
void SparseArray<T>::resize(int n)
{
  if ( packedSize() == 0 )
    return;

  int k;
  for (k=0; k<n; k++)
    {
      if ( nonNull[k] >= n )
	break;
    }

  // remove everything in the packed vector after and at k
  nonNull.resize(k);

  Array<T>::resize(n);
}

const SparseArray<T> & SparseArray<T>::operator +=(const SparseArray<T> & rhs)
{
  sizeCheck(*this, rhs);

  int k;
  nonNull = merge(nonNull, rhs.nonNull);
  for (k=0; k<packedSize(); k++)
    {
      (*this)[ nonNull[k] ] += rhs[ nonNull[k] ];
    }

  return *this;
}

const SparseArray<T> & SparseArray<T>::operator +=(const SparseArray<T> & rhs)
{
  sizeCheck(*this, rhs);

  int k;
  nonNull = merge(nonNull, rhs.nonNull);
  for (k=0; k<packedSize(); k++)
    {
      (*this)[ nonNull[k] ] += rhs[ nonNull[k] ];
    }

  return *this;
}

const SparseArray<T> & SparseArray<T>::operator +=(const SparseArray<T> & rhs)
{
  sizeCheck(*this, rhs);

  int k;
  nonNull = merge(nonNull, rhs.nonNull);
  for (k=0; k<packedSize(); k++)
    {
      (*this)[ nonNull[k] ] += rhs[ nonNull[k] ];
    }

  return *this;
}

const SparseArray<T> & SparseArray<T>::operator +=(const SparseArray<T> & rhs)
{
  sizeCheck(*this, rhs);

  int k;
  nonNull = merge(nonNull, rhs.nonNull);
  for (k=0; k<packedSize(); k++)
    {
      (*this)[ nonNull[k] ] += rhs[ nonNull[k] ];
    }

  return *this;
}

