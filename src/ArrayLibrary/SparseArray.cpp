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

