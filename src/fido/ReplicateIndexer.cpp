// Written by Oliver Serang 2009
// see license for more information

template <typename T>
Array<Set> ReplicateIndexer<T>::replicates(unsigned int (*hashFunction) (const T & key), const Array<T> & rhs) {
  Array<Set> repSets;
  HashTable<T> ht(hashFunction, 10001);
      
  for (int k = 0; k < rhs.size(); k++) {
    if ( ht.add(rhs[k]) ) {
      repSets.add( Set::SingletonSet(k) );
    } else {
      repSets[ ht.lookup(rhs[k]) ].add(k);
    }
  }

  return repSets;
}

