/*
 * Copyright (c) 2008-2009 MacCoss Lab, University of Washington
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */




#ifndef LAZYMATRIX_H_
#define LAZYMATRIX_H_


#include <vector>
#include <deque>
#include <iostream>
#include <assert.h>
#include "crawutils.H"

#ifdef HAVE_CRAWDAD
#include "logging.H"
#endif

using namespace std;

enum sparse_level {
	LM_FULL_DATA,
	LM_RETAIN_DATA,
	LM_SPARSE_DATA
};



class LazyMatrix
{
  friend LazyMatrix * flip_lazy_matrix( LazyMatrix & lm );	
  friend class msmat;
  typedef float value_type;
  typedef vector< vector<value_type> * > data_type;
  typedef LazyMatrix self;	
  //static const int LM_FULL_DATA = 0;
  //static const int LM_RETAIN = 1;
  //static const int LM_SPARSE = 2;
  //typedef value_type * iterator;
  //typedef const value_type* const_iterator;
  /* typedef Value_T* row_type;
  typedef stride_iter<value_type*> col_type;
  typedef const value_type* const_row_type;
  typedef stride_iter<const value_type*> const_col_type;
  */

protected:
	mutable data_type m;
	int nrows;
	int ncols;

	//extremely expedient hack for copied matrices,
	//there is a better way to do this...
	int * data_copies;
	sparse_level sp_level;
	FILE * matrix_fh;
    off_t data_offset;
    void _get_row_lazy ( vector<float> & r_out , int i ) const;
    void _get_row_retain ( vector<float> & r_out, int i );
    void _get_row_full ( vector<float> & r_out, int i ) const;


public:

	LazyMatrix();
  
    	
  LazyMatrix(int r, int c, sparse_level l, FILE * fh, off_t data_offset = 0);

  LazyMatrix( const self & x );

  void copy_m (const self & x);
  void copy_m_flip( const self & x ); 
  LazyMatrix( const self & x, bool flip );

  
  ~LazyMatrix();
  
  void to_full_data(bool force = false);
  void to_full_data(const value_type v, bool force = false);

  void common_init();
	
  void init_data_cntr ();

  void set_row( vector<float> & r_values, int i );

  void set_col( vector<float> & c_values, int j );

  //not const as the retain option can change the object
  void get_row( vector<float> & r_out , int i );


  vector<float> * get_row_p ( int i ) const {
    return m[i];
  }

  void get_col ( vector<float> & c_out, int j );


  void add_val ( float v );
  void mult_val ( float v );
  void div_val ( float v );
  void column_shift ( int shift );
  void linear_shift ( int shift );
  void set_min_val ( float v = 0.0f );
  void set_max_val ( float v );

  inline int rows() const { return nrows; }
  inline int cols() const { return ncols; }
  inline int size() const { return nrows * ncols; }
  inline size_t byte_size() const { return (size_t)( size() * sizeof(data_type)) ; }

  
  
  /*
  row_type row_begin(int n) { return m[n].begin(); }
  row_type row_end(int n) { return m[n].end(); }
  const_row_type row_begin(int n) const { return &m[n * cols( )]; }
  const_row_type row_end(int n) const { return row_begin(n ) + cols(
  ); }
  */
  
  self & operator = (const self & x) { 
    cerr << "LazyMatrix & operator & const" << endl;
    m = x.m; nrows = x.nrows; 
    ncols = x.ncols; sp_level = x.sp_level; return *this; 
  }

  vector<value_type> * at(int i) const {   return m.at(i); };
  //vector<value_type> * operator[](int i) const { return m[i]; };
  vector<value_type> * operator()(int i) const { return m[i]; };
  inline float operator() (int i, int j) const { 
    return (*(m[i]))[j];
  }
  inline float get_val(int i, int j) const { 
    return (*this)(i,j);
  }
  void set_val( int i , int j, float v) {
	  (*(m[i]))[j] = v;
  }

  const vector<value_type> * operator[](int n) const { return m[n]; }

  void remove_start_rows(int n );

  void remove_end_rows(int n);
  void remove_rows( int lh_idx, int rh_idx);

  
  void remove_start_cols( int n );

  void remove_stop_cols( int n );
  void remove_cols ( int start_idx, int stop_idx);

  /*self& operator+=(const self& x) { m += x.m; return *this; }
  self& operator-=(const self& x) { m -= x.m; return *this; }
  self& operator+=(value_type x) { m += x; return *this; }
  self& operator-=(value_type x) { m -= x; return *this; }
  self& operator*=(value_type x) { m *= x; return *this; }
  self& operator/=(value_type x) { m /= x; return *this; }
  self& operator%=(value_type x) { m %= x; return *this; }
  self operator-( ) { return -m; }
  self operator+( ) { return +m; }
  self operator!( ) { return !m; }
  self operator~( ) { return ~m; }
  */

  // friend operators
  /*friend self operator+(const self& x, const self& y) { return self(x) += y; }
  friend self operator-(const self& x, const self& y) { return self(x) -= y; }
  friend self operator+(const self& x, value_type y) { return self(x) += y; }
  friend self operator-(const self& x, value_type y) { return self(x) -= y; }
  friend self operator*(const self& x, value_type y) { return self(x) *= y; }
  friend self operator/(const self& x, value_type y) { return self(x) /= y; }
  friend self operator%(const self& x, value_type y) { return self(x) %= y; }
  */	
};


///A Lazy Matrix that is meant to retain rows recently accessed. Not fully implemented.
class BufferedLazyMatrix : public LazyMatrix {
  private :
    int top_populated_idx;
    int bottom_populated_idx;
public :
  //TODO -- populate the functions below for file reading usage 
  void common_init();
  void get_row( vector<float> & r_out , int i );
  
  
};

LazyMatrix * flip_lazy_matrix( LazyMatrix & lm );


#endif /*LAZYMATRIX_H_*/
