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
#include "LazyMatrix.H"


LazyMatrix * flip_lazy_matrix( LazyMatrix & lm ) {
	LazyMatrix * new_lm = new LazyMatrix(lm.cols(), lm.rows(), LM_FULL_DATA , lm.matrix_fh );
	for ( uint j = 0 ; j < lm.rows() ; j++ ) {
		for ( uint i = 0 ; i < lm.cols() ; i++ ) {
			new_lm->set_val(i,j, lm(j,i) );
		}
	}
	return new_lm;
};

LazyMatrix::LazyMatrix() : m(), nrows(0), ncols(0), 
			   sp_level(LM_FULL_DATA), matrix_fh(NULL) {
  common_init();
}
  
    	
LazyMatrix::LazyMatrix(int r, int c, sparse_level l, FILE * fh, off_t data_start_offset) : m(r), nrows(r),
								  ncols(c), sp_level(l), matrix_fh(fh), data_offset(data_start_offset) {
  common_init();
  if ( sp_level == LM_FULL_DATA ) {
    for ( int i = 0 ; i < nrows ; i++ ) {
      m[i] = new vector<value_type>(ncols);
    }	
  }
  else {
    for ( int i = 0 ; i < nrows; i++ ) {
      m[i] = NULL;
    }
  } 
};

LazyMatrix::LazyMatrix( const self & x ) : m(x.m), nrows(x.nrows), ncols(x.ncols),
              sp_level(x.sp_level) , matrix_fh(x.matrix_fh),  data_offset(x.data_offset) {
     common_init();
     //copy the pointer, not the value
     data_copies = x.data_copies;
     //increment ptr value, both objects now have +1 ref counts
     *data_copies++;
}

  void LazyMatrix::copy_m (const self & x) {
    for ( uint i = 0 ; i < x.m.size() ; i++ ) {
      copy( (x.m[i])->begin() , (x.m[i])->end() , 
	    m[i]->begin() );
    }
  }

  void LazyMatrix::copy_m_flip( const self & x ) {
      for ( uint i = 0 ; i < x.nrows ; i++ ) {
  
      }
  }
/* gad, this should be a clone method */
  LazyMatrix::LazyMatrix( const self & x, bool flip ) :
	  m(x.m.size()),
    nrows(x.nrows), ncols(x.ncols), sp_level(x.sp_level), matrix_fh(x.matrix_fh)  , data_offset(x.data_offset)
  {
    common_init();
    *data_copies = 1;
    to_full_data(/*force*/ true);
    copy_m(x);
 
  }
  
  LazyMatrix::~LazyMatrix() {
	  int dc = *(data_copies) - 1;
	  *data_copies = dc;
	  if ( *data_copies <= 0 ) {
	    for ( uint i = 0 ; i < m.size() ; i++ ) {
	      if ( m[i] != NULL ) {
		delete m[i];
		m[i] = NULL;
	      }
	      
	    }
	  }  
  };
  
  //todo bool flag for being in full data state
  void LazyMatrix::to_full_data(bool force) {
      if ( force ||  !( sp_level == LM_FULL_DATA) ) { 
          sp_level = LM_FULL_DATA;
          for ( int i = 0 ; i < m.size() ; i++ ) {
              if ( m[i] == NULL ) {
                  m[i] = new vector<value_type>(cols());
                  //TODO -- actually read from file here...
              }
              else if ( m[i]->size() == 0 ) {
                  m[i]->resize(cols());
              }
          }
      }
  }
void LazyMatrix::to_full_data(const value_type v, bool force) {
      if ( force || !( sp_level == LM_FULL_DATA) ) { 
	sp_level = LM_FULL_DATA;
	  for ( int i = 0 ; i < m.size() ; i++ ) {
	    if ( m[i] == NULL ) {
	      m[i] = new vector<value_type>(cols(),v);
              //TODO -- actually read from file here...
	    }
	    else if ( m[i]->size() == 0 ) {
              m[i]->resize(cols(),v);
	    }
	  }
	  
      }
  }
  //todo test reference counting, remove it if not necessary
  void LazyMatrix::common_init() {
  	 init_data_cntr();
  }
	
  void LazyMatrix::init_data_cntr () {
  	 data_copies = new int;
  	 *data_copies = 1;
  }

  void LazyMatrix::set_row( vector<float> & r_values, int i ) {
	  assert(r_values.size() == cols());
	  if ( m[i] == NULL ) {
		  m[i] = new vector<float>(r_values.size());
	  }
	  copy(r_values.begin(),r_values.end(),m[i]->begin());
  }

  void LazyMatrix::set_col( vector<float> & c_values, int j ) {
	  if ( sp_level != LM_FULL_DATA ){
		  cerr << "converting matrix to full_data mode" << endl;
		  to_full_data();
	  }
	  for ( int i = 0 ; i < nrows ; i++ ) {
		  set_val(i,j,c_values[i]);
	  }
  }

  void LazyMatrix::_get_row_lazy ( vector<float> & r_out , int i ) const {
      
      if ( m[i] == NULL ) {
          FSEEK(matrix_fh, (off64_t)this->data_offset + i * cols() * sizeof(value_type) , SEEK_SET);
#ifdef DEBUG
	  //          std::cerr << "lazy reading row # " << i << ", offset:" << FTELL(matrix_fh) << std::endl;
#endif
          fread(&(r_out[0]), sizeof(float) , cols() , matrix_fh);
      }
      else {
          _get_row_full( r_out, i );
      }
      //seek to the correct position in the file
  }
  void LazyMatrix::_get_row_retain ( vector<float> & r_out, int i ) {
        if ( m[i] == NULL ) {
             m[i] = new vector<float>(cols());
             FSEEK(matrix_fh, (off64_t)this->data_offset + i * cols() * sizeof(value_type) , SEEK_SET);
#ifdef DEBUG
          std::cerr << "lazy reading row # " << i << ", offset:" << FTELL(matrix_fh) << std::endl;
#endif
             std::vector<float>::iterator b = this->m[i]->begin();
             fread(&(r_out[0]), sizeof(float) , cols() , matrix_fh);
             copy(r_out.begin(), r_out.end(), m[i]->begin() );
        }
        else {
           _get_row_full( r_out, i);
        }
  }
  void LazyMatrix::get_row ( vector<float> & r_out, int i ) {
      if (sp_level == LM_FULL_DATA) {
            _get_row_full(r_out,i);
      }
      else if ( sp_level == LM_RETAIN_DATA ) {
           _get_row_retain(r_out,i);
      }
      else if ( sp_level == LM_SPARSE_DATA ) {
             _get_row_lazy(r_out,i);
      }
      else {
#ifdef CRAW_LOGGING
	      LOGH(LOG_V1,std::cerr) << "undefined sp_level" << std::endl;
#endif
      }

   }

  void LazyMatrix::_get_row_full( vector<float> & r_out , int i ) const{
	  assert(r_out.size() == cols() );
	  copy(m[i]->begin(),m[i]->end(),r_out.begin());
  }

  void LazyMatrix::get_col ( vector<float> & c_out, int j ) {
      if ( sp_level != LM_FULL_DATA ){
		  cerr << "converting matrix to full_data mode" << endl;
		  to_full_data();
	  }
    assert(c_out.size() == nrows);
    for ( int i = 0 ; i < nrows ; i++ ) {
      c_out[i] = (*this)(i,j);
    }
  }
  
  

  void LazyMatrix::remove_rows(int start_idx, int stop_idx) {
     assert( (start_idx >= 0 && start_idx <= rows() && start_idx <= stop_idx ) && 
                  ( stop_idx >= 0 && stop_idx <= rows() ) );
    data_type::iterator start = m.begin() + start_idx;
	data_type::iterator stop = m.begin() + stop_idx;
	const int n_copies = *data_copies;
	if ( n_copies > 1 ) { //TODO this seems wrong...
	for ( data_type::iterator i = start ;
	      i <  stop ; i++ ) {
	      	if ( *i != NULL ) {
	      	  delete *i ;
	      	  *i = NULL;
	      	}	
	      }
	}
	m.erase(start,stop);
	nrows = nrows - ( stop_idx - start_idx + 1 );
  };

  

  void LazyMatrix::remove_start_rows(int n ) {
  //assert( ((nrows - n) > 0 ) );
    //data_type::iterator begin = m.begin();
    data_type::iterator begin = m.begin();
	const int n_copies = *data_copies;
	if ( n_copies > 1 ) {
	for ( data_type::iterator i = m.begin() ;
	      i < begin + n ; i++ ) {
	      	if ( *i != NULL ) {
	      	  delete *i ;
	      	  *i = NULL;
	      	}	
	      }
	}
	m.erase(begin, begin + n);
		//and reduce row number counter
	nrows = nrows - n;
		
  };

  void LazyMatrix::remove_end_rows(int n) {
  //    assert( ((nrows - n) > 0)  );
	data_type::iterator end = m.end();
	const int n_copies = *data_copies;
	//TODO -- why worry about ref counting here??
	//TODO -- do I need copy-on-write ? :(
	if ( n_copies > 1 ) {
	for ( data_type::iterator i = m.end();
	    i >= end - n; i-- ) {
	    	if ( *i != NULL ) {
	    	  delete *i;
	    	  *i = NULL;
	    	}
	    }
	}
	m.erase(end - n, end);
	nrows = nrows - n;    
  };
  

  void LazyMatrix::remove_cols( int start_idx, int stop_idx ) {
     assert( (start_idx >= 0 && start_idx <= rows() && start_idx <= stop_idx ) && 
                  ( stop_idx >= 0 && stop_idx <= rows() ) );
	  for ( int i = start_idx; i < stop_idx ; i++ ) {
		  m[i]->erase(m[i]->begin() + start_idx, m[i]->begin() + stop_idx);
	  }
	  ncols = cols() - ( start_idx - stop_idx + 1);
  }


  void LazyMatrix::remove_start_cols( int n ) {
	  for ( int i = 0 ; i < rows() ; i++ ) {
		  m[i]->erase(m[i]->begin(), m[i]->begin() + n);
	  }
	  ncols = cols() - n;

  }
  void LazyMatrix::remove_stop_cols( int n ) {
	  for ( int i =0 ; i < rows() ; i++ ) {
		  m[i]->erase(m[i]->end() - n , m[i]->end());
	  }
	  ncols = cols() - n;
  }

  void LazyMatrix::add_val ( float v ) {
     to_full_data();
	 for ( uint i = 0 ; i < nrows ; i++ ) {
        vector<float> * r = get_row_p(i);
		for ( uint e = 0 ; e < r->size(); e++ ) {
           (*r)[e] += v;
		}
	 }
  }
  void LazyMatrix::mult_val ( float v ) {
     to_full_data();
	 for ( uint i = 0 ; i < nrows ; i++ ) {
        vector<float> * r = get_row_p(i);
		for ( uint e = 0 ; e < r->size(); e++ ) {
           (*r)[e] *= v;
		}
	 }
  }
  void LazyMatrix::div_val ( float v ) {
     to_full_data();
	 for ( uint i = 0 ; i < nrows ; i++ ) {
        vector<float> * r = get_row_p(i);
		for ( uint e = 0 ; e < r->size(); e++ ) {
           (*r)[e] /= v;
		}
	 }
  }
  void LazyMatrix::set_min_val( float v ) {
    for ( uint i = 0 ; i < nrows ; i++ ) {
      vector<float> * r = get_row_p(i);
      for ( uint e = 0 ; e < r->size(); e++ ) {
	if ( (*r)[e] < v  ) {
	  (*r)[e] = v;
	}
      }
    
    }
  }

  void LazyMatrix::set_max_val( float v ) {
    for ( uint i = 0 ; i < nrows ; i++ ) {
      vector<float> * r = get_row_p(i);
      for ( uint e = 0 ; e < r->size(); e++ ) {
	if ( (*r)[e] > v  ) {
	  (*r)[e] = v;
	}
      }
    
    }
  }

//interpret the shift as a modification to apply
//i.e. shift == 5 means shift all rows by +5
    
/* TODO FIX THIS TO ACTUALLY SHIFT LINEARLY, AND ADD THAT FUNCTIONALITY TO CRAW_C TO TEST FIRST */

void LazyMatrix::column_shift( int shift ) {
  vector<value_type> tmp_row(cols(),(value_type)0);
  for ( uint i = 0 ;i < rows() ; i++ ) {
    for ( uint j = 0 ; j < cols() ; j++ ) {
      uint dest_idx = j + shift;
      if ( dest_idx < 0 ) continue;
      if ( dest_idx >= cols() ) break;
      tmp_row[dest_idx] = (*this)(i,j);
    }
    set_row(tmp_row,i);

    // a little inefficient, but this ensures that we zero-out data
    // beyond the shift constraints
    tmp_row.assign(tmp_row.size(),(value_type)0);

  }
}

void LazyMatrix::linear_shift( int shift ) {
  deque< vector<value_type> * > buf(0);
  if ( shift > 0 ) {
    //fill the buffer with the first shift rows
    for ( uint i = 0; i < shift ; i++ ) {
      buf.push_back(m[i]);
      m[i] = new vector<value_type>(cols(),(value_type)0);
    }
    for ( uint i = shift; i < nrows; i++ ) {
      buf.push_back(m[i]);
      m[i] = buf[0];
      buf.pop_front();
    }
    for ( uint i = 0 ; i < buf.size() ; i++ ) {
      delete(buf[i]);
    }
  }

  else if ( shift < 0 ) {
    uint offset = abs(shift);
    for ( uint i = 0 ; i < offset ; i++ ) {
      delete(m[i]);
      m[i] = m[i+offset];
    }
    for ( uint i = offset ; i < nrows - offset; i++ ) {
      m[i] = m[i+offset];
    }
    for ( uint i = nrows - offset ; i < nrows; i++ ) {
      m[i] = new vector<value_type>(cols(),(value_type)0);
    }
  }
}  



