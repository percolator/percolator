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

#include "msmat.H"
#include "msmat_base.H"
#include "msmat_header_parser.H"
#include "crawutils.H"
#ifdef CRAW_LOGGING
#include "logging.H"
#endif
#include <stdio.h>
#include <list>
#include <map>
#include <algorithm>
#include <set>





/* msmat */

std::vector<transition_info> transition_manager::get_trans_by_precursor( MSMAT_MZ_TYPE precursor ) {
       for ( int i = 0 ; i < precursors.size() ; i++ ) 
       {
           if ( fabs(precursor - precursors[i]) < TRANSITION_MATCH_MZ_TOL ) {
               return get_trans_by_precursor(i);
           }
       }
       //if precursor not found, throw an error
       throw("did not find precursor...");
    }
std::vector<transition_info> transition_manager::get_trans_by_precursor_idx (int idx) {
        std::vector<transition_info> ti( products_by_transitions[idx].size() );
        for ( int i = 0 ; i < ti.size() ; i++ ) {
           ti[i] = transitions[ products_by_transitions[idx][i] ];
        }
        return ti;
    }
transition_manager::transition_manager ( std::vector<transition_info> transitions ) {
      this->transitions = transitions;
      typedef std::map< MSMAT_MZ_TYPE, std::vector<int> > tbp_t;
      tbp_t trans_by_precs;
      for ( int  i = 0 ; i < transitions.size() ; i++ ) {
          tbp_t::iterator k = trans_by_precs.find(transitions[i].precursor_mz);
          if ( k == trans_by_precs.end() ) {
              std::vector<int> idxs;
             idxs.push_back(i);
             trans_by_precs[transitions[i].precursor_mz] = idxs;
          }
          else {
              k->second.push_back(i);
          }
      }
    }


std::pair<MSMAT_MZ_TYPE, MSMAT_MZ_TYPE> transition_info::trans_from_str ( const std::string & s ) {
    int dash_pos = s.find('-');
    std::string prec_str = s.substr(0,dash_pos);
    std::string trans_str = s.substr(dash_pos+1, s.size() - (dash_pos+1) );
    MSMAT_MZ_TYPE prec_val = atof(prec_str.c_str());
    MSMAT_MZ_TYPE trans_val = atof(trans_str.c_str());
    return std::pair<MSMAT_MZ_TYPE, MSMAT_MZ_TYPE>(prec_val, trans_val);
}

void msmat::init_transitions_from_labels() {
   this->transitions.clear();
   if ( this->labels.size() != this->mzs.size() ) {
       throw("labels size does not match mzs size...");
   }
   this->transitions.reserve(this->labels.size());
   for ( int i = 0 ; i < labels.size() ; i++ ) {
      this->transitions.push_back( transition_info(this->labels[i],i) );
   }
}

void msmat::init_transition_manager() {
    this->t_manager = transition_manager(this->transitions);
}


msmat::msmat(FILE * msmat_f ) {

  init_from_msmat_file(msmat_f);
  calc_summary_data();
};

/* todo -- check that all the input is sorted */
void msmat::set_mzs( const vector<float> & in_mzs ) {
  mzs.clear();
  mzs.resize(in_mzs.size());
  copy(in_mzs.begin(),in_mzs.end(),mzs.begin());

}


void msmat::tabular_input  ( std::ifstream & is ) {
     this->data->to_full_data();
     is.seekg(0);
     std::vector<float> chrom_data;
     chrom_data.resize(this->get_num_rts());
     float t;
     int num_chroms = this->get_num_mzs();
     int num_rts = this->get_num_rts();
     for ( int mz_idx = 0 ; mz_idx < this->get_num_mzs() ; mz_idx++ ) {
         for ( int c_idx = 0 ; c_idx < this->get_num_rts() ; c_idx++ ) {
             is >> t;
             chrom_data.at(c_idx) = t;
         }
	 this->set_chrom( chrom_data, mz_idx );
     }
  }

void msmat::to_full_data ( bool force ) const {
  if ( loaded_full_data == false || force ) {
    
    this->data->to_full_data(force);
    if ( this->msmat_file != NULL && ( loaded_full_data == false || force ) ) {
      msmat * volatile_this = const_cast<msmat*>(this); 
      volatile_this->init_matrix_from_file( this->msmat_file, LM_FULL_DATA );
      loaded_full_data = true;
    }
  }
}

void msmat::set_min_val ( float min_I ) {
  this->data->to_full_data();
  this->data->set_min_val(min_I); // replaces anything below min_I with min_I in the LazyMatrix for this msmat
  this->summary_data_ok = false;
}

void msmat::set_max_val ( float max_I ) {
  this->data->to_full_data();
  this->data->set_max_val(max_I);
  this->summary_data_ok = false;
}

void msmat::set_mzs ( float * in_mzs, int num_mzs ) {
  /* less efficient since I cannot use copy, I believe */

  mzs.clear();
  mzs.resize(num_mzs);
  copy(in_mzs,in_mzs+num_mzs,mzs.begin());

}
void msmat::set_rts ( const vector<float> & in_rts ) {
  rts.clear();
  rts.resize(in_rts.size());
  copy(in_rts.begin(),in_rts.end(),rts.begin());

}
void msmat::set_rts ( float * in_rts, int num_rts ) {
  rts.clear();
  rts.resize(num_rts);
  copy(in_rts, in_rts+num_rts, rts.begin());
}

void msmat::set_ort_wrt_map( flt_map & ow_map ) {
  ort_wrt_map = ow_map;
}


void msmat::set_ort_wrt_map( const vector<float> & ort,
                                  const vector<float> & wrt ) {
                                   assert(ort.size() == wrt.size());
                                   /* ugly */
                                   for ( size_t i = 0 ; i < ort.size() ; i++ ) {
                                     ort_wrt_map[ort[i]] = wrt[i];
                                   }
}


void msmat::zero_chrom ( float mz, float start_rt , float stop_rt ) {
  int chrom_idx = mz_to_mzidx(mz);
  zero_chrom(chrom_idx, start_rt, stop_rt);
  this->summary_data_ok = false;
}

void msmat::zero_chrom( int chrom_idx, int start_scan_idx, int stop_scan_idx ) {
  
  vector<float> tmpchr(get_num_rts(),0.0f);
  get_chrom(chrom_idx,tmpchr);
  for ( int i = start_scan_idx ; i <= stop_scan_idx ; i++ ) {
    tmpchr[i] = 0.0f;
  }
  set_chrom(tmpchr, chrom_idx); 
  this->summary_data_ok = false;   
}



void msmat::zero_chrom ( int chrom_idx , float start_rt , float stop_rt ) {
  int start_rt_idx = 0;
  int stop_rt_idx = get_num_rts() - 1; 
  if ( start_rt == -1 && stop_rt == -1 ) {
    
  }
  else {
    if ( start_rt != -1.0f ) {
      start_rt_idx = get_scan_idx_by_rt(start_rt);
    }
    if ( stop_rt != -1.0f ) {
      stop_rt_idx = get_scan_idx_by_rt(stop_rt);
    }
  }
  zero_chrom(chrom_idx, start_rt_idx, stop_rt_idx);
  this->summary_data_ok = false;
}



void msmat::retain_dist_set( int chrom_idx, std::vector<int> & v ) {
    std::set<int> m(v.begin(), v.end());
    for ( int i = 0 ; i < this->get_num_rts() ; i++ ) {
        if ( m.count(i) == 0 )  {
           this->zero_pixel( chrom_idx, i );
        }
    }
    this->summary_data_ok = false;
}

void msmat::retain_chrom_set ( int chrom_idx , std::vector< std::pair< int,  int > > & v ) {
    int set_idx = 0;

    std::vector< std::pair<int,int> > retain_range = crawutils::unique_ranges( v , true);
    std::vector< std::pair<int,int> > delete_range;
    if ( retain_range[0].first > 0 ) {
        delete_range.push_back(std::pair<int,int>(0,retain_range[0].first-1) );
    }

    int del_start = retain_range[0].second + 1; 
    int del_stop;
    for ( int i = 1 ; i < retain_range.size() ; i++ ) {
        del_stop = retain_range[i].first;
        delete_range.push_back(std::pair<int,int>(del_start, del_stop) );
        del_start = retain_range[i].second + 1;
    }
    del_stop = get_num_rts() - 1;
    delete_range.push_back(std::pair<int,int>(del_start,del_stop) );

    for ( int i = 0 ; i < delete_range.size() ; i++ ) {
       this->zero_chrom(chrom_idx,delete_range[i].first, delete_range[i].second);
    }
    this->summary_data_ok = false;
}

void msmat::retain_chrom ( int chrom_idx, int start_scan_idx, int stop_scan_idx ) {
  vector<float> tmpchr(get_num_rts(),0.0f);
  for ( int i = 0 ; i < start_scan_idx ; i++ ) {
    tmpchr[i] = 0.0f;
  }
  for ( int i = stop_scan_idx + 1 ; i < get_num_rts() ; i++ ) {
    tmpchr[i] = 0.0f;
  }
  set_chrom(tmpchr,chrom_idx);
  this->summary_data_ok = false;
}

/* TODO -- check arithmetic with f_off type */

int msmat::check_size_from_position() const {
  off64_t cur_pos, end_pos, zero_pos;
  int ret_val;
  zero_pos = 0;
  cur_pos = FTELL(msmat_file);
  FSEEK(msmat_file, zero_pos, SEEK_END);
  end_pos = FTELL(msmat_file);

  if ( ( end_pos - cur_pos ) == get_num_mzs() * get_num_rts() * sizeof(float) ) {
    ret_val = 1;
  }
  else {
    ret_val = 0;
  }
  FSEEK(msmat_file,cur_pos,SEEK_SET);
  return ret_val;
}

void msmat::init_from_msmat_file( FILE * msmat_f ) {	
  init_from_msmat_file(msmat_f, LM_FULL_DATA);
  if ( msmat_file != NULL ) {
    
  }
}


  void msmat::init_self_from_header() {
    //FSEEK(msmat_file,0,SEEK_SET);
    hp->read_header(this,msmat_file, false, NULL);
    this->data_start_offset = hp->get_data_start();
    if ( ! this->check_size_from_position() ) {
      cerr << "msmat file is not positioned properly" << endl;
      exit(-1);
      //TODO -- raise exception
    }
  }

void msmat::init_from_msmat_file (FILE * msmat_f, sparse_level l) {

  init_msmat();	
  msmat_file = msmat_f;
  sp_level = l;
  init_self_from_header();
  init_data();
  init_matrix_from_file(msmat_f,sp_level);
}


int msmat::get_scan_idx_by_rt ( float rt ) const {
  //TODO -- check that rts is sorted when it is assigned
  vector<float>::const_iterator nr_val = find_nearest( rts, rt);
  /* now we see if we can determine the index of this iterator */
  int p_diff = nr_val - rts.begin();
  assert(rts[p_diff] == *nr_val);
  return p_diff;
}

int msmat::get_chrom_idx_by_mz(float mz) const {
  vector<float>::const_iterator nr_val = find_nearest( mzs, mz );
  /* now we see if we can determine the index of this iterator */
  int p_diff = nr_val - mzs.begin();
  assert(mzs[p_diff] == *nr_val);
  return p_diff;
}


//erase goes from [ start, stop )
void msmat::trim_mzs_list_bound ( int start_idx, int stop_idx ) {
  if ( labels.size() == mzs.size() ) {
    labels.erase(labels.begin() + start_idx, labels.begin() + stop_idx);
  }
  mzs.erase(mzs.begin() + start_idx, mzs.begin() + stop_idx);
}

void msmat::trim_mzs_data_bound ( int start_idx, int stop_idx ) {
  _trim_start_mzs(start_idx);
  _trim_stop_mzs(stop_idx);
}

void msmat::trim_mzs_start( float keep_start_mz ) {
  mzs_type::iterator lk = lower_bound(mzs.begin(), mzs.end(), keep_start_mz);
  int lk_nmzs = lk - mzs.begin();  
  _trim_start_mzs(lk_nmzs);
  _trim_start_mzs_list(lk_nmzs);
}
void msmat::trim_mzs_stop( float keep_stop_mz ) {
  mzs_type::iterator uk = lower_bound(mzs.begin(), mzs.end(), keep_stop_mz);
  int uk_nmzs = uk - mzs.begin();
  _trim_stop_mzs(mzs.size() - 1 - uk_nmzs);
  _trim_stop_mzs_list(mzs.size() - 1 - uk_nmzs);

}


void msmat::trim_mzs( float keep_start_mz , float keep_stop_mz ) {


}

  void msmat::init_msmat()
{
  bin_size = -1.0;
  data_start_offset = -1;
  //TODO make sure this is not already checked
  data = NULL;
  msmat_file = NULL;
  summary_data_ok = false;
  //hp.attach_msmat(*this);
  sp_level = LM_FULL_DATA;
  hp = new msmat_header_parser();
  loaded_full_data = false;
  f_info = NULL;
  //init_data();
}

/* finds the appropriate bin by converting mz to mzidx by searching through the 
   actual mzs */


 int msmat::mz_to_mzidx( float mz ) const {
     if ( mz < mzs[0] ) {
         throw(MZRangeError(MZRangeError::MZLowError) );
     }
     if ( mz > mzs[mzs.size()-1] ) {
         throw(MZRangeError(MZRangeError::MZHighError) );
     }

     std::vector<float>::const_iterator lb = std::upper_bound( this->mzs.begin() , this->mzs.end() , mz );
     int iter_delt = lb - this->mzs.begin() - 1;
     return iter_delt;
  }

int msmat::rt_to_rtidx ( float rt ) const {
    if ( rt < rts[0] ) {
         throw(RTRangeError(RTRangeError::RTLowError) );
     }
     if ( rt > rts[rts.size()-1] ) {
         throw(RTRangeError(RTRangeError::RTHighError) );
     }

     std::vector<float>::const_iterator lb = std::upper_bound( this->rts.begin() , this->rts.end() , rt );
     int iter_delt = lb - this->rts.begin() - 1;
     if ( iter_delt == rts.size() - 1 ) 
     {
         if ( fabs(rts[iter_delt] - rt) > 0.0001 ) {
            throw(RTRangeError(RTRangeError::RTHighError));
         }
         else {
           return iter_delt;
         }
     }
     float lh_rt, rh_rt;
     lh_rt = rts[iter_delt];
     rh_rt = rts[iter_delt+1];
     if ( (rt - lh_rt) < (rh_rt - rt ) ) {
       return iter_delt;
     }
     else {
       return iter_delt + 1;
     }
}


/*
  int msmat::mz_to_mzidx( float mz ) const {
     float mz_low, mz_high;
     mz_low = mzs[0];
     mz_high = mzs[mzs.size() - 1];
     if ( mz < mz_low )
        return -1;
     if ( mz > mz_high )
       return mzs.size() - 1;
     else {
        int idx = (int)floor(( mz - mz_low ) / get_bin_size());
        return idx;
     }
  }
*/

/*
  int msmat::rt_to_rtidx( float rt ) const {
     float rt_low, rt_high;
     rt_low = rts[0];
     rt_high = rts[rts.size() - 1];
     if ( rt < rt_low ) 
         throw(RTRangeError(RTRangeError::RTLowError));
     if ( rt > rt_high ) 
         throw(RTRangeError(RTRangeError::RTHighError));
     else {
        int idx = (int)floor(( rt - rt_low ) / get_bin_size());
        return idx;
     }
  }
*/

/*

void msmat::trim_mzs ( mz_type mz_min, mz_type mz_max ) {
  mzs_type::iterator ub_start = upper_bound(mzs.begin(),mzs.end(),mz_min);
  mzs_type::iterator lb_stop  = lower_bound(mzs.begin(),mzs.end(),mz_max);
  
  int n_start_mzs = ub_start - mzs.begin();
  int n_stop_mzs  = mzs.end() - lb_stop;

  if ( n_start_mzs > 0 ) {
    trim_mz_start(n_start_mzs);
  }
  if ( n_stop_mzs > 0 ) {
    trim_mz_stop(n_stop_mzs);
  }

}

*/

void msmat::trim_start (int nscans) {
  _trim_start(nscans);
}
void msmat::trim_start( rt_type rt ) {
  rts_type::iterator ub = upper_bound(rts.begin(), rts.end(), rt );
  int nscans = ub - rts.begin();
  if ( nscans == 0 ) {
    cerr << "Complete trimming of run at trim_start " << rt << endl;
  }
  trim_start(nscans);
}

void msmat_scans::_trim_start( int nscans ) {
  /* first deal with non-scan data */
  _trim_start_rts(nscans);
  /* now with the LazyMatrix * data */
  data->remove_start_rows(nscans);
  this->summary_data_ok = false;
}

void msmat::trim_stop( int nscans ) {
  _trim_stop(nscans);
}
void msmat::trim_stop( rt_type rt ) {
  rts_type::iterator lb = lower_bound(rts.begin(), rts.end(), rt);
  int nscans = rts.end() - lb;
  if ( nscans == 0 )  {
    cerr << "Complete trimming of run at trim_stop " << rt << endl;
  }

  trim_stop(nscans);
}

void msmat_scans::_trim_stop( int nscans ) {
  /* */
  _trim_stop_rts(nscans);
  data->remove_end_rows(nscans);
  this->summary_data_ok = false;
}

void msmat_scans::_trim_start_mzs( int nscans ) {
  data->remove_start_cols(nscans);
  this->summary_data_ok = false;
}
void msmat_scans::_trim_stop_mzs( int nscans ) {
  data->remove_stop_cols(nscans);
  this->summary_data_ok = false;
}


void msmat_chroms::_trim_start( int nscans ) {
  _trim_start_rts(nscans);
  /* we keep all rows in LazyMatrix -- what we do is erase from
  * each LazyMatrix Vector */
  data->remove_start_cols(nscans);
}
void msmat_chroms::_trim_stop ( int nscans ) {
  _trim_stop_rts(nscans);
  data->remove_stop_cols(nscans);

}
void msmat_chroms::_trim_start_mzs( int nscans ) {
  data->remove_start_rows(nscans);
}
void msmat_chroms::_trim_stop_mzs( int nscans ) {
  data->remove_end_rows(nscans);
}


void msmat::_trim_start_rts( int nscans ) {
  vector<float> new_rts(get_num_rts() - nscans);
  copy(rts.begin() + nscans, rts.end(), new_rts.begin());
  rts.swap(new_rts);
}

void msmat::_trim_stop_rts ( int nscans ) {
  vector<float> new_rts(get_num_rts() - nscans);
  copy(rts.begin(), rts.end() - nscans, new_rts.begin());
  rts.swap(new_rts);
}

void msmat::_trim_start_mzs_list ( int nscans ) {
  vector<float> new_mzs(get_num_mzs() - nscans);
  copy(mzs.begin() + nscans, mzs.end(), new_mzs.begin());
  mzs.swap(new_mzs);
}

void msmat::_trim_stop_mzs_list ( int nscans ) {
  vector<float> new_mzs(get_num_mzs() - nscans);
  copy(mzs.begin(), mzs.end() - nscans, new_mzs.begin());
  mzs.swap(new_mzs);
}


void msmat::trim_start_stop( int start_scans, int stop_scans) {
  trim_start(start_scans);
  trim_stop(stop_scans);
}
void msmat::trim_start_stop( rt_type start_rt, rt_type stop_rt) {
  trim_start(start_rt);
  trim_stop(stop_rt);
}





/* protected data init functions */

void msmat::init_data() {
  matrix_dim d = get_data_dim();
  data = new LazyMatrix(d.row_num,d.row_size,get_sparse_level(), this->msmat_file, data_start_offset );  
  cerr << "allocated matrix:" << d.row_num << "," << d.row_size << endl;
  cerr << "matrix size:" << data->m.size() << endl;
}

void msmat::copy_from_msmat_nodata( const msmat & msmat ) {

  bin_size = msmat.bin_size;
  mzs = msmat.mzs;
  rts = msmat.rts;
  msmat_file = NULL;
  //TODO does this need to be copied deeply?
  ort_wrt_map = msmat.ort_wrt_map;
  sp_level = msmat.sp_level;
  labels = msmat.labels;
  loaded_full_data = false;
  this->original_file_name = msmat.original_file_name;
  this->original_prefix    = msmat.original_prefix;

}



/* copy_from_msmat_mzs --
   Copy data from t to this msmat based on a range of mz indexes passed in
*/
void msmat::copy_from_msmat_mzs (msmat & t , const std::vector<int> mzidxs) {
  //assert(mzs.size() == t.mzs.size());
  this->copy_from_msmat_nodata( t );
  this->retain_mzs_nodata( mzidxs );
  //check that we have only the correct number of mzs data now
  assert(mzs.size() == mzidxs.size()); 
  //initialize our data to be the correct size
  this->init_data();
  vector<float> v(this->get_num_mzs());
  for ( uint i = 0 ; i < mzidxs.size() ; i++ ) {
    t.get_chrom(mzidxs[i], v);
    this->set_chrom(v, i);
  }
  this->summary_data_ok = false;
}


void msmat::retain_mzs_nodata( const vector<int> mzidxs ) {
  //new_mzs is made from our current set of mzs
  vector<float> new_mzs(mzidxs.size());
  for ( uint i = 0; i < mzidxs.size() ; i++ ) {
    new_mzs[i] = mzs[mzidxs[i]];
  }
  mzs.swap(new_mzs);
  this->summary_data_ok = false;
}


void msmat::copy_from_msmat( const msmat & msmat) {
  copy_from_msmat_nodata(msmat);
  data = msmat.data;
  summary_data_ok = msmat.summary_data_ok;
  if ( ! summary_data_ok && ( data != NULL ) ) {
    calc_summary_data();
  }
}

void msmat::copy_from_msmat_deep( const msmat & msmat ) {
  //this->to_full_data();
  //data->to_full_data();
  copy_from_msmat_nodata(msmat);
  data = new LazyMatrix(*(msmat.data),false);
  summary_data_ok = msmat.summary_data_ok;
  if ( ! summary_data_ok && ( data != NULL ) ) {
    calc_summary_data();
  }
}

void msmat::clear_data() {
  delete(data);
}

void msmat::calc_bp_chrom() {
  int scan_idx, num_scans, scan_length;
  scan_idx = 0;
  num_scans = get_num_rts();
  scan_length = get_num_mzs();
  bp_chrom.resize(num_scans);
  vector<float> scan(get_num_mzs());
  for ( ; scan_idx < num_scans; scan_idx++ ) {

    get_scan(scan_idx,scan);
    float scan_max = *(max_element(scan.begin(), scan.end()));
    bp_chrom[scan_idx] = scan_max;
  }
}

void msmat::calc_tic_chrom() {
  int scan_idx, num_scans, scan_length;
  scan_idx = 0;
  num_scans = get_num_rts();
  scan_length = get_num_mzs();
  tic_chrom.resize(num_scans);
  vector<float> scan(get_num_mzs());
  for ( ; scan_idx < num_scans; scan_idx++ ) {
    //get the ptr to the first element

    get_scan(scan_idx,scan);

    float * scan_ptr = &(scan[0]);
    float scan_sum = sum_list<float>(scan_ptr, scan_length);
    tic_chrom[scan_idx] = scan_sum;
  }
}



/*double msmat::get_mode_rt_interval() {
size_t i;
const float granularity = 1e-3;
const float lb = 0.0;
const float rb = 1.0;
const float nbins_f = (rb - lb) / granularity;
int nbins = (int)round(nbins_f);
vector<int> bins(nbins);
for ( i = 0; i < bins.size(); i++ ) {
bins[i] = 0;
}

for ( i = 1; i < rts.size() ; i++ ) {
float interval = rts[i] - rts[i-1];
if ( interval > rb ) {
LOGH(LOG_V1,std::cerr) << "scan interval of > 1 min. code needs updating" << std::endl; exit(1);
}
int bin_idx = (int)(interval / granularity);
bins[bin_idx]++;
}

vector<int>::iterator me = max_element(bins.begin(), bins.end());
int me_idx = me - bins.begin();
float interval = me_idx * granularity;
return (interval); 
}*/


void msmat::get_mode_rt_intervals(float granularity, float lb, float rb, std::vector<uint> & bins ) const {

  for ( uint i = 1; i < rts.size() ; i++ ) {
    float interval = rts[i] - rts[i-1];
    int bin_idx;
    if ( interval > rb ) {
      bins[bins.size() - 1]++;
    }
    else {
      bin_idx = (uint)(interval / granularity );
      bins[bin_idx]++;
    }
  }

}
double msmat::get_mode_rt_interval() const {

  const float granularity = 1e-3f;
  const float lb = 0.0;
  const float rb = 1.0;
  const float nbins_f = (rb - lb) / granularity;
  int nbins = (int)round(nbins_f);
  vector<uint> bins(nbins,0);

  get_mode_rt_intervals(granularity, lb, rb, bins);

  vector<uint>::iterator me = max_element(bins.begin(), bins.end());
  int me_idx = me - bins.begin();
  float interval = me_idx * granularity;
  return (interval); 
}

double msmat::get_mode_rt_interval(double roundto) const {
  double interval = get_mode_rt_interval();
  double rnd_interval = round_to_dbl(interval, roundto);
  return rnd_interval;
}

void msmat::pad_start_stop( int start_mapped_coord, 
                                int stop_mapped_coord ) {
                                  vector<float>null_vect(get_num_mzs(),0.0f);
                                  for ( int i = 0 ; i < start_mapped_coord ; i++ ) {
                                    vector<float> nv(null_vect);
                                    set_scan(nv,i);
                                  }
                                  for ( uint i = stop_mapped_coord+1 ; i < get_num_rts()  ; i++  ) {
                                    vector<float> nv(null_vect);
                                    set_scan(nv,i);
                                  }
  this->summary_data_ok = false;
}
/*
def pad_start_stop(self,start_mapped_coord,stop_mapped_coord,template_rts,null_scan) :
'''pad self with zero intensity scans upto and after the mapped coordinates
mentioned above, including correct retention times'''
for i in range(start_mapped_coord) :
null = null_scan.copy()
null.null_intensities()
null.scanid = i + 1
null.rtime = template_rts[i]
self.scans[i] = null
self.rts[i] = null.rtime
print >> sys.stderr , "adding a zeroed scan at " , i

for i in range(stop_mapped_coord+1 ,len(self.scans) - 1) :
null = null_scan.copy()
null.null_intensities()
null.scanid = i + 1
null.rtime = template_rts[i]
self.scans[i] = null
self.rts[i] = null.rtime
print >> sys.stderr , "adding a zeroed scan at " , i
*/
















void msmat::interpolate_new_scan( const vector<float> & lh_scan,
                                      const vector<float> & rh_scan, vector<float> & interpolated_scan, float xdelt) const {

                                        assert(lh_scan.size() == rh_scan.size());
                                        assert(interpolated_scan.size() == lh_scan.size());
                                        for ( uint i = 0; i < lh_scan.size() ; i++ ) {
                                          interpolated_scan[i] = lh_scan[i] + ( rh_scan[i] - lh_scan[i]) * xdelt;
                                        }
}






void msmat::sample_values( std::vector<float> & vals , const int num_vals, 
				float min_threshold, float max_threshold ) {
   throw("sample_values not implemented");



}


void msmat::warp_average_new_output (const vector<float> & t_coords,
                                          const vector<float> & a_coords, 
                                          const msmat & templ, 
                                          const msmat & align,
                                          int linear_scans_offset ) {
                                              this->summary_data_ok = false;
                                  int start_temp_coord = (int)t_coords[0];
                                  int stop_temp_coord = (int)t_coords[t_coords.size()-1];
                                    //pad the start and stop of the alignment run with zeroed data when there is no
                                    //mapped template scan
                                  pad_start_stop( start_temp_coord, stop_temp_coord);

                                    vector<float> interpolated_scan(mzs.size());
                                    vector<float> lh_scan(mzs.size());
                                    vector<float> rh_scan(mzs.size());

                                    for ( uint i = 0 ; i < t_coords.size() ; i++ ) {
                                      int t_coord = (int)t_coords[i];
				      /* TODO --When applying a linear shift, some of the logic for the retention times at the
					  beginning and end of the path need to be adjusted. I.E. as it stands,
					  the adjustment using linear_scans_offset could result in align scan 0 or align_scan[last]
					  being mapped multiple times -- should actually refuse to map those where the alignment scan
					  is invalid -- but what effect would that have in Diffs.py, if an ort cannot be found is just
	                                  the mapping for that run lost?												   
													   
					  */
                                      float al_coord = (float)std::min( std::max((int)(a_coords[i] - linear_scans_offset),0),
                                                                          (int)align.get_num_rts() - 1);
                                                         


                                      ort_wrt_map.insert(pair<float,float>( align.rt_from_scannum(al_coord),
									    templ.rt_from_scannum(t_coords[i]) ) );


                                      /* interpolate the appropriate align scan */
                                      if ( a_coords[i] >= align.get_num_rts() - 1 ) {
                                        align.get_scan(align.get_num_rts() - 1, interpolated_scan);
                                      }
                                      else {
                                        int lh_idx = (int)a_coords[i];
                                        int rh_idx = lh_idx + 1;
                                        float xdelt = a_coords[i] - lh_idx;
                                        //derive new scan
                                        align.get_scan(lh_idx, lh_scan);
                                        align.get_scan(rh_idx, rh_scan);
                                        interpolate_new_scan(lh_scan,
                                          rh_scan,
                                          interpolated_scan,
                                          xdelt);
                                      }
                                      set_scan(interpolated_scan,(int)round(t_coord)); //HACK


                                      //put info into ort_wrt map
                                      //add info scans
                                      //add info rts

                                    }

}



void msmat::warp_to_template( const vector<float> & t_coords,
                                  const vector<float> & a_coords, 
                                  const msmat & templ, 
                                  const msmat & align,
                                  int linear_scans_offset ) {
                                      this->summary_data_ok = false;
                                    //linear_scans_offset is the # of scans by which the
                                    //alignment run was mapped
                                    int start_temp_coord = (int)t_coords[0];
                                    int stop_temp_coord = (int)t_coords[t_coords.size()-1];
                                    //pad the start and stop of the alignment run with zeroed data when there is no
                                    //mapped template scan
                                    pad_start_stop( start_temp_coord, stop_temp_coord);

                                    vector<float> interpolated_scan(mzs.size());
                                    vector<float> lh_scan(mzs.size());
                                    vector<float> rh_scan(mzs.size());

                                    for ( uint i = 0 ; i < t_coords.size() ; i++ ) {
                                      int t_coord = (int)round(t_coords[i]);
				      /* TODO --When applying a linear shift, some of the logic for the retention times at the
					  beginning and end of the path need to be adjusted. I.E. as it stands,
					  the adjustment using linear_scans_offset could result in align scan 0 or align_scan[last]
					  being mapped multiple times -- should actually refuse to map those where the alignment scan
					  is invalid -- but what effect would that have in Diffs.py, if an ort cannot be found is just
	                                  the mapping for that run lost?												   
													   
					  */


                                      float al_map_coord = std::min( std::max((a_coords[i] - float(linear_scans_offset)), 0.0f),
								     (float)align.get_num_rts() - 1.0f );
                                                         
				      float al_coord = a_coords[i];


				      float al_rt, tm_rt;
				      al_rt = align.rt_from_scannum(al_map_coord);
				      tm_rt = templ.rt_from_scannum(t_coords[i]);
				      
				      /* yes, the alignment time (what we copy from) is original, and the template time
					 (what we copy to) is warped */
                                      ort_wrt_map.insert(pair<float,float>( al_rt, tm_rt ));
				      std::cout << i << " " << t_coords[i] << " " << a_coords[i] << " " << al_coord << " " << tm_rt << " " << al_rt <<  std::endl;				      
				      


                                      /* interpolate the appropriate align scan */
                                      if ( a_coords[i] >= align.get_num_rts() - 1 ) {
                                        align.get_scan_interp((float)(align.get_num_rts() - 1), interpolated_scan);
                                      }
                                      else {
					//align.get_scan_interp((al_coord),interpolated_scan);
                                        int lh_idx = (int)a_coords[i];
                                        int rh_idx = lh_idx + 1;
                                        float xdelt = a_coords[i] - lh_idx;
                                        //derive new scan
                                        align.get_scan(lh_idx, lh_scan);
                                        align.get_scan(rh_idx, rh_scan);
                                        interpolate_new_scan(lh_scan,
                                          rh_scan,
                                          interpolated_scan,
                                          xdelt);
                                      }
                                      set_scan(interpolated_scan,(int)round(t_coord)); //HACK


                                      //put info into ort_wrt map
                                      //add info scans
                                      //add info rts

                                    }

}


void msmat_scans::linear_shift( int shift ) {
  data->to_full_data();
  data->linear_shift( shift );
  summary_data_ok = false;
}

void msmat_chroms::linear_shift( int shift ) {
  data->to_full_data();
  data->column_shift(shift);
  summary_data_ok = false;
}





void msmat::time_interpolate(double interval, msmat & new_matrix) {
  /* 1. calculate size of new matrix in scans	
  * 2. allocate matrix, and new rts
  * 3. populate new matrix by interpolating old one's data
  */

  /* set the new start and stop to be within the current bounds */


  /* guarantee that the new rt boundaries are within the current time boundaries,
  and will be a consistent scale from run to run using the same interval */
  double new_start, new_end;
  double t;
  this->summary_data_ok = false;
  for (  t = 0.0 ; t < rts[0] ; t += interval ) {
  }
  new_start = t;
  for ( t = new_start ; t < rts[rts.size()-1] ; t+= interval) {
  }
  new_end = t-interval;



  /*float new_start = round_to_dbl(rts[0], interval);
  if ( new_start < rts[0] ) {
  new_start = round_to_dbl(rts[1], interval);
  }

  float new_end = round_to_dbl(rts[rts.size() - 1], interval);
  if ( new_end > rts[rts.size() - 1]) {
  new_end = round_to_dbl(rts[rts.size() -2 ], interval);
  }
  */

  double new_len = new_end - new_start;
  int num_new_scans = (int)round(new_len / interval);

  cerr << "new_start, new_end, new_len, num_new_scans:" << new_start << " " << new_end << " " << 
    new_len << " " << num_new_scans << endl;
  // create new retention times
  vector<float> new_rts(num_new_scans);

  new_rts[0] = new_start;
  for ( uint i = 1 ; i < new_rts.size() ; i++ ) {
    new_rts[i] = new_rts[i-1] + interval;
  }

  new_matrix.set_rts(new_rts);
  new_matrix.init_data();

  //copy the attributes of new_matrix from this function
  //TODO MONDAY -- DEAL WITH THE EFFECTS ON NEW MATRIX, AND ON TRIMMING..

  cerr << "from an old start, stop of " << rts[0] << rts[rts.size() -1] << " using an interval of ";
  cerr << interval << " new start, stop are" << new_start << "," << new_end << endl;
  cerr << "difference in final length" << fabs(new_rts[new_rts.size() - 1] - new_end) << endl;

  matrix_dim dim = get_data_dim();

  uint old_rts_idx = 0;
  uint new_rts_idx = 0;
  uint jmp_idx = 0;
  uint last_lower_old = 0;
  while (true) {
    /* we want the flanking elements from rts for each new member of new_rts */
    if ( new_rts_idx >= new_rts.size() ) {
      break;
    }
    float new_rt = new_rts[new_rts_idx];
    if ( new_rt < rts[old_rts_idx]) {
      cerr << "new_rt, old_rts_idx, rts[old_rts_idx], rts.size()" << new_rt << " " << old_rts_idx;
      cerr << " " << rts[old_rts_idx] << " " << rts.size() << endl;
    }
    assert(new_rt >= rts[old_rts_idx]);
    while( rts[old_rts_idx+jmp_idx] < new_rt) {	
      if ( (old_rts_idx + jmp_idx) < (rts.size() - 1) ) 
        jmp_idx++;
      else 
        break;
    }
    if ( jmp_idx > 0 ) {
      assert(rts[old_rts_idx+jmp_idx] >= new_rt && 
        rts[old_rts_idx+jmp_idx-1] < new_rt);
    }
    else {
      assert(rts[old_rts_idx+jmp_idx] >= new_rt);
    }
    if ( new_rt < rts[old_rts_idx+jmp_idx]) {
      last_lower_old = old_rts_idx+jmp_idx-1;
    }
    else if ( new_rt == rts[old_rts_idx+jmp_idx]) {
      last_lower_old = old_rts_idx + jmp_idx;
    }
    else {
      assert(rts[old_rts_idx+jmp_idx] >= new_rt );
    }
    assert(new_rt >= rts[last_lower_old] && new_rt <= rts[last_lower_old+1]);
    vector<float> new_scan = interpolate_scan_by_rt(new_rt, last_lower_old,last_lower_old+1);
    new_matrix.set_scan(new_scan,new_rts_idx);
    old_rts_idx = last_lower_old;
    jmp_idx = 0;
    new_rts_idx++;
  }
  new_matrix.calc_bp_chrom();
  new_matrix.calc_tic_chrom();
}; 

vector<float> msmat::interpolate_scan_by_idx ( float new_idx ){
  int lh_idx, rh_idx;
  this->summary_data_ok = false;
  lh_idx = (int)floor(new_idx);
  rh_idx = (int)ceil(new_idx);
  float dx = new_idx - lh_idx;
  vector<float> lh_scan(get_num_mzs()), rh_scan(get_num_mzs());
  get_scan(lh_idx,lh_scan);
  get_scan(rh_idx,rh_scan);
  vector<float> new_scan(lh_scan.size());
  for ( uint i = 0 ; i < lh_scan.size() ; i++ ) {
    new_scan[i] = dx * ( rh_scan[i] - lh_scan[i]) + lh_scan[i];
  }
  return new_scan;
}


//void msmat::position_fh_for_dataread 

void msmat::init_matrix_from_file ( FILE * msmat_fh, sparse_level l) {

  //ignore sparse level for now..

  FSEEK(msmat_fh, data_start_offset, SEEK_SET);
  if ( ! this->check_size_from_position() ) {
    cerr << "msmat file is at incorrect position" << FTELL(msmat_fh) << endl;
    exit(-1);
  }

  

  if ( l == LM_FULL_DATA ) {

    for ( int i = 0 ; i < data->rows(); i++ ) {
#ifdef DEBUG 
      //std::cerr << "ms1mat_fh_pos , row# " << i << ":" << FTELL(msmat_fh) << std::endl;
#endif
      //cerr << data->at(i)->end() - data->at(i)->begin() << endl;
      vector<float> * v = data->at(i);
      float * beg = &((*v)[0]);
      if ( fread( beg, sizeof(float), v->size(), msmat_fh) != 
        data->at(i)->size() ) {
	  #ifdef CRAW_LOGGING
          LOGH(LOG_V1,std::cerr) << "incorrect read from msmat file\n"; 
	  #endif
	  exit(1);
      }
    }
    
    loaded_full_data = true;
    
  }
  else {
    //lazy version, do not read yet
  }
  
  
}

void msmat::write_matrix_to_file ( FILE * msmat_fh ) {
  if ( sp_level != LM_FULL_DATA ) {
    #ifdef CRAW_LOGGING
    LOGH(LOG_V2,std::cerr) << "matrix was not full for writing " << std::endl;
    LOGH(LOG_V1,std::cerr) <<  "matrix needs to be full to write to a file" << std::endl ; exit(1);
    #endif
  }
  for ( int i = 0 ; i < data->rows(); i++ ) {
    vector<float> * v = data->at(i);
    float * beg = &((*v)[0]);
    if ( fwrite(beg, sizeof(float), v->size(), msmat_fh) != 
      v->size() ) {
        #ifdef CRAW_LOGGING
        LOGH(LOG_V2,std::cerr) << "incorrect msmat write " << std::endl;
        LOGH(LOG_V2,std::cerr) << "incorrect write from msmat file\n"; exit(1);
	#endif
    }
  }
}

  void msmat::calc_summary_data() {
      if ( ! summary_data_ok ) {
          calc_bp_chrom();
          calc_tic_chrom();
          summary_data_ok = true;
      }
  }

  vector<float> msmat::get_rts() const {
    vector<float> tmp_rts(rts.size());
    copy(rts.begin(), rts.end(), tmp_rts.begin());
    return tmp_rts;
  }
  vector<float> msmat::get_mzs() const {
    vector<float> tmp_mzs(mzs.size());
    copy(mzs.begin(), mzs.end(), tmp_mzs.begin());
    return tmp_mzs;
  }


msmat * msmat::copy_no_data() { 
  msmat_scans * m = new msmat_scans();
  m->copy_from_msmat_nodata(*this);
  return (msmat*) m;
}

msmat * msmat_chroms::copy_no_data() {
    msmat_chroms * m = new msmat_chroms();
    m->copy_from_msmat_nodata(*this);
    return m;
}

void msmat::write_to_file( FILE * msmat_fh ) {
  if ( ! this->summary_data_ok )  {
    calc_summary_data();
  }
  hp->write_header(this,msmat_fh);				
  #ifdef CRAW_LOGGING
  LOGH(LOG_V2,std::cerr) << "wrote header " << std::endl;
  #endif
  write_matrix_to_file(msmat_fh);
  #ifdef CRAW_LOGGING
  LOGH(LOG_V2,std::cerr) << "wrote matrix " << std::endl;
  #endif
}

double msmat::total_current ( float start_rt, float stop_rt, float threshold ) const {
  uint start_scan_idx, stop_scan_idx;
  if ( start_rt > 0 || stop_rt > 0 ) {
    start_scan_idx = this->get_scan_idx_by_rt( start_rt );
    stop_scan_idx = this->get_scan_idx_by_rt( stop_rt );
    return total_current(start_scan_idx,stop_scan_idx,threshold);
  }
  else {
    return total_current(threshold);
  }
}

double msmat::total_current (  uint start_scan , uint last_scan,  float threshold ) const {
  vector<float> t(get_num_cols(),0.0f);
  vector<float> s(get_num_cols(),0.0f);
  for ( uint i = start_scan ; i <= last_scan ; i++ ) {
    data->get_row(s,i);
    for ( uint i = 0 ; i < s.size() ; i++ ) {
      t[i] += s[i];
    }
  }
  double total = 0.0;
  for ( uint i = 0 ; i < t.size() ; i++ ) {  
    total += t[i];
  }
  return total;
}

double msmat::get_avg() const {
  return get_avg(0);
}

double msmat::get_median() const {
  vector<float> data_sorted(this->get_num_rows() * this->get_num_cols());
  int cnt = 0;
  for ( int i = 0 ; i < this->get_num_rows() ; i++ ) {
    for ( int j = 0 ; j < this->get_num_cols() ; j++ ) {
      data_sorted[cnt] = this->data->get_val(i,j);
      cnt++;
    }
  }
  return crawstats::median(data_sorted);  
}

double msmat::get_avg( int offset ) const  {
  assert(abs(offset) <= get_num_rts() / 2 );
  vector<double> avgs(get_num_rts() - abs(offset));
  vector<float> scan_vec(get_num_mzs());
  int num_rts = get_num_rts();
  int avg_idx = 0;
  for ( int i = offset ; i < (offset + num_rts) ; i++ ) {
    if ( i < 0 || i >= get_num_rts()  ) {
      continue;
    }
    else {
      this->get_scan(i,scan_vec);
      avgs[avg_idx++] = crawutils::sum_vect(  scan_vec ) / get_num_mzs();
    }
  }
  return (crawutils::sum_vect(avgs) / avgs.size());
}


double msmat::calc_twomat_denom( const msmat & other,
                                     double self_mean, double other_mean, int offset )  const{
                                       double this_denom = calc_mat_denom( self_mean, offset );
                                       double other_denom = other.calc_mat_denom( other_mean, offset * -1 );
                                       return sqrt(this_denom * other_denom);
}

double msmat::calc_twomat_num( const msmat & other, 
                                   double self_mean, double other_mean, int offset ) const {
                                     //offset of -2 looks like
                                     //S: 01234**
                                     //O: **01234
                                     //so, for this one, we iterate through the size, and then add offset to the index
                                     //of the other run
                                     double res = 0.0;
                                     assert(get_num_rts() == other.get_num_rts());
                                     assert(get_num_mzs() == other.get_num_mzs());
                                     vector<float> this_scan(get_num_mzs());
                                     vector<float> other_scan(get_num_mzs());
                                     for ( int self_idx = 0 ; self_idx < get_num_rts() ; self_idx++ ) {
                                       int other_idx = self_idx - offset;
                                       if ( other_idx < 0 || other_idx >= other.get_num_rts() ) {
                                         continue;
                                       }
                                       //tmp variable for comparing scans
                                       double sct = 0.0;
                                       get_scan(self_idx,this_scan);
                                       other.get_scan(other_idx,other_scan);
                                       uint num_mzs = get_num_mzs();
                                       for ( uint i = 0 ; i < num_mzs ; i++ ) {
                                         sct += (this_scan[i] - self_mean)*(other_scan[i] - other_mean);
                                       }
                                       res += sct;
                                     }
                                     return res;
}



double msmat::calc_mat_denom( double mean, int offset )  const{
  double d = 0.0;
  vector<float> scan_ref(get_num_mzs());
  int rts_num = get_num_rts();
  for ( int i = offset ; i < rts_num + offset ; i++ ) {
    if ( i < 0 || i >= get_num_rts() ) {
      continue;
    }
    else {
      this->get_scan(i,scan_ref);
      for ( uint j = 0 ; j < scan_ref.size() ; j++ ) {
        float v = scan_ref[j];
        d = d + pow((v - mean) ,2);
      }
    }
  }
  return d;
}


/* 2d correlation coefficent of the matrix shifted with another -- 
   obviously this would be faster as an FFT */
double msmat::correlation_w_shift ( const msmat & other, int shift ) const {
  const msmat & self = *this;
  /* formula :   sum ( sum ( ( ( mi - mean i ) * (sub_mj - mj_mean ) ))) ) / 
  sqrt(i_denom * j_denom )   */
  double self_avg = get_avg(shift);
  double other_avg = other.get_avg( shift * -1 );

  double denom = calc_twomat_denom( other, self_avg, other_avg, shift );
  if ( denom == 0.0 ) {
    return 0.0;
  }
  double num     =  calc_twomat_num(other, self_avg, other_avg, shift);
  return num / denom;
}


double msmat::total_current ( float threshold ) const {
  uint start_idx = 0;
  uint stop_idx  = get_num_rows() - 1;
  return total_current( start_idx, stop_idx, threshold);
}
double msmat::total_current () const {
  float threshold = 0.0f;
  return total_current(threshold);
}

void msmat::mult_val( float f ) {
  this->data->mult_val(f);
}

void msmat::div_val( float f ) {
  this->data->div_val(f);
}

void msmat::add_val( float f ) {
  this->data->add_val(f);
}




/* -------------------------------------------------------------------------------- 
msmat_scans 
-------------------------------------------------------------------------------- */

  msmat * msmat_chroms::clone () {
     msmat * m = new msmat_chroms();
     this->data->to_full_data();
     m->copy_from_msmat_deep(*this);
     return m;
  }

  msmat_chroms * msmat_chroms::extract_mzidx_range( int r_start, int r_stop ) {
    assert(r_stop < mzs.size() );
    vector<int> mz_range(r_stop - r_start);
    int foo = 5;
    for ( int i = 0;  i < r_stop - r_start ; i++ ) {
      mz_range[i] = r_start + i;
    }
    return _extract_mzidx_range(mz_range);  
  };


  


matrix_dim msmat_scans::get_data_dim() const {
  matrix_dim d = matrix_dim();
  d.row_num = get_num_rts();
  d.row_size = get_num_mzs();
  return d;
}

matrix_dim msmat_chroms::get_data_dim() const {
  matrix_dim d = matrix_dim();
  d.row_num = get_num_mzs();
  d.row_size = get_num_rts();
  return d;
}

void msmat_scans::set_scan( vector<float> & scan_data, int scan_idx ) {
  set_scan( scan_data, scan_idx , this->data);
}

void msmat_scans::set_scan( vector<float> & scan_data, int scan_idx, LazyMatrix * d ) {
  assert(scan_data.size() == d->cols() );
  d->set_row(scan_data,scan_idx);
}

void msmat_scans::zero_pixel(int chrom_idx, int scan_idx) {
   this->data->set_val(scan_idx,chrom_idx,0.0f);
}


void msmat_scans::set_chrom( vector<float> &  chrom_data, int chrom_idx ) {
  set_chrom(chrom_data, chrom_idx, this->data);
}
void msmat_scans::set_chrom( vector<float> &  chrom_data, int chrom_idx , LazyMatrix * d) {
  assert(chrom_data.size() == d->rows() );
  d->set_col(chrom_data,chrom_idx);

}

void msmat_chroms::set_scan( vector<float> & scan_data, int scan_idx ) {
  set_scan(scan_data, scan_idx, this->data);
}

void msmat_chroms::set_scan( vector<float> & scan_data, int scan_idx, LazyMatrix * d ) {
  assert(scan_data.size() == d->rows() );
  d->set_col(scan_data,scan_idx);
};


void msmat_chroms::set_chrom( vector<float> &  chrom_data, int chrom_idx ) {
  set_chrom(chrom_data, chrom_idx, this->data);
}



void msmat_chroms::set_chrom( vector<float> &  chrom_data, int chrom_idx, LazyMatrix * d ) {
  assert(chrom_data.size() == d->cols());
  d->set_row(chrom_data,chrom_idx);
}


void msmat_scans::get_scan( int scan_idx, vector<float> & inv ) const {
  data->get_row(inv, scan_idx);
}


void msmat_scans::get_chrom( int mz_idx, vector<float> & inv ) const {
  msmat * volatile_this;
  //change the to_full_data to just modify the internal data object, keeping this const.. ?
  this->to_full_data();
  inv.resize(data->rows());
  data->get_col(inv, mz_idx);
}

vector<float> msmat_scans::get_scan( int scan_idx ) const {
  vector<float> v(data->cols());
  data->get_row(v,scan_idx);
  return v;
}


vector<float>  msmat_scans::get_chrom( int mz_idx ) const {
  vector<float> v(data->rows());
  data->get_col(v,mz_idx);
  return v;
}

void msmat_chroms::get_scan( int scan_idx, vector<float> & inv ) const {
  this->to_full_data();
  inv.resize(data->rows());
  data->get_col(inv,scan_idx);
}

void msmat_chroms::get_chrom( int mz_idx, vector<float> & inv ) const {
  inv.resize(data->cols());
  data->get_row(inv, mz_idx);
}

vector<float> msmat_chroms::get_scan( int scan_idx ) const {
  vector<float> v(data->rows());
  data->get_col(v,scan_idx);
  return v;
}

vector<float>  msmat_chroms::get_chrom( int mz_idx ) const {
  vector<float> v(data->cols());
  data->get_row(v,mz_idx);
  return v;
}

void msmat_chroms::zero_pixel(int chrom_idx, int scan_idx) {
   this->data->set_val(chrom_idx,scan_idx,0.0f);
}


/*
1. As you swap indices, swap scans
2. OR, record the order of swaps from fisher-yates shuffle ( more generic )


*/

void msmat::scan_scramble() {
  std::vector<int> scan_idxs( get_num_rts());
  for ( int i = 0 ; i < get_num_rts() ; i++ ) {
    scan_idxs[i] = i;
  }
  std::vector< std::pair<int, int> > scram_cmds = crawutils::fyates_shuffle_commands(get_num_rts());
  std::vector<float> tmp_scan1(get_num_mzs());
  std::vector<float> tmp_scan2(get_num_mzs());
  for ( int i = 0 ; i < scram_cmds.size() ; i++ ) {
    get_scan( scram_cmds[i].first, tmp_scan1 );
    get_scan( scram_cmds[i].second, tmp_scan2 );
    set_scan( tmp_scan1, scram_cmds[i].second );
    set_scan( tmp_scan2, scram_cmds[i].first );
  }

}

/*
* test_ms_matrix ( int argc, char ** argv ) 
* 1. open and write msmat file
* 2. open and smooth msmat file
* 3. smooth and resample msmat file
*/


void convert_type ( msmat * in_mat, msmat * out_mat ) {

}

msmat * test_msmat ( int argc , char ** argv ) {
  //assert(argc == 2);
  msmat * msmat;
  //msmat * sg_msmat;
  char * fname = argv[1];
  char fname_o[128];
  char fname_ti[128];
  char fname_sg[128];
  crawutils::file_info f_info(fname);

  msmat_open( &msmat, f_info);

  msmat->report_matrix_type();

  return(msmat);

  strcpy(fname_o,fname);
  strcpy(fname_ti,fname);
  strcpy(fname_sg,fname);
  strcat(fname_o,".out.msmat");
  strcat(fname_ti,".ti.msmat");
  strcat(fname_sg,".sg.msmat");

  FILE * o_file, *ti_file, *sg_file;

  o_file  = fopen(fname_o, "wb");
  ti_file = fopen(fname_ti,"wb");
  sg_file = fopen(fname_sg,"wb");

  msmat->write_to_file(o_file);
  fclose(o_file);
  cerr << "wrote identical msmat back" << endl;

  //msmat_scans interpolated;
  //interpolated.copy_from_msmat_nodata(*msmat);
  //msmat->time_interpolate(0.01,interpolated);
  //cerr << "interpolated type:";
  //interpolated.report_matrix_type();
  //cerr << endl;
  //interpolated.write_to_file(ti_file);
  //fclose(ti_file);
  /* check on deletions */
  //delete(msmat);

  /*	msmat_scans * sg_msmat = new msmat_scans;
  sg_msmat->copy_from_msmat_nodata(interpolated);
  interpolated.chrom_smooth(11,2,*sg_msmat);
  sg_msmat->write_to_file(sg_file);
  fclose(sg_file);*/
  //sg_msmat = msmat; 
  //return sg_msmat;
  return msmat;
}





/* 
* 
* 	


FILE * msmat_out_file;
FILE * msmat_interp_out_file;
*/

int msmat_open( msmat ** test_msmat_p, file_info & m_info , sparse_level s ) {

  struct stat statbuf;

  FILE * msmat_file;
  const char * fname_ref;
  string tmps = m_info.full_path();
  string old_filename = m_info.original_filename;
  string old_prefix = m_info.get_prefix();
  const char * fname = tmps.c_str();
  char tmp_fname[1024];
  tmp_fname[0] = '\0';     
  //std::cerr << "stat of " << fname << std::endl;
  if ( stat(fname,&statbuf) == 0 ) {
    fname_ref = fname;
  }
  else {
    //std::cerr << "strcpy( " << tmp_fname << "," << fname << ")" << std::endl;
    strcpy(tmp_fname,fname);
    //std::cerr << "tmp_fname is " << tmp_fname << std::endl;
    strcat(tmp_fname,".msmat");
    //std::cerr << "tmp_fname is " << tmp_fname << std::endl;
    //std::cerr << "stat of " << tmp_fname << std::endl;
    if ( stat(tmp_fname,&statbuf) == 0 ) {
      m_info.add_extension("msmat");
      fname_ref = (const char*)tmp_fname;
    }
    else {
      char * dir_buf = (char*)malloc(sizeof(char) * 512);


      if ( (dir_buf = GETCWD(dir_buf,512)) == NULL ) {
        perror("error with _getcwd");
      }
      else {
        std::cerr << "could not find msmat file in directory " << dir_buf << std::endl;
      }
      perror("stat error:\n");
      #ifdef CRAW_LOGGING
      LOGH(LOG_V1,std::cerr) << "Fatal Error: cannot find msmat file for " << tmp_fname << std::endl;
      #endif
      exit(1);
    }
  }


  if ( (msmat_file = fopen(fname_ref,"rb")) == NULL ) {

    cerr << "failed to open file " << fname << endl;
    perror("");
    return(0);
  }
  *test_msmat_p = msmat_from_fh(msmat_file, s);
  (*test_msmat_p)->f_info = new crawutils::file_info(m_info); 
  (*test_msmat_p)->set_original_prefix(old_prefix);
  (*test_msmat_p)->set_original_file_name(old_filename);

  /* lets print the base peak and tic chromatograms, compare via a plotting program */
  //output_bp_chrom( **test_msmat_p, "test.out" );

  return 1;
}


void output_bp_chrom( msmat & m, char * fname  ) {
  ofstream outf(fname);		/* one column == rts, other column = bp_chrom */
  const vector<float> rts = m.get_rts();
  const vector<float> & bp_chrom = m.get_bp_chrom(); 
  assert ( m.get_bp_chrom().size() == rts.size() );
  for ( uint i = 0 ; i < rts.size() ; i++ ) {
    outf << rts[i] << "\t"  << bp_chrom[i] << "\n";
  }
  if ( fname != NULL) {
    outf.close();
  }
}


//should make this a static function, then make the constructor private
msmat * msmat_from_fh ( FILE * msmat_file , sparse_level s ) {
  msmat_header_parser hp = msmat_header_parser();
  matrix_type_enum mat_type = hp.get_class_type(msmat_file);
  if ( mat_type == CHROMS_MATRIX ) {
    return new msmat_chroms(msmat_file,s);
  }
  else if ( mat_type == SCANS_MATRIX ) {
    return new msmat_scans(msmat_file,s);
  }
  return NULL;
}

msmat * msmat_from_thevoid 
( const std::vector<float> & mzs, const std::vector<float> & rts, std::vector < std::vector< float > > & chroms , bool as_scans ) {
  msmat * out_mat;
  if ( mzs.size() != chroms.size() ) {
    /* barf */
  }
  for ( int i = 0 ; i < chroms.size() ; i++ ){
    if ( chroms[i].size() != rts.size() ) {
      /* barf */
    }
  }

  if ( as_scans == true ) {
    out_mat = new msmat_scans();
  }
  else {
    out_mat = new msmat_chroms();
  }

  out_mat->set_mzs(mzs);
  out_mat->set_rts(rts);
  out_mat->to_full_data();
  for ( int i = 0 ; i < chroms.size() ; i++ ){
    out_mat->set_chrom( chroms[i], i );
  }
  return out_mat;
  
}






vector<float> msmat::interpolate_scan_by_rt( float new_rt, int lh_idx, int rh_idx ) {
  float lh_rt, rh_rt;
  lh_rt = rts[lh_idx];
  rh_rt = rts[rh_idx];
  assert((new_rt >= lh_rt) && (new_rt <= rh_rt));
  float gap = rh_rt - lh_rt;
  float dx = ( new_rt - lh_rt ) / gap;
  vector<float> lh_scan(get_num_mzs()), rh_scan(get_num_mzs());
  get_scan(lh_idx,lh_scan);
  get_scan(rh_idx,rh_scan);
  vector<float> new_scan(lh_scan.size());
  for ( uint i = 0 ; i < lh_scan.size() ; i++ ) {
    new_scan[i] = dx * ( rh_scan[i] - lh_scan[i]) + lh_scan[i];
  }
  return new_scan;
}

void msmat::_decimate_rts( int decimate_size ) {
  int new_num_rts = rts.size() / decimate_size;
  vector<float> tmp_rts(new_num_rts);
  int old_idx = 0;          
  for ( int i = 0 ; i < new_num_rts ; i++, old_idx += decimate_size ) {
    assert(old_idx < rts.size());
    tmp_rts[i] = rts[old_idx];
  }
  rts.resize(new_num_rts);
  copy(tmp_rts.begin(), tmp_rts.end(), rts.begin());
}


void msmat::_decimate_mzs( int decimate_size ) {
  int new_num_mzs = mzs.size() / decimate_size;
  vector<float> tmp_mzs(new_num_mzs);
  int old_idx = 0;          
  for ( int i = 0 ; i < new_num_mzs ; i++, old_idx += decimate_size ) {
    assert(old_idx < mzs.size());
    tmp_mzs[i] = mzs[old_idx];
  }
  mzs.resize(new_num_mzs);
  copy(tmp_mzs.begin(), tmp_mzs.end(), mzs.begin());
}

void msmat::decimate_mzs_by_avg(int decimate_size) {
  assert(decimate_size > 0);
  _decimate_mzs(decimate_size);

  vector< vector<float>  > tmp_chroms(decimate_size);
  matrix_dim d = this->get_data_dim();
  LazyMatrix * new_data = new LazyMatrix(d.row_num,d.row_size,LM_FULL_DATA, this->msmat_file, data_start_offset );  
  for ( uint i = 0; i < tmp_chroms.size() ; i++ ) {
    tmp_chroms[i].resize(this->get_num_rts());
  }
  int old_mz_idx = 0;
  for ( int new_cnt = 0; 
    new_cnt < mzs.size() ;
    new_cnt++, old_mz_idx += decimate_size ) {
      for ( int i = 0 ; i < decimate_size ; i++ ) {
        this->get_chrom(old_mz_idx + i,tmp_chroms[i]);
      }
      vector<float> chrom_avg(this->get_num_rts());
      crawutils::average_mat_2d( tmp_chroms, chrom_avg);
      set_chrom(chrom_avg, new_cnt, new_data);
  }
  //check the destructor to see if this is kosher?
  delete this->data;
  this->data = new_data;
}


/*  reduces the matrix by resampling every group of decimate_size scans
and averaging
*/
void msmat::decimate_rts_by_avg(int decimate_size) {
  assert(decimate_size > 0);
  _decimate_rts(decimate_size);

  vector< vector<float>  > tmp_scans(decimate_size);
  matrix_dim d = this->get_data_dim();
  LazyMatrix * new_data = new LazyMatrix(d.row_num,d.row_size, LM_FULL_DATA, this->msmat_file, data_start_offset );  
  for ( uint i = 0; i < tmp_scans.size() ; i++ ) {
    tmp_scans[i].resize(this->get_num_mzs());
  }
  int old_rt_idx = 0;
  for ( int new_cnt = 0; 
    new_cnt < rts.size() ;
    new_cnt++, old_rt_idx += decimate_size ) {
      for ( int i = 0 ; i < decimate_size ; i++ ) {
        this->get_scan(old_rt_idx + i,tmp_scans[i]);
      }
      vector<float> scan_avg(this->get_num_mzs());
      crawutils::average_mat_2d( tmp_scans, scan_avg);
      set_scan(scan_avg, new_cnt, new_data);
  }
  //check the destructor to see if this is kosher?
  delete this->data;
  this->data = new_data;
}

void msmat::decimate_mzs_by_sample(int decimate_size) {
  assert(decimate_size > 0);
  _decimate_mzs(decimate_size);
  matrix_dim d = this->get_data_dim();
  LazyMatrix * new_data = new LazyMatrix(d.row_num,d.row_size, LM_FULL_DATA, this->msmat_file, data_start_offset );  
  int old_mz_idx = 0;
  vector<float> t_chrom(this->get_num_rts());
  for ( int new_cnt = 0; 
    new_cnt < mzs.size() ;
    new_cnt++, old_mz_idx += mzs.size() ) {
      get_chrom(old_mz_idx,t_chrom);
      set_chrom(t_chrom, new_cnt, new_data);
  }
  //check the destructor to see if this is kosher?
  delete this->data;
  this->data = new_data;
}

void msmat::decimate_rts_by_sample(int decimate_size) {
  assert(decimate_size > 0);
  _decimate_rts(decimate_size);
  matrix_dim d = this->get_data_dim();
  LazyMatrix * new_data = new LazyMatrix(d.row_num,d.row_size, LM_FULL_DATA, this->msmat_file, data_start_offset );  
  int old_rt_idx = 0;
  vector<float> t_scan(this->get_num_mzs());
  for ( int new_cnt = 0; 
    new_cnt < rts.size() ;
    new_cnt++, old_rt_idx += rts.size() ) {
      get_scan(old_rt_idx,t_scan);
      set_scan(t_scan, new_cnt, new_data);
  }
  //check the destructor to see if this is kosher?
  delete this->data;
  this->data = new_data;
}


bool msmat::equals_dimensions( const msmat & other ) const {
  matrix_dim my_dim = get_data_dim();
  matrix_dim other_dim = other.get_data_dim();
  if ( my_dim == other_dim ) {
    return true;
  }
  return false;
}

bool msmat::equals_data ( const msmat & other, double tolerance )  const {
  vector<float> r(data->cols());
  vector<float> other_r(other.data->cols());
  for ( int i = 0 ; i < get_num_rows() ; i++ ) {
    data->get_row(r,i);
    other.data->get_row(other_r,i);
    assert(r.size() == other_r.size());
    for ( uint j = 0 ; j < r.size() ; j++ ) {
      if ( fabs(r[j] - other_r[j]) > tolerance ) {
        return false;
      }
    }
  }
  return true;
}


double msmat::sum_abs_diffs ( const msmat & other ) const {
  if ( ! equals_dimensions(other) ) {
    return -1.0f;
  }
  //technically I should check that the other run is also of the same class
  double total_d = 0.0;
  size_t nrows = get_num_rows();
  vector<float> r(get_num_cols());
  vector<float> other_r(get_num_cols());
  for ( int i = 0 ; i < get_num_rows() ; i++ ) {
    data->get_row(r,i);
    other.data->get_row(other_r,i);
    assert(r.size() == other_r.size());
    for ( uint j = 0 ; j < r.size() ; j++ ) {
      total_d += fabs(r[j] - other_r[j]);
    }

  }
  return total_d;

}

float msmat::rt_from_scannum( int scannum ) const {
  return rts.at(scannum);
}
float msmat::rt_from_scannum( float scannum ) const {
  int base_scannum;
  if ( scannum > rts.size() - 1 || scannum < 0 ) {
    std::cerr << "Error with scan number interpolation: scannum, rts.size()" << scannum << " , " << rts.size() << std::endl;
    base_scannum = rts.size() - 1;
  }
  else {
      base_scannum = (int)scannum;
  }
  float delt = scannum - base_scannum;
  if ( delt == 0.0 ) {
    return rts.at(base_scannum);
  }
  else {
    float rt_delt = (rts.at(base_scannum+1) - rts.at(base_scannum)) * delt;
    return rts.at(base_scannum) + rt_delt;
  }
}


void msmat::tabular_output ( std::ostream & o , int mz_idx_start , int mz_idx_stop ,
				  int rt_idx_start , int rt_idx_stop , bool labels, char delim ) {
  if ( mz_idx_stop == - 1 ) {
    mz_idx_stop = this->get_num_mzs() - 1;
  } 
  if  ( rt_idx_stop == - 1 ) {
    rt_idx_stop = this->get_num_rts() - 1;
  }
  assert( mz_idx_stop > mz_idx_start);
  assert( rt_idx_stop > rt_idx_start);
  const int num_mzs = get_num_mzs();
  const int slice_rt_size = rt_idx_stop - rt_idx_start + 1;
  const int slice_mz_size = mz_idx_stop - mz_idx_start + 1;
  vector<float> tmp_chrom(get_num_rts());
  vector<float> out_chrom(rt_idx_stop - rt_idx_start + 1);
  vector<float> rts_slice(slice_rt_size);
  vector<float> mzs_slice(slice_mz_size);
  //note copy is [start,stop) not [start,stop], hence adding 1 to the second arg.
  std::copy(this->rts.begin() + rt_idx_start , this->rts.begin() + rt_idx_stop + 1,rts_slice.begin());
  std::copy(this->mzs.begin() + mz_idx_start , this->mzs.begin() + mz_idx_stop + 1,mzs_slice.begin());


  if ( labels ) {
    o << "M/Z" << delim;
    for ( int i = 0 ; i < rts_slice.size() - 1 ; i++ ) {
      o << rts_slice[i] << delim;
    }
    o << rts_slice[rts_slice.size() - 1] << std::endl;
  }

  for ( int i = 0 ; i < mzs_slice.size() ; i++ ) {
    int mz_idx = mz_idx_start + i;

    this->get_chrom(mz_idx,tmp_chrom);
    copy(tmp_chrom.begin() + rt_idx_start , tmp_chrom.begin() + rt_idx_stop + 1 , out_chrom.begin());
    if ( labels ) {
      cerr << "mzs_slice[" << i << "] :" << mzs_slice[i];
      o << 0.0 + mzs_slice[i] << delim;
    }
    for ( int j = 0 ; j < out_chrom.size() - 1 ; j++ ) {
      o << out_chrom[j] << delim;
    }
    o << out_chrom[out_chrom.size() - 1] << std::endl;
    o.flush();
  }
}
void msmat::tabular_output ( std::ostream & o ,  float mz_start, float mz_stop , float rt_start , float rt_stop, bool labels, char delim ) {
  int mz_start_idx, mz_stop_idx, rt_start_idx, rt_stop_idx;
  std::cerr << "hello!" << std::endl;
  mz_start_idx = this->mz_to_mzidx(mz_start);
  mz_stop_idx = this->mz_to_mzidx(mz_stop);
  rt_start_idx = this->rt_to_rtidx(rt_start);
  rt_stop_idx = this->rt_to_rtidx(rt_stop);
  std::cerr << "converted!" << std::endl;
#ifdef CRAW_LOGGING
  LOGH(LOG_DEBUG,std::cerr) << mz_start << " " << mz_stop << " " << rt_start << " " << rt_stop << std::endl;
  LOGH(LOG_DEBUG,std::cerr) << mz_start_idx << " " << mz_stop_idx << " " << rt_start_idx << " " << rt_stop_idx << std::endl;
#endif
  tabular_output(o,mz_start_idx,mz_stop_idx,rt_start_idx,rt_stop_idx,labels,delim);
}

void msmat::adjust_ort_wrt_map_linear ( float rt ) {
  
    flt_map::iterator ort_k = this->ort_wrt_map.begin();
    flt_map::iterator ort_end = ort_wrt_map.end();
    while ( ort_k != ort_end ) {
       (*ort_k).second = (*ort_k).second + rt;
    }
}
void msmat::adjust_ort_wrt_map_linear ( int nscans ) {
    this->adjust_ort_wrt_map_linear( nscans * (rts[1] - rts[0]) );
}

/* ---------------END MSMAT -------- */

//TODO -- make the above step more efficient


void msmat_scans::get_scan_interp( float scan_idx, vector<float> & inv ) const {
  if ( scan_idx > get_num_rts() - 1 ) {
    return get_scan((int)floorf(scan_idx),inv);
  }
  
  if ( (scan_idx >= roundf(scan_idx) - 1E-3) && (scan_idx <= roundf(scan_idx) + 1E+3 ) ) {
    return get_scan((int)roundf(scan_idx),inv);
  }
  vector<float> s1(inv.size());
  vector<float> s2(inv.size());
  data->get_row(s1,(int)floorf(scan_idx));
  data->get_row(s2,(int)ceilf(scan_idx));
  float delta = scan_idx - ceilf(scan_idx);
  for ( int i = 0 ; i < s1.size() ;i++ ) {
    inv[i] = s1[i] + ( s2[i] - s1[i] ) * delta;
  }
}

vector<float>  msmat_scans::get_scan_interp( float scan_idx ) const {
  vector<float> out_v(get_num_mzs());
  get_scan_interp(scan_idx,out_v);
  return out_v;
}

void msmat_chroms::get_scan_interp( float scan_idx, vector<float> & inv ) const {
  if ( scan_idx > get_num_rts() - 1 ) {
    return get_scan((int)floorf(scan_idx),inv);
  }
  
  if ( (scan_idx >= roundf(scan_idx) - 1E-3) && (scan_idx <= roundf(scan_idx) + 1E+3 ) ) {
    return get_scan((int)roundf(scan_idx),inv);
  }
  vector<float> s1(inv.size());
  vector<float> s2(inv.size());
  data->get_row(s1,(int)floorf(scan_idx));
  data->get_row(s2,(int)ceilf(scan_idx));
  float delta = scan_idx - ceilf(scan_idx);
  for ( int i = 0 ; i < s1.size() ;i++ ) {
    inv[i] = s1[i] + ( s2[i] - s1[i] ) * delta;
  }
}

vector<float>  msmat_chroms::get_scan_interp( float scan_idx ) const {
  vector<float> out_v(get_num_mzs());
  get_scan_interp(scan_idx,out_v);
  return out_v;
}



void ttest_score_matrix::unmask_search ( mz_type mz, rt_type start_rt, rt_type stop_rt ) {
  int start_rt_idx = this->rt_to_rtidx(start_rt);
  int stop_rt_idx = this->rt_to_rtidx(stop_rt);
  unmask_search( mz, start_rt, stop_rt );
}

void ttest_score_matrix::unmask_search ( mz_type mz, int start_rt_idx, int stop_rt_idx ) {
  int mz_idx = this->mz_to_mzidx(mz);
  this->unmask_search( mz_idx, start_rt_idx, stop_rt_idx);
}

void ttest_score_matrix::unmask_search ( int mz_idx, int start_rt_idx, int stop_rt_idx ) {
  for ( int i = start_rt_idx ; i <= stop_rt_idx ; i++ ) {
    //note that ttest_score_matrix is derived from chroms, so I can
    //set mz_idx as the chromatogram value
    this->data->set_val(mz_idx,i,TTEST_SEARCH);
  }  
}


#ifdef HAVE_CRAWDAD


void msmat::chrom_smooth( ChromSmoother & sg ) {
 vector<float> tmp_in_chrom(get_num_rts());
 vector<float> tmp_out_chrom(get_num_rts());
 for ( uint i = 0 ; i < this->get_num_mzs() ; i++ ) {
     get_chrom(i,tmp_in_chrom);
     sg.smooth_vect(tmp_in_chrom,tmp_out_chrom);
     set_chrom(tmp_out_chrom,i);
 }
}


void msmat::spike_filter( ChromSmoother & sg ) {
 vector<float> tmp_in_chrom(get_num_rts());
 vector<float> tmp_out_chrom(get_num_rts());
 for ( uint i = 0 ; i < this->get_num_mzs() ; i++ ) {
     get_chrom(i,tmp_in_chrom);
     sg.spike_filter(tmp_in_chrom,tmp_out_chrom);
     set_chrom(tmp_out_chrom,i);
 }
}





void msmat::split_by_num_chroms ( int slice_size, vector< msmat_scans * > & v ) {
  for ( int i = 0 ; i < get_num_mzs() ; i += slice_size ) {
    msmat_scans * c = new msmat_scans();
    vector<int> mz_idxs(0);
    crawutils::create_range( i , std::min((int)i+slice_size,(int)get_num_mzs() ) , 1 ,mz_idxs );
    c->copy_from_msmat_mzs( *this , mz_idxs );
    v.push_back(c);
  }
}
void msmat::split_by_num_segments ( int num_segments, vector< msmat_scans * > & v ) {

  int segment_size = get_num_mzs() / num_segments;
  split_by_num_chroms( segment_size, v );
}



double msmat::avg_diag_score ( msmat & other , SpectraCompFunc comp_f ) const {
  return total_diag_score(other,comp_f) / get_num_rts();
}

double msmat::avg_diag_score ( msmat & other, const int start_scan_idx, const int stop_scan_idx, SpectraCompFunc comp_f ) const {
  return total_diag_score(other,start_scan_idx,stop_scan_idx,comp_f) / get_num_rts();
}


double msmat::total_diag_score( msmat & other, SpectraCompFunc comp_f ) const {
  int last_idx = get_num_rts() - 1;
  return total_diag_score( other, 0 , last_idx, comp_f );
}

double msmat::total_diag_score( msmat & other , const int start_scan_idx , const int stop_scan_idx, SpectraCompFunc comp_f ) const {
  assert(get_num_rts() == other.get_num_rts());
  double t = 0.0;
  vector<float> this_scan(get_num_mzs()) , other_scan(get_num_mzs());
  
  for ( int s = start_scan_idx ; s < stop_scan_idx + 1 ; s++ ) {
    get_scan(s,this_scan);
    other.get_scan(s,other_scan);
    t += comp_f( this_scan, other_scan , crawutils::SQRSUM_NULL, crawutils::SQRSUM_NULL);
    if ( s % 100 == 0 ) {
      std::cerr << ".";
    }
    std::cerr << std::endl;
  }
  return t;
}

std::pair<float,float> msmat::weighted_diag_score ( msmat & other , 
      const int start_scan_idx, const int stop_scan_idx, const SpectraCompareConfig & c ) const {
   assert(get_num_rts() == other.get_num_rts());
   float total = 0.0;
   float weights = 0.0;
   vector<float> this_scan(get_num_mzs()) , other_scan(get_num_mzs());
   for ( int s = start_scan_idx ; s < stop_scan_idx + 1 ; s++ ) {
    get_scan(s,this_scan);
    other.get_scan(s,other_scan);
    float sim_score = c.comp_func(this_scan, other_scan, -1.0, -1.0);
//    if ( sim_score > 1.0f ) {
//      std::cerr << sim_score << std::endl;
//    }
    float weight   = c.weight_func(this_scan, other_scan);

    if ( sim_score != 0.0f ) {
      total += sim_score * weight;
      weights += weight;
    }
    

  }
  return std::pair<float,float>(total,weights);                
}

std::pair<float,float> msmat::bp_weighted_diag_score ( msmat & other, const int start_scan_idx, const int stop_scan_idx, SpectraCompFunc comp_f ) const {

  float total = 0.0;
  float weights = 0.0;
  vector<float> this_scan(get_num_mzs()) , other_scan(get_num_mzs());
  std::cerr << "stop_scan_idx" << stop_scan_idx << std::endl;
  float logMaxI = 0.0f;
  for ( int s = start_scan_idx ; s < stop_scan_idx + 1 ; s++ ) {
    get_scan(s,this_scan);
    other.get_scan(s,other_scan);

    float sim = comp_f(this_scan, other_scan, crawutils::SQRSUM_NULL, crawutils::SQRSUM_NULL) * logMaxI;
    if ( sim > 0 ) {
      total += sim;
      weights += logMaxI;
    }

  }
  return std::pair<float,float>(total,weights);  
}

void msmat::mask_based_on_ppids ( ppid_row_parser & p, float min_retain_i ) {
    /* chrom_idx -> pair_num -> ( start_rt, stop_rt ) */
    this->summary_data_ok = false;
    std::vector< std::vector < std::pair< int, int > > > ppids_by_chrom_rt(get_num_mzs());
    std::vector< float > tmp_chrom(get_num_rts(),0.0f);
    int num_ppids = p.rows.size();
    for ( int i = 0; i < num_ppids ; i++ ) {
       bool skip_flag = false;
       ppid_row r = p.rows[i];
       int chrom_idx, start_rt_idx, stop_rt_idx, charge;
       //base_chrom_idx = this->mz_to_mzidx(r.mono_iso_mz);
       try {
	 start_rt_idx = this->rt_to_rtidx(r.rt_start);
	 stop_rt_idx = this->rt_to_rtidx(r.rt_stop);
       } catch ( RTRangeError r ) {
         skip_flag = true;
       }
       if ( skip_flag ) continue;
       for ( int m = 0 ; m < MASK_ISOTOPE_PEAKS ; m++ ) {
         float mi_start_rt_idx, mi_stop_rt_idx;
         mi_start_rt_idx = start_rt_idx;
         mi_stop_rt_idx  = stop_rt_idx;
          float cand_mz = r.observed_mono_mz + ( 1.0 / r.charge ) * m;
          int chrom_idx;
          if ( cand_mz < this->mzs[0] || cand_mz >= this->mzs[this->mzs.size() - 1] ) {
            continue;
          }
          try { 
              chrom_idx = this->mz_to_mzidx(cand_mz);
          } catch ( MZRangeError m ) {
             skip_flag = true;
          }
          if ( skip_flag ) continue;
          
          
          if ( min_retain_i > 0.0f ) {

      get_chrom(chrom_idx,tmp_chrom);
            for ( int i = start_rt_idx ; i >= 0 ; i-- ) {
              if (  tmp_chrom[i] >= min_retain_i ) {
                mi_start_rt_idx = i;
              }
              else {
                break;
              }
            }
            for ( int i = stop_rt_idx ; i <= get_num_rts() ; i++ ) {
              if ( tmp_chrom[i] >= min_retain_i ) {
                mi_stop_rt_idx = i;
              }
              else {
                break;
              }
            }
          }
      ppids_by_chrom_rt[chrom_idx].push_back( std::pair<int,int>(mi_start_rt_idx,mi_stop_rt_idx ));
       }
    }
    assert(ppids_by_chrom_rt.size() == get_num_mzs());
    for ( int chrom_idx = 0 ; chrom_idx < ppids_by_chrom_rt.size() ; chrom_idx++ ) {
        int num_ppids = ppids_by_chrom_rt[chrom_idx].size();
        if ( num_ppids > 0 ) {
           this->retain_chrom_set(chrom_idx, ppids_by_chrom_rt[chrom_idx]);
        }
        else {
           this->zero_chrom(chrom_idx);
     }
    }
}

void msmat::mask_based_on_hardklor ( hardklor_row_parser & p ) {
    /* chrom_idx -> pair_num -> ( start_rt, stop_rt ) */
    this->summary_data_ok = false;
    //initialize with a vector for _each_ chromatogram - gets a bit memory intensive, oh well
    std::vector< std::vector < int > > retain_by_chrom_rt(mzs.size());
    int num_dists = p.pep_lines.size();
    for ( int i = 0; i < num_dists ; i++ ) {
       bool skip_flag = false;
       float rt_idx;
       hardklor_pep_row r = p.pep_lines[i];
       int chrom_idx, start_rt_idx, stop_rt_idx, charge;
       //base_chrom_idx = this->mz_to_mzidx(r.mono_iso_mz);
       try {
           rt_idx = this->rt_to_rtidx(r.retention_time);
       } catch ( RTRangeError r ) {
         skip_flag = true;
       }
       if ( skip_flag ) continue;
       /* output the next five isotopic flags */
       for ( int m = 0 ; m < MASK_ISOTOPE_PEAKS ; m++ ) {

           double cand_mz = r.base_mz + ( 1.0 / r.charge_state ) * m;
           int chrom_idx;
           if ( cand_mz < this->mzs[0] || cand_mz >= this->mzs[this->mzs.size() - 1] ) {
               continue;
           }
           try { 
               chrom_idx = this->mz_to_mzidx(cand_mz);
           } catch ( MZRangeError m ) {
               skip_flag = true;
           }
/*          if ( skip_flag ) continue;
           get_chrom(chrom_idx,tmp_chrom);
           for ( int i = start_rt_idx ; i >= 0 ; i-- ) {
               if (  tmp_chrom[i] >= min_retain_i ) {
                   mi_start_rt_idx = i;
               }
               else {
                   break;
               }
           }
           for ( int i = stop_rt_idx ; i <= get_num_rts() ; i++ ) {
               if ( tmp_chrom[i] >= min_retain_i ) {
                   mi_stop_rt_idx = i;
               }
               else {
                   break;
               }
           }
           ppids_by_chrom_rt[chrom_idx].push_back( std::pair<int,int>(mi_start_rt_idx,mi_stop_rt_idx ));*/
          retain_by_chrom_rt[chrom_idx].push_back(rt_idx);
       }
    }
    for ( int chrom_idx = 0 ; chrom_idx < retain_by_chrom_rt.size() ; chrom_idx++ ) {
        int num_rts = retain_by_chrom_rt[chrom_idx].size();
        if ( num_rts > 0 ) {
           //this->retain_chrom_set(chrom_idx, ppids_by_chrom_rt[chrom_idx]);
           this->retain_dist_set( chrom_idx, retain_by_chrom_rt[chrom_idx]);
        }
        else {
           this->zero_chrom(chrom_idx);
        }
    }
}



    coord_opts::coord_opts( const crawusage::optset & o ) {
      start_rt = stop_rt = start_mz = stop_mz = -1.0f;
      opts = o;
      std::string opt_val;
      if ( opts.get_val("start_rt",opt_val) ) {
         start_rt = atof(opt_val.c_str());
      }
      if ( opts.get_val("stop_rt",opt_val) ) {
         stop_rt = atof(opt_val.c_str());
      }
      if ( opts.get_val("start_mz",opt_val) ) {
         start_mz = atof(opt_val.c_str());
      }
      if ( opts.get_val("stop_mz",opt_val) ) {
         stop_mz = atof(opt_val.c_str());
      }
    }

    void coord_opts::update_from_msmat ( const msmat & m ) {
        if ( start_rt == -1.0f ) {
           start_rt = m.rts[0];
        }
        if ( stop_rt == -1.0f ) {
	  stop_rt = m.rts[m.get_num_rts()-1];
        }
        if ( start_mz == -1.0f ) {
           start_mz = m.mzs[0];
        }
        if ( stop_mz == -1.0f ) {
           stop_mz = m.mzs[m.get_num_mzs()-1];
        }
    }

#endif
