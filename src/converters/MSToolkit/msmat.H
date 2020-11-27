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
/* $Id: msmat.H,v 1.7 2009/02/23 07:20:30 gfinney Exp $ */

#ifndef _INCLUDE_MS1_MATRIX
#define _INCLUDE_MS1_MATRIX 

#define MASK_ISOTOPE_PEAKS 5

#define MSMAT_MZ_TYPE float
#define MSMAT_RT_TYPE float

// STL data structures
#include <string>
#include <vector>
#include <map>
//#include <algorithm>

// I/O
#include <stdio.h>
#include <iostream>
#include <fstream>

// Misc.
//#include <assert.h>
#include <sys/types.h>
#include <sys/stat.h>

#ifndef _MSC_VER
#include <unistd.h>
#endif



#include "msmat_base.H"
#include "LazyMatrix.H"
#include "crawutils.H"



#ifdef HAVE_CRAWDAD
#include "SpectraCompareConfig.H"
#include "PPIDReader.H"
#include "filters/ChromSmoother.H"
#include "CrawUsage.H"
#endif


//TODO SG implementation -- find matlab code in downloads folder
//savitzky-golay implementation -- need to find citation -- look 

using namespace crawutils;

static const float TTEST_NO_SEARCH = 1.0f;
static const float TTEST_SEARCH = 0.0f;

class msmat_scans;
class msmat_chroms;
class msmat_header_parser;



class MZRangeError {
public :
   enum MZError { MZLowError , MZHighError };
   MZError errflag;
   MZRangeError( MZError m ) {
      errflag = m;
   }
};


class RTRangeError {
public :
   enum RTError { RTLowError , RTHighError };
   RTError errflag;
   RTRangeError( RTError r ) {
      errflag = r;
   }
};



class matrix_dim {
  public :
	matrix_dim(size_t rs, size_t rn) {
		row_size = rs;
		row_num = rn;
	}
	matrix_dim() {
		row_size = 0;
		row_num = 0;
	}
    size_t row_size, row_num;
    ///returns true when row_size and row_num are equal
    bool operator==(matrix_dim & d)  {
    	if( row_size == d.row_size && row_num == d.row_num ) { return true; }
    	else { return false ; }
    }
  
};


static MSMAT_MZ_TYPE TRANSITION_MATCH_MZ_TOL = 1E-6;



struct transition_info {


   static std::pair<MSMAT_MZ_TYPE, MSMAT_MZ_TYPE> trans_from_str ( const std::string & s );

   MSMAT_MZ_TYPE precursor_mz;
   MSMAT_MZ_TYPE product_mz;
   int mz_idx;

   bool matches ( const std::pair< MSMAT_MZ_TYPE, MSMAT_MZ_TYPE> & p ) const {
      if ( (fabs(precursor_mz - p.first ) < TRANSITION_MATCH_MZ_TOL ) &&
          (fabs(product_mz   - p.second ) < TRANSITION_MATCH_MZ_TOL ) ) {
        return true;
      }
      return false;
   }
  
   void init() {
     precursor_mz = -1.0;
     product_mz   = -1.0;
     mz_idx = -1;
   };

   transition_info () {
     init();
   }

   transition_info( std::string & trans_str, int mz_idx ) {
       std::pair< MSMAT_MZ_TYPE, MSMAT_MZ_TYPE> mz_p = trans_from_str(trans_str);
       precursor_mz = mz_p.first;
       product_mz   = mz_p.second;
       mz_idx       = mz_idx; 
   }

};

///stores information about transitions, and transitions organized by precursor and product
struct transition_manager {
    std::vector<transition_info> transitions;
    std::vector<MSMAT_MZ_TYPE> precursors;
    std::vector< std::vector< int > > products_by_transitions;

    std::vector<transition_info> get_trans_by_precursor( MSMAT_MZ_TYPE precursor );
    std::vector<transition_info> get_trans_by_precursor_idx (int idx);
    transition_manager ( std::vector<transition_info> transitions );
    transition_manager() {};

};


class msmat {
  friend class msmat_header_parser;
  #ifdef HAVE_CRAWDAD
  friend class coord_opts;
  #endif


//protected:
public:
  
  float bin_size;

  off_t data_start_offset;
  FILE * msmat_file;

  vector<transition_info> transitions;
  transition_manager t_manager;
  vector<float> rts;
  vector<std::string> labels;
  flt_map ort_wrt_map;
  bool summary_data_ok;
  sparse_level sp_level;
  msmat_header_parser * hp;
  std::string original_file_name;
  std::string original_prefix;
  std::map< std::string, float > mzstr_id_map;
  mutable bool loaded_full_data;
  crawutils::file_info * f_info;

  
public:
  /* init and maintenance functions */

  vector<float> bp_chrom, tic_chrom;	

  // this should really be protected but for the moment I have
  // MSMATSigDist modifying it
  LazyMatrix * data;  


  typedef MSMAT_RT_TYPE rt_type;
  typedef MSMAT_MZ_TYPE mz_type;
  typedef vector<rt_type> rts_type;
  typedef vector<mz_type> mzs_type;
  vector<mz_type> mzs;

  /* -- END DATA -- */


public:

    void init_transitions_from_labels();
    void init_transition_manager();

  void to_full_data ( bool force = false ) const ;
  
  void set_min_val ( float v = 0.0f );
  void set_max_val ( float v );

  void init_self_from_header();

  void set_sparse_level( sparse_level sp ) {
    sp_level = sp;
  }
  sparse_level get_sparse_level () {
    return sp_level;
  }
  void close_file() {
    if (msmat_file != NULL) {
      fclose(msmat_file);
    }	
  }    


  void init_matrix_from_file ( FILE * msmat_fh, sparse_level l = LM_FULL_DATA);
  void write_matrix_to_file ( FILE * msmat_fh );
  void reset_from_header() {
    if ( data != NULL ) {
      delete(data);
      data = NULL;
    }
    init_self_from_header();
  };

  

  double total_current ( float start_rt, float stop_rt, float threshold ) const;
  double total_current (  uint start_scan , uint last_scan,  float threshold ) const ;
  double total_current ( float threshold ) const;
  double total_current () const;

  void mult_val( float f );
  void div_val( float f );
  void add_val( float f );
  void retain_mzs_nodata( const vector<int> mzidxs );
  void copy_from_msmat_mzs( msmat & t, const vector<int> mzidxs );

  void decimate_rts_by_sample( int decimate_size );
  void decimate_rts_by_avg      (  int decimate_size );
  void decimate_mzs_by_sample( int decimate_size );
  void decimate_mzs_by_avg      (  int decimate_size );
  void _decimate_mzs ( int decimate_size);
  void _decimate_rts( int decimate_size );

  float rt_from_scannum( int scannum ) const;
  float rt_from_scannum( float scannum ) const;

  void tabular_output ( std::ostream & o , int mz_idx_start = 0 , int mz_idx_stop = -1 ,
			int rt_idx_start = 0 , int rt_idx_stop = -1 , bool labels=true, char delim='\t' );

  void tabular_output ( std::ostream & o , float mz_start , float mz_stop,
			float rt_start ,float rt_stop , bool labels=true, char delim='\t' );

///reads data in from a file, replacing the current data. This is not designed to 
///be a way to _create_ a file from scratch, but rather to do offline processing and
///to load in the raw data
  void tabular_input  ( std::ifstream & is );
  void float_pseudo_matrix_input ( std::istream & i ) {
  
  }

//protected :
public:
  /* internal calculations */

  void clear_data();
  void interpolate_new_scan( const vector<float> & lh_scan,
			     const vector<float> & rh_scan, 
			     vector<float> & interpolated_scan, 
			     float xdelt ) const;
  void calc_bp_chrom();
  void calc_tic_chrom();
  void calc_summary_data();

  /* modification functions */
  //TODO -- figure out why not-callable virtual functions fail
  virtual void _trim_start( int nscans ) = 0 ; //{throw "virt hell";};
  virtual void _trim_stop ( int nscans ) = 0 ;
  virtual void _trim_start_mzs( int  ) = 0 ; //{throw "virt hell";};
  virtual void _trim_stop_mzs ( int  ) = 0;
  void _trim_start_rts ( int nscans );
  void _trim_stop_rts  ( int nscans );
  void _trim_start_mzs_list ( int nscans ) ;
  void _trim_stop_mzs_list ( int nscans ) ;

  
  float ort_to_wrt ( float lookup_ort );
  float wrt_to_ort ( float lookup_ort );

  void get_mode_rt_intervals(float granularity, float lb, float rb, std::vector<uint> & cnt ) const;
  double get_mode_rt_interval() const;
  double get_mode_rt_interval(double roundto) const;

  virtual ~msmat() {} ;
  void init_msmat();

  msmat() {
    init_msmat();
    //init_data();
  };

  /* deep copy constructor */
  msmat( const msmat & o ) {
    init_msmat();
    this->copy_from_msmat_deep(o);
  }

  //ensure that this remains a virtual base class
  virtual msmat * clone() = 0;

  msmat * copy_no_data();
  /// returns index in the msmat file for an mz. Note that it returns the left-hand (i.e. lt ) value, rather than the closest
  int mz_to_mzidx ( float mz ) const;
  int rt_to_rtidx ( float rt ) const;

  float mzidx_to_mz ( int mz_idx ) const {
     return mzs[mz_idx];
  };
  float rtidx_to_rt ( int rt_idx ) const {
     return rts[rt_idx];
  };


  //HACK -- we assume that a) bins are evenly spaced and
  //HACK -- b) that rts[1] - rts[0] reflects that spacing
  float get_rt_interval() const { return rts[1] - rts[0]; };

  void init_from_msmat_file (FILE * msmat_f, sparse_level l);
  void init_from_msmat_file( FILE * msmat_f );
  msmat( FILE* msmat_f );
  void copy_from_msmat_nodata( const msmat & msmat);
  void copy_from_msmat( const msmat & msmat );
  void copy_from_msmat_deep( const msmat & msmat );
  void set_mzs ( const vector<float> & mzs );
  void set_mzs ( float * in_mzs, int num_mzs );
  void set_rts ( const vector<float> & rts );
  void set_rts ( float * in_rts, int num_rts );
  void set_ort_wrt_map ( flt_map & ow_map );
  void set_ort_wrt_map ( const vector<float> & ort,
			 const vector<float> & wrt );
  string get_original_prefix() const {
     return this->original_prefix;
  }
  void set_original_prefix( std::string s ) {
     this->original_prefix = s;
  }
  string get_original_file_name() const {
     return original_file_name;
  }
  void set_original_file_name( string s) {
     original_file_name = s;
  }
  size_t get_num_mzs () const {
    return mzs.size();
  }
  size_t get_num_rts () const {
    return rts.size();
  }
  const flt_map & get_ort_wrt_map () const{
    return ort_wrt_map;
  }
  void set_bp_chrom ( vector<float> i_bp_chrom ) {
    this->bp_chrom = i_bp_chrom;
  }
  void set_tic_chrom ( vector<float>  i_tic_chrom) {
    this->tic_chrom = i_tic_chrom;
  }
  const vector<float> & get_bp_chrom() const {
    return this->bp_chrom;
  }
  const vector<float> & get_tic_chrom() const {
    return this->tic_chrom;
  }
  vector<float> get_rts() const;
  vector<float> get_mzs() const;
  const vector<float> & get_rts_vect() const {
    return rts;
  }
  const vector<float> & get_mzs_vect() const {
    return mzs;
  }
  float get_bin_size() const {
    return bin_size;
  }
  float get_scandata_len() const {
	  return rts[rts.size() - 1] - rts[0];
  }
  int check_size_from_position() const;
  vector<float> interpolate_scan_by_rt( float new_rt, int lh_idx, int rh_idx );
  vector<float> interpolate_scan_by_idx( float new_idx );
  void write_to_file ( FILE * msmat_fh );
  int get_scan_idx_by_rt(float rt) const;
  int get_chrom_idx_by_mz(float mz) const;
  void trim_mzs ( float keep_start_mz , float keep_stop_mz);
  void trim_mzs_start ( float keep_start_mz );
  void trim_mzs_stop ( float keep_stop_mz );
  void trim_mzs_data_bound ( int start_idx, int stop_idx );
  void trim_mzs_list_bound ( int start_idx , int stop_idx);
  void trim_start( int nscans );
  void trim_stop( int nscans);
  void trim_start( rt_type rt );
  void trim_stop ( rt_type rt );
  void trim_start_stop( int start_scans, int stop_scans);
  void trim_start_stop( rt_type start_rt, rt_type stop_rt);
  void time_interpolate(double interval, msmat & new_matrix);
  //  void chrom_smooth( ChromSmoother & sg );
  //  void spike_filter ( ChromSmoother & spike_filter );

  
  void init_data();
  void pad_start_stop( int start_mapped_coord, 
		int stop_mapped_coord );

  double get_avg() const;
  double get_avg ( int offset ) const;
  double get_median () const;
  double calc_twomat_denom( const msmat & other,double self_mean, double other_mean, int offset) const;
  double calc_mat_denom ( double mean, int offset) const;
  double correlation_w_shift ( const msmat & other, int shift) const; 
  double calc_twomat_num( const msmat & other, double self_mean, double other_mean, int offset ) const;
  
  bool equals ( const msmat & other ) const;
  bool equals_header ( const msmat & other ) const;
  bool equals_dimensions ( const msmat & other ) const;
  bool equals_data   ( const msmat & other, double tolerance=0.000001 ) const;
  double sum_abs_diffs ( const msmat & other ) const;


  virtual void zero_pixel ( int chrom_idx, int scan_idx ) = 0;
  void zero_chrom ( int chrom_idx , float start_rt = -1 , float stop_rt = -1);
  void zero_chrom ( int chrom_idx , int start_scan_idx , int stop_scan_idx );
  void zero_chrom ( float mz , float start_rt = -1 , float stop_rt = -1);
  void retain_chrom ( int chrom_idx, int start_scan_idx, int stop_scan_idx );
  void retain_chrom_set ( int chrom_idx , std::vector< std::pair< int , int > > & v );
  void retain_dist_set ( int chrom_idx , std::vector< int > &  v );
 


  void warp_to_template( const vector<float> & t_coords,
			 const vector<float> & a_coords, 
			 const msmat & templ, 
			 const msmat & align, int linear_scans_offset = 0);
  void warp_average_new_output (const vector<float> & t_coords,
					    const vector<float> & a_coords, 
					    const msmat & templ, 
					    const msmat & align,
					    int linear_scans_offset );

  void adjust_ort_wrt_map_linear ( int nscans );
  void adjust_ort_wrt_map_linear ( float rt );



  

  virtual void linear_shift(int shift) = 0;
  virtual matrix_dim get_data_dim() const { 
    cerr << "called virtual get_data_dim()!" << endl;
    return matrix_dim(this->get_num_rts(), this->get_num_mzs());
  };
  virtual void set_scan( vector<float>& scan_data, int scan_idx ) = 0;
  virtual void set_scan ( vector<float> & scan_data, int scan_idx , LazyMatrix * d ) = 0;
  virtual void set_chrom ( vector<float> & chrom_data, int chrom_idx ) = 0;
  virtual void set_chrom ( vector<float> & chrom_data, int chrom_idx , LazyMatrix * d) = 0;

  virtual void get_scan(int scan_idx, vector<float> & inv) const = 0;
  virtual void get_chrom(int scan_idx, vector<float> & inv) const = 0;
  
  virtual vector<float> get_scan(int scan_idx) const = 0;
  virtual void get_scan_interp(float scan_coord, vector<float> & inv) const = 0;
  virtual vector<float> get_scan_interp(float scan_coord) const = 0 ;
  virtual vector<float> get_chrom(int chrom_idx) const = 0;
  //  virtual vector<float> * get_scan( int scan_idx ) = 0;
  //  virtual vector<float> * get_chrom( int scan_idx ) = 0;

  virtual size_t get_num_rows() const = 0;
  virtual size_t get_num_cols() const = 0;
  virtual void report_matrix_type() const = 0;
  virtual matrix_type_enum get_matrix_type() const = 0;

  //  void mask_based_on_ppids ( ppid_row_parser & p, float min_retain_i = 0.0f );
  //  void mask_based_on_hardklor ( hardklor_row_parser & p );
  void sample_values ( std::vector<float> & vals, const int num_vals, float min_threshold=0.0f, float max_threshold=0.0f );    
  void scan_scramble();

  


#ifdef HAVE_CRAWDAD

  
  double avg_diag_score ( msmat & other, SpectraCompFunc s = crawutils::spectra_cosine_angle  ) const ;
  double avg_diag_score ( msmat & other, const int start_scan_idx, const int stop_scan_idx, SpectraCompFunc s = crawutils::spectra_cosine_angle ) const;
  double total_diag_score( msmat & other, SpectraCompFunc s = crawutils::spectra_cosine_angle ) const ;
  double total_diag_score( msmat & other , const int start_scan_idx , const int stop_scan_idx, SpectraCompFunc = crawutils::spectra_cosine_angle ) const;

  void chrom_smooth( ChromSmoother & sg );
  void spike_filter ( ChromSmoother & spike_filter );

  void mask_based_on_ppids ( ppid_row_parser & p, float min_retain_i = 0.0f );
  void mask_based_on_hardklor ( hardklor_row_parser & p );

  void split_by_num_chroms ( int mz_idx_size, vector< msmat_scans * > & v );
  void split_by_num_segments ( int num_segments, vector<msmat_scans *> & v );


  ///weighted_diag_score : compares an msmat to another using the spectral comparison in SpectraCompareConfig,
  ///returns the total score and the total set of weights for the scores
std::pair<float,float> weighted_diag_score ( msmat & other , 
      const int start_scan_idx, const int stop_scan_idx, const SpectraCompareConfig & c ) const;
  std::pair<float,float> bp_weighted_diag_score ( msmat & other, const int start_scan_idx, const int stop_scan_idx, SpectraCompFunc comp_f = crawutils::spectra_cosine_angle ) const;

#endif


  	   


};



class msmat_scans: public msmat  {

public:
	matrix_dim get_data_dim() const;
  
  //  void set_scan ( float * scan_data, int scan_idx );
  //  void set_chrom ( float * chrom_data, int scan_idx);
  

    

  msmat * clone () {
     msmat * m = new msmat_scans();
     this->data->to_full_data();
     m->copy_from_msmat_deep(*this);
     return m;
  }
  void linear_shift(int shift);
  void set_scan  ( vector<float> & scan_data,  int scan_idx);
  inline void set_scan  ( vector<float> & scan_data,  int scan_idx, LazyMatrix * d);
  void set_chrom ( vector<float> & chrom_data, int mz_idx );
  inline void set_chrom ( vector<float> & chrom_data, int mz_idx , LazyMatrix * d);
  
  void get_chrom(int scan_idx, std::vector<float> & inv) const;
  void get_scan(int scan_idx, std::vector<float> & inv) const;
  void get_scan_interp(float scan_coord, vector<float> & inv) const;
  vector<float> get_scan_interp(float scan_coord) const ;

  std::vector<float> get_chrom( int chrom_idx ) const;
  std::vector<float> get_scan( int scan_idx) const;
  //vector<float> * get_scan(  int scan_idx );
  //vector<float> * get_chrom( int chrom_idx );
  
  virtual void zero_pixel ( int chrom_idx, int scan_idx );

	
	msmat_scans(FILE * msmat_file, sparse_level s = LM_FULL_DATA ) {
		init_from_msmat_file(msmat_file,s);
		//TODO make sure file and constructor match
		
	}
	msmat_scans() : msmat() {
	  int a = 3;
	};
	
    msmat * copy_no_data();
	
	size_t get_num_rows() const {
		size_t num_rts = get_num_rts();
		return num_rts;
	}
	size_t get_num_cols() const {
		return get_num_mzs();
	}
    bool has_summary_data() const {
        if ( bp_chrom.size() > 0 && tic_chrom.size() > 0 ) {
          return true;
        }
        else return false;
    }

	~msmat_scans() {
	  if ( data != NULL ) {
	    delete(data);
        data = NULL;
	  }
      if ( this->f_info != NULL ) {
        delete(f_info);
      }
#ifdef DEBUG
      std::cerr << "msmat_scans dtor" << std::endl;
#endif
	};
/*	
	msmat_scans( msmat_chroms & c ) {
			copy_from_msmat_nodata(c);
			//TODO check that this is a FULL MATRIX TYPE
			data = flip_lazy_matrix(*(c.data));
	}
*/

  

	void report_matrix_type() const {
	  cerr << "msmat_scans type!" << endl;
	};
	matrix_type_enum get_matrix_type() const {
	  return SCANS_MATRIX;
	};



private :
	void _trim_start( int nscans);
	void _trim_stop( int nscans);
  void _trim_start_mzs(int);
  void _trim_stop_mzs(int);
		
};

class msmat_chroms : public msmat {
public :
  msmat * clone ();
  matrix_dim get_data_dim() const;


  msmat * copy_no_data();

  void linear_shift(int shift);
	  
  void set_scan  ( vector<float> & scan_data,  int scan_idx);



  void set_scan  ( vector<float> & scan_data,  int scan_idx, LazyMatrix * d);
  void set_chrom ( vector<float> & chrom_data, int scan_idx );
  void set_chrom ( vector<float> & chrom_data, int scan_idx , LazyMatrix * d);

  void zero_pixel ( int chrom_idx, int scan_idx );
	  
  vector<float> get_scan(int scan_idx) const;
  vector<float> get_chrom(int chrom_idx) const;
  void get_scan_interp(float scan_coord, vector<float> & inv) const;
  vector<float> get_scan_interp(float scan_coord) const ;

        //vector<float> * get_scan(int scan_idx);
        //vector<float> * get_chrom(int chrom_idx);

  void get_scan(int scan_idx, vector<float> & inv) const;
  void get_chrom(int mz_idx, vector<float> & inv) const;

  vector<float> get_scan_vec( int scan_idx ) const;
  vector<float> get_chrom_vec( int chrom_idx ) const;
	


  msmat_chroms * extract_mzidx_range( int r_start, int r_stop );
  
  msmat_chroms * _extract_mzidx_range ( vector<int> mz_range ) {
    msmat_chroms * new_mat = new msmat_chroms();
    return new_mat;
    
  }
  
  msmat_chroms() : msmat() {
    int a  = 3;
  };
  
  
  msmat_chroms(FILE * msmat_file, sparse_level l = LM_FULL_DATA ) {
    init_from_msmat_file(msmat_file,l);
    //TODO make sure file and constructor match
    
  }

  

  size_t get_num_rows() const {
    size_t num_mzs = get_num_mzs();
    return num_mzs;
  }
  size_t get_num_cols() const {
    return get_num_rts();
  }
	
  ~msmat_chroms() { 
    if ( data != NULL ) {
      delete(data); 
      data = NULL;
    }
    if ( this->f_info != NULL ) { 
      delete(f_info);
    }
#ifdef DEBUG
    std::cerr << "msmat_chroms dtor" << std::endl;
#endif
  };
	
  void report_matrix_type()  const {
    cerr << "msmat_chroms type!" << endl;
  }
  matrix_type_enum get_matrix_type() const {
    return CHROMS_MATRIX;
  }


  ///scrambles the position of each scan in a run


	
private :
  void _trim_start(int nscans);
  void _trim_stop(int nscans);	
  void _trim_start_mzs(int);
  void _trim_stop_mzs(int);
	

	
};

/*  A justifaction for subclassing ttest_score_matrix from
    msmat_chroms, rather than ScoreMatrix -- ttest_score_matrix
    needs to be cognizant of m/z and retention time boundaries, which
    information is available to the msmat_chroms class. Will just
    need to be wary of any changes to msmat_chroms class which
    could complicate uses of the msmat_class.
*/

///A class for storing ttest result values from two groups. 
class ttest_score_matrix : public  msmat_chroms  {

public :

  ttest_score_matrix ( const msmat & m ) {
    copy_from_msmat_nodata(m);
    this->init_data();
    this->data->to_full_data(TTEST_NO_SEARCH);
    
    /* set all values in the matrix to the TTEST_NO_SEARCH value */
  };
  inline void unmask_search ( mz_type mz, rt_type start_rt, rt_type stop_rt);
  inline void unmask_search ( mz_type mz, int start_rt_idx, int stop_rt_idx);
  void unmask_search ( int mz_idx, int start_rt_idx, int stop_rt_idx);




  
  
};






#ifdef HAVE_CRAWDAD
class coord_opts {
private :
    crawusage::optset opts;
public:
    float start_rt,stop_rt,start_mz,stop_mz;
    coord_opts () { 
       start_rt = stop_rt = start_mz = stop_mz = -1.0f;
    }
    coord_opts( const crawusage::optset & o );
    void update_from_msmat ( const msmat & m );

};
#endif



int msmat_open(msmat ** msmat_f, crawutils::file_info & m_info, sparse_level sp = LM_FULL_DATA );
static int test_msmat_time_interp( msmat * msmat_f );
static int test_msmat_sg_test( msmat * msmat_fh);
static int test_msmat_output ( msmat * msmat_fh, char * out_fname);
static msmat * test_msmat(int argc, char ** argv);


static void output_bp_chrom( msmat & m , char * fname);
static msmat *  msmat_from_fh ( FILE * msmat_file, sparse_level s = LM_FULL_DATA );
static msmat * msmat_from_thevoid ( const std::vector<float> & mzs, const std::vector<float> & rts, std::vector < std::vector< float > > & chroms , bool as_scans = true );





#endif

