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
#ifndef _INCLUDE_MS1_MATRIX_HEADER_PARSER
#define _INCLUDE_MS1_MATRIX_HEADER_PARSER
/* would like to put into a separate file but had build problems */

#include "crawutils.H"
#include "msmat_base.H"
#include <stdio.h>
//#include "msmat.H"

class msmat;
class msmat_header_parser {
  friend class msmat;
 public:

  int read_header(msmat * m,FILE* header_fh, bool t, int * mat_type);
  int write_header(msmat * msmat, FILE* header_fh);
  //int valid_header();
  int valid_field(char * key_field);
  //void attach_msmat(msmat & in_msmat) {
  //    msmat = &(in_msmat);
  //}
  matrix_type_enum get_class_type ( FILE * header_fh );
  msmat_header_parser() {
      data_start_offset = -1;
  }
  off_t get_data_start() {
    return data_start_offset;
  }

  void load_map_from_str( flt_map & t_t_map , const char * field_data , int field_len );
  void load_strmap_from_str ( std::map< std::string , float > & str_map , const char * field_data , int field_len );
  void load_string_list_from_str ( vector<std::string> & in, const char * field_data, int field_len );

 private:

  int process_field(msmat * m,FILE * fname, char * value_start);
  void set_field(msmat * m,char * field_name, const char * field_data, int data_len);
  void write_header_field( FILE * header_fh, const char * field_name, 
			   size_t elements, size_t element_size, 
			   const void * element_value );
  void write_header_field_float_map ( FILE * header_fh, const char * field_name,
				      flt_map & float_map );
  void write_header_field_string_list ( FILE * header_fh, const char * field_name,
				  const vector<std::string> & s );
  void write_header_field_strmap( FILE * header_fh, const char * field_name, 
            const std::map< std::string, float > & str_map );


  //  msmat * msmat;
  off_t data_start_offset;

};
#endif
