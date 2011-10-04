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
#ifdef CRAW_LOGGING
#include "logging.H"
#endif
#include "msmat_header_parser.H"
#include "msmat.H"

using namespace std;

class msmat;



int msmat_header_parser::read_header(msmat * m, FILE * header_fh, bool get_class, int * class_type) {
  /* note that we need to avoid using line-reading  functions since
  * some  of the binary data may contain \n or \r\n
  */
  /* Game plan
  * 1. check the header field
  * 2. go through each field, check that it contains one of the valid string
  * 3. */

  int status = 1;
  fseek(header_fh,0,SEEK_SET);
  char cur_line[128];
  char key_field[128];

  fgets(cur_line,128,header_fh);
  int last_line_start = ftell(header_fh);
  /* should be at next line */;

  while ( 1 ) {
    fgets(cur_line,128,header_fh);
    if ( strncmp(cur_line,HEADER_END,HEADER_END_LEN) == 0) {
      this->data_start_offset = last_line_start + HEADER_END_LEN + 1 ;
      return status;
    }
    strncpy( key_field, cur_line, 128); 
    char * colon_pos;
    if ( ( colon_pos = strchr(key_field,':') ) == NULL ) {
      cerr << "Failed to find semicolon in header field" << key_field << endl; 
      exit(-1);
    }
    if ( get_class ) {
      if ( class_type != NULL && ( strncmp(key_field,"array_type",10) == 0  ) ) {
        /* read the digit, then the number of bytes , comma etc.. */
        char * array_field = colon_pos + 1;
        /* TODO -- fix this horrrible hack */
        if ( strncmp(array_field,"5,scans",7) == 0 ) {
          *class_type = (int)SCANS_MATRIX;
        }
        else if ( strncmp(array_field,"6,chroms",8) == 0 ) {
          *class_type = (int)CHROMS_MATRIX;
        }
        else {
          cerr << "Invalid array_type " << array_field << endl;
          exit(-1);
        }
      }
      return 1;
    }
    else {
      *colon_pos = '\0'; 
      if ( valid_field(key_field) ) {


        /* pointer to value start and key field */

        int semi_start_pos = last_line_start + (int)( colon_pos - key_field );
        fseek(header_fh,semi_start_pos,SEEK_SET);
        char should_be_semi = (char)fgetc(header_fh);
        assert(should_be_semi == ':');		
        if ( (last_line_start = process_field(m,header_fh, key_field )) < 0 ) {
          cerr << "error processing field " << key_field << endl;
          exit(-1);
        }

      }
    }

  }
  return 1;
}



int msmat_header_parser::write_header(msmat * msmat, FILE  * header_fh ) {
  int status = 1;

  /*
  Good link for writing binary data (packer / unpacker)
  http://www.codeproject.com/file/packersfx.asp

  example for writing std::vector w/ fwrite
  fwrite(&filesList[0], sizeof(pdata), filesList.size(), fpArchive);

  1. write the header elements
  */

  vector<float> tmp_rts = msmat->get_rts();
  vector<float> tmp_mzs = msmat->get_mzs();

  const char * array_type_name;
  array_type_name = ( msmat->get_matrix_type() == SCANS_MATRIX ? "scans" : "chroms");

  fprintf(header_fh,"%s\n",HEADER_START);
  write_header_field( header_fh, "array_type", strlen(array_type_name), sizeof(char), array_type_name);
  write_header_field( header_fh, "mzs", msmat->get_num_mzs(), sizeof(float), &(msmat->mzs[0]) );
  write_header_field( header_fh, "rts", msmat->get_num_rts(), sizeof(float), &(msmat->rts[0]) );
  if ( msmat->summary_data_ok ) {
    write_header_field( header_fh, "bp_chrom",msmat->get_num_rts(), sizeof(float), &(msmat->bp_chrom[0]));
    write_header_field( header_fh, "tic_chrom",msmat->get_num_rts(), sizeof(float), &(msmat->tic_chrom[0]));
  }

  write_header_field( header_fh, "bin_size", 1, sizeof(float), 
    &(msmat->bin_size));
  write_header_field_float_map( header_fh, "ort_wrt_map", msmat->ort_wrt_map);
  write_header_field_string_list( header_fh, "mzlabel_list", msmat->labels);

  

  fprintf(header_fh,"%s\n",HEADER_END);
  return status;
};



int msmat_header_parser::valid_field( char * field ) {
  for ( int i = 0 ; i < num_valid_fields ; i++ ) {
    if  ( strcmp(field, valid_fields[i]) == 0 ) {
      return 1;
    }
  }
  cerr << "invalid field:" << field << endl;
  return 0;
}

matrix_type_enum msmat_header_parser::get_class_type( FILE * header_fh ) {
  int m_type_tmp;
  (matrix_type_enum)read_header(NULL,header_fh, true, &m_type_tmp );
  matrix_type_enum ret_val = (matrix_type_enum)m_type_tmp;
  return ret_val;
}
/* end public */

/* process_field -- given a field we use a lookup field given the key
* first value should be the length of the field in bytes as decimals , then a comma, then \n
* 
*/
int msmat_header_parser::process_field(msmat * m, FILE * header_fh, char * field_name) {
  /* get the field length */
  int field_len = -1;

  /* fscanf returns number of fields parsed */
  if ( fscanf(header_fh,"%d,",&field_len) != 1 ) {
    char line_buf[256];
    cerr << "error trying to find data from line -- failed";
    fgets(line_buf,256,header_fh);
    cerr << "line:" << line_buf << endl;
    exit(-1);
  }

  if ( field_len < 0 ) {
    cerr << "field_len of " << field_len << " processing " << field_name << endl;
    exit(-1);
  }
  else if ( field_len == 0 ) {
    cerr << "field_len of " << field_len << " processing " << field_name << endl;
  }
  else {
    char * field_data_input = (char*)calloc(field_len, sizeof(char) );
    fread(field_data_input,sizeof(char),field_len,header_fh);
    const char * field_data = (const char *)field_data_input;
    set_field(m,field_name,field_data, field_len);
    free(field_data_input);
  }
  char should_be_newline = (char)fgetc(header_fh);
  if ( should_be_newline != '\n') {
    char tmp_data[65];
    size_t fpos = ftell(header_fh);
    fseek(header_fh,fpos-32,SEEK_SET);
    fread(tmp_data,sizeof(char),64,header_fh);
    fprintf(stdout,"Nearby data: %s\n",tmp_data);
    tmp_data[34] = 0;
    fprintf(stdout,"Flanking: %s\n",tmp_data + 31);
    exit(-1);
  }
  return (ftell(header_fh));
}

void msmat_header_parser::set_field(msmat * msmat, char * field_name, const char * field_data, int data_len) {
  /* lookup table for what to do with data based on field_name */
  /* just do a monstrous statement for now */
  if ( strcmp(field_name,"array_type") == 0) {
    if ( strcmp(field_name,"chroms_matrix") == 0 ) {
      //msmat->matrix_type = CHROMS_MATRIX;
    }
    else if ( strcmp(field_name,"scans_matrix") == 0 ) {	
      //msmat->matrix_type = SCANS_MATRIX;
    }
  }
  else if ( strcmp(field_name, "rts") == 0 ) {
    vector<float> tmp_rts(data_len / sizeof(float) );
    load_array<float>(&tmp_rts[0],field_data,data_len);
    msmat->set_rts(tmp_rts);
  }
  else if ( strcmp( field_name, "mzs") == 0) {	
    vector<float> tmp_mzs(data_len / sizeof(float));
    load_array<float>(&tmp_mzs[0],field_data,data_len);
    msmat->set_mzs(tmp_mzs);

  }
  else if ( strcmp( field_name, "bin_size") == 0) {
    load_type<float>(&(msmat->bin_size),field_data,data_len);
  }
  else if ( strcmp( field_name, "num_mz") == 0 ) {
    //load_int_type<int>(&(msmat->num_mzs),field_data,data_len);
  }
  else if ( strcmp ( field_name, "num_scans") == 0 ) {
    //load_int_type<int>(&(msmat->num_rts),field_data,data_len);
  }
  else if ( strcmp ( field_name, "event_trail") == 0 ) {
    //load_not_implemented()
  }
  else if ( strcmp (field_name, "bp_chrom" ) == 0 ) {
    vector<float> bp_chrom(data_len / sizeof(float));
    load_array<float>(&bp_chrom[0],field_data,data_len);
    msmat->bp_chrom.resize(bp_chrom.size());
    copy(bp_chrom.begin(), bp_chrom.end(), msmat->bp_chrom.begin());
  }
  else if ( strcmp (field_name, "tic_chrom" ) == 0 ) {
    vector<float> tic_chrom(data_len / sizeof(float));
    load_array<float>(&tic_chrom[0],field_data,data_len);
    msmat->tic_chrom.resize(tic_chrom.size());
    copy(tic_chrom.begin(), tic_chrom.end(), msmat->tic_chrom.begin());
  }

  else if ( strcmp ( field_name, "ort_wrt_map") == 0) {
    load_map_from_str(msmat->ort_wrt_map, field_data, data_len);
  }
  //else if ( strcmp ( field_name, "mzstr_id_map") == 0) {
  //  load_strmap_from_str(msmat->mzstr_id_map, field_data, data_len);
  //}
  else if ( strcmp ( field_name, "mzlabel_list") == 0 ) {
    load_string_list_from_str( msmat->labels, field_data, data_len);
    msmat->init_transitions_from_labels();
    msmat->init_transition_manager();
  }
  else {
    cerr   << "did not recognize field name: " << field_name << endl;
  }
}

void msmat_header_parser::write_header_field ( FILE * header_fh, const char * field_name,
                                                   size_t elements, size_t element_size, const void * element_value ) {
                                                     fprintf(header_fh,"%s:%zd,",field_name, elements * element_size);
                                                     if ( elements > 0 ) {
                                                       fwrite(element_value, element_size, elements, header_fh);
                                                     }
                                                     fprintf(header_fh,"\n");
}

void msmat_header_parser::load_map_from_str( flt_map & t_t_map , const char * field_data , int field_len ) {
  //check that the size of the field makes sense
  if ( field_len % ( 2 * sizeof(flt_map_type ) ) != 0 ) {
#ifdef CRAW_LOGGING
    LOGH(LOG_V1,std::cerr) << "invalid length of flt_map field" << std::endl; exit(1);
#endif
  }
  uint num_items = (field_len / 2) / sizeof(flt_map_type);
  uint items_length = field_len / 2;

  flt_map_type * k = new flt_map_type[num_items];
  flt_map_type * v = new flt_map_type[num_items];

  load_array( k , field_data , items_length);
  load_array( v , field_data + items_length , items_length);
  for ( uint i = 0 ; i < num_items ; i++ ) {
    t_t_map.insert(pair<flt_map_type, flt_map_type>( k[i], v[i] ));
  }
  delete k;
  delete v;
};


void msmat_header_parser::load_strmap_from_str( std::map< std::string , float > & str_map , const char * field_data , int field_len ) {


 /*	    
	    c1 = header_data.index(chr(1))
	    key_len = int(header_data[:c1])
	    key_last = c1 + key_len
	    key_data = header_data[c1+1:key_last]
	    assert(header_data[key_last] == chr(0))
	    keys = key_data.split(chr(0))
	    
	    c2 = header_data.index(chr(1),key_last)
	    val_len = int(header_data[key_last+1:c2])
	    assert(c2+1+val_len == len(header_data))
	    val_bytes = header_data[c2+1:c2+1+val_len]
	    val_array = Numeric.fromstring(val_bytes,typecode='f')
	    assert(len(val_array) == len(keys))
	    m = {}
	    for k_idx in range(len(keys)) :
		    m[keys[k_idx]] = val_array[k_idx]
	    return m
    
*/

  int c1_loc, c2_loc,key_fields_len, key_last_loc, key_data_len, key_len, val_len;
  char key_len_data[128], val_len_data[128];
  std::string field_data_str(field_data);
  c1_loc = -1;
  
  c1_loc = field_data_str.find(char(1));

  if ( c1_loc == -1 ) {
     #ifdef CRAW_LOGGING
    LOGH(LOG_V1,std::cerr) << "invalid length of flt_map field" << std::endl; exit(1);
    #endif
    exit(1);
  }
  
  std::string key_len_str(field_data_str,0,c1_loc);
  key_len = atoi(key_len_str.c_str());
  assert(key_len > 0);
  int key_last = c1_loc + key_len;
  assert(field_data[key_last] == char(0));
  std::string key_data(field_data_str,c1_loc+1,key_len);
  std::vector< std::string> key_toks;
  
  //go through key_data, splitting on nulls, into strings
  
  int key_data_idx;
  std::string sbuf;
  for ( int i = c1_loc+1 ; i < key_last+1 ; i++ ) {
      if ( field_data[i] == '\0' ) {
        key_toks.push_back(sbuf);
        sbuf.clear();
      }
      else {
        sbuf.push_back(field_data[i]);
      }
  }
  
  
  c2_loc = field_data_str.find(char(1),key_last);
  std::string val_len_str(field_data_str,key_last+1,(c2_loc-(key_last+1)));
  val_len = atoi(val_len_str.c_str());


  assert(c2_loc+1+val_len == field_data_str.size());
  
  std::string value_data(field_data_str,c2_loc+1,val_len);
  const char * value_data_chr = value_data.c_str();
  float * values = new float[val_len];
  load_array( values , value_data_chr , val_len);
  for ( uint i = 0 ; i < key_toks.size() ; i++ ) {
    str_map.insert(pair< std::string, float  >( key_toks[i],values[i]) );
  }
  delete values;
}
  
  
  
 



void msmat_header_parser::write_header_field_float_map ( FILE * header_fh, const char * field_name,
                                                             flt_map & float_map ) {
                                                               int num_elements = float_map.size();
                                                               size_t map_array_len = num_elements * 2 * sizeof(float);
                                                               fprintf(header_fh,"%s:%zd,",field_name,map_array_len);
                                                               if  ( num_elements > 0 ) {
                                                                 vector<float> map_keys(num_elements);
                                                                 vector<float> map_vals(num_elements);
                                                                 int ele_count = 0;
                                                                 flt_map::const_iterator i = float_map.begin();
                                                                 for ( ; i != float_map.end() ; i++, ele_count++ ) {
                                                                   float key = i->first;
                                                                   float val = i->second;
                                                                   map_keys[ele_count] = key;
                                                                   map_vals[ele_count] = val;

                                                                 }
                                                                 assert(ele_count == num_elements);
                                                                 fwrite( (const void*)&(map_keys[0]) , sizeof(float), num_elements, header_fh);
                                                                 fwrite( (const void*)&(map_vals[0]) , sizeof(float), num_elements, header_fh);
                                                               }
                                                               fprintf(header_fh,"\n");
}



  void msmat_header_parser::write_header_field_string_list ( FILE * header_fh, const char * field_name,
							    const vector<std::string> & s ) {
    int num_elements = s.size();
    uint total_len = 0;
    for ( uint i = 0 ; i < s.size() ; i++ ) {
      total_len += s[i].length();
    }
    total_len += s.size(); // for null characters
    fprintf(header_fh,"%s:%d,",field_name,total_len);
    for ( uint i = 0 ; i < s.size() ; i++ ) {
      fwrite(s[i].c_str(), sizeof(char), s[i].size() + 1, header_fh); 
    }
    fprintf(header_fh,"\n");
  }

void msmat_header_parser::load_string_list_from_str ( vector<std::string> & in, const char * field_data, int field_len ) {
  crawutils::read_string_array( in, field_data, field_len );
}

/* End msmat_header_parser */

