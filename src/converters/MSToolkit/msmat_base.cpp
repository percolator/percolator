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
#include "msmat_base.H"
#include <vector>
#include <map>
#include <string>
#include <stdio.h>

using namespace std;

int write_float_vect ( vector<float> & data, FILE * ofh ) {
  size_t cnt;
  cnt = fwrite((void*)&data[0], sizeof(float), data.size(), ofh);
  if ( cnt != data.size() ) {
    return 0;
    //TODO -- raise exception
  }
  return 1;
}

int read_float_vect ( vector<float> & data, FILE * ofh , int num_fields ) {
  size_t cnt;
  data.resize(num_fields);
  cnt = fread((void*)&data[0], sizeof(float), data.size(), ofh);
  if ( cnt != data.size() ) {
    return 0;
    //TODO -- raise exception
  }
  return 1;
}

float is_msmat( const char * fname) {
  
	
  //float version  = 0.0 ;
  float status = 0;
  char line[LINE_LEN_LIM];
  int linelen;
  FILE * ms_fh = fopen(fname,"rb");
  if (fgets(line,LINE_LEN_LIM,ms_fh) == NULL ) {
    //TODO raise error
  }
  linelen = strlen(line);
  line[linelen-1] = '\0';
  //make sure we don't have \r\n characters

  if ( line[linelen-2] == '\r' ) {
	  line[linelen-2] = '\0';
  }

  if ( strncmp(line,VALID_HEADER,MAX_HEADER_LEN) == 0 ) {
    //TODO -- properly deal with version number
    status = 1.0;
  #ifdef MACCOSS_LAB_LEGACY
  } else if ( strcmp(line, HEADER_START_ARCHAIC) == 0 ) {
    status = -1.0;
  } else if ( strcmp(line, HEADER_START_ANCIENT) == 0 ) {
    status = -1.0;
  #endif
  }
  
  return status;
  
}



