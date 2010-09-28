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
#ifndef _INCLUDE_MS1_MATRIX_BASE
#define _INCLUDE_MS1_MATRIX_BASE 0

#include <string>
#include <stdio.h>
#include <vector>
#include <assert.h>

#include "crawutils.H"

//define MACCOSS_LAB_LEGACY if you need to be able to try
//and read old .msmat files
#define MACCOSS_LAB_LEGACY 



#define HEADER_END "_END_MSMAT_HEADER_"
#define HEADER_END_LEN 18
#define HEADER_START "_MSMAT_V1_"
#define VALID_HEADER "_MSMAT_"
#define MAX_HEADER_LEN 20
#define LINE_LEN_LIM 256

#ifdef MACCOSS_LAB_LEGACY
#define HEADER_START_ARCHAIC "HEADER_START"
#define HEADER_START_ANCIENT "HEADER_START_NEW"
#endif




using namespace std;

int test_msmat_base();


typedef enum matrix_type {
  CHROMS_MATRIX,
  SCANS_MATRIX
} matrix_type_enum ;



int write_float_vect ( vector<float> & data, FILE * ofh );
int read_float_vect ( vector<float> & data, FILE * ofh , int num_fields );

/* returns a float indicating the version number 
0 - not valid
<0 - obsolete version */
float is_msmat( const char * fname);


const static int num_valid_fields = 12;
const static char * valid_fields[] = { 		  
  "array_type",		
  "rts",
  "mzs",
  "num_mz",
  "num_scans",
  "event_trail",
  "ort_wrt_map",
  "mzstr_id_map",
  "bin_size",
  "bp_chrom",
  "tic_chrom",
  "mzlabel_list",
};


#endif

