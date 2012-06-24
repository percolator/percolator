/*******************************************************************************
 Copyright 2006-2009 Lukas Käll <lukas.kall@cbr.su.se>

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

 *******************************************************************************/

#ifndef READER_H
#define READER_H

#include "FragSpectrumScanDatabase.h"
#include <string>
#include <set>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <vector>
#include "Globals.h"
#include "Enzyme.h"
#include "MassHandler.h"
#include "DataSet.h"
#include "FeatureNames.h"
#include "percolator_in.hxx"
#include "parseoptions.h"
#include <boost/filesystem.hpp>
#include <assert.h>
#if defined (__WIN32__) || defined (__MINGW__) 
#include <direct.h>
#include <io.h>
#include <stdio.h>
#define  mkdir( D, M )   _mkdir( D )
#include <fcntl.h>
#include <errno.h>
#endif

using namespace std;

class Reader
{

public:
  
  Reader(ParseOptions po);
  virtual ~Reader();
  
  //NOTE template function need to be declared and implemented in header
  template<class T>
  void translateFileToXML(const std::string fn,
     ::percolatorInNs::featureDescriptions & fds,
     ::percolatorInNs::experiment::fragSpectrumScan_sequence & fsss, bool isDecoy,
     vector<T*>& databases, unsigned int lineNumber_par)
  {
  
    //TODO check its is a metafile or not and if it is valid formed
    if(checkValidity(fn))
    {
      // there must be as many databases as lines in the metafile containing sqt
      // files. If this is not the case, add a new one
      if(databases.size()==lineNumber_par)
      {
	// initialize databese     
	T* database = new T(fn);
      
	//NOTE this is actually not needed in case we compile with the boost-serialization scheme
	//indicate this with a flag and avoid the creating of temp files when using boost-serialization
	if(!po.boost_serialization)
	{
	  // create temporary directory to store the pointer to the database
	  string tcf = "";
	  char * tcd;
	  string str;
      
	  //TODO it would be nice to somehow avoid these declararions and therefore avoid the linking to
	  //filesystem when we dont use them
	  try
	  {
	    boost::filesystem::path ph = boost::filesystem::unique_path();
	    boost::filesystem::path dir = boost::filesystem::temp_directory_path() / ph;
	    boost::filesystem::path file("converters-tmp.tcb");
	    tcf = std::string((dir / file).string()); 
	    str =  dir.string();
	    tcd = new char[str.size() + 1];
	    std::copy(str.begin(), str.end(), tcd);
	    tcd[str.size()] = '\0';
	    if(boost::filesystem::is_directory(dir))
	    {
	      boost::filesystem::remove_all(dir);
	    }
	
	    boost::filesystem::create_directory(dir);
	  } 
	  catch (boost::filesystem::filesystem_error &e)
	  {
	    std::cerr << e.what() << std::endl;
	  }
	
	  tmpDirs.resize(lineNumber_par+1);
	  tmpDirs[lineNumber_par]=tcd;
	  tmpFNs.resize(lineNumber_par+1);
	  tmpFNs[lineNumber_par]=tcf;
	  database->init(tmpFNs[lineNumber_par]);
	}
	else
	{
	  database->init("");
	}
      
	databases.resize(lineNumber_par+1);
	databases[lineNumber_par]=database;
	assert(databases.size()==lineNumber_par+1);
      }
      if (VERB>1){
	std::cerr << "reading " << fn << std::endl;
      }
    
      getMaxMinCharge(fn);
      read(fn,fds, fsss, isDecoy,databases[lineNumber_par]);
    
    } else {
      // we hopefully found a meta file
      unsigned int lineNumber=0;
      std::string line2;
      std::ifstream meta(fn.data(), std::ios::in);
      while (getline(meta, line2)) {
	if (line2.size() > 0 && line2[0] != '#') {
	  //NOTE remove the whitespaces
	  line2.erase(std::remove(line2.begin(),line2.end(),' '),line2.end());
	  translateFileToXML(line2, fds, fsss, isDecoy,databases, lineNumber);
	  lineNumber++;
	}
      }
      meta.close();
    }
  }

  string getRidOfUnprintables(std::string inpString);
  
  virtual void read(const std::string fn,
    ::percolatorInNs::featureDescriptions& fds,
     ::percolatorInNs::experiment::fragSpectrumScan_sequence& fsss,
      bool is_decoy,FragSpectrumScanDatabase* database){};
      
  virtual bool checkValidity(std::string file){};
  
  virtual void getMaxMinCharge(std::string fn){};
  
  void push_backFeatureDescription(
    percolatorInNs::featureDescriptions::featureDescription_sequence&,
    const char*);

  void computeAAFrequencies(const string& pep,
    percolatorInNs::features::feature_sequence & f_seq);
  
  double calculatePepMAss(std::string pepsequence,double charge);

  void initMassMap(bool useAvgMass);
  
private:
  
   std::vector<char*> tmpDirs;
   std::vector<std::string> tmpFNs;

protected:
  
   static const std::string aaAlphabet;
   static const std::string ambiguousAA;
   static const std::string modifiedAA;
   int maxCharge;
   int minCharge;
   ParseOptions po;
   std::map<char, double> massMap_;
};

#endif



