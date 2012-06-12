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

 
/*****************************************************************************
  
    Implementation of method to estimate protein FDR as described in :
    
    http://prottools.ethz.ch/muellelu/web/LukasReiter/Mayu/
 
    Mayu
    Lukas Reiter
    Manfred Claassen

    Mayu Software
    Lukas Reiter - Hengartner Laboratory
    lukas.reiter@molbio.uzh.ch
    Institute of Molecular Biology
    Winterthurerstrasse 190
    University of Zürich - Irchel
    CH-8057 Zürich
    +++++++++++++++++++++++++++++++++++++++++
    Located at:
    Institute for Molecular Systems Biology
    Aebersold Laboratory
    Wolfgang-Pauli-Str. 16
    ETH Hönggerberg, HPT C 75
    CH-8093 Zürich
    Tel: +41 44 633 39 45

****************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <iostream>
#include <algorithm>
#include <ProteinFDRestimator.h>
#include <math.h>
#include <cmath>

struct RetrieveValue
{
  template <typename T>
  typename T::second_type operator()(T keyValuePair) const
  {
    return keyValuePair.second;
  }
};

/** external functions used to estimate the expexted value of the hypergeometric distribution **/
double stirling_log_factorial(double n)
{
  double PI = 3.141592653589793;
  return ( log(sqrt(2*PI*n)) + n*log(n) - n);
}

double exact_log_factorial(double n)
{
  double log_fact = 0;
  for(int i = 2; i <= n; i++)
    log_fact += log(i);
  return log_fact;
}

double log_factorial(double n)
{
  if(n < 1000)
    return exact_log_factorial(n);
  else
    return stirling_log_factorial(n);
}

double log_binomial(double n,double k)
{
  return (log_factorial(n) - log_factorial(k) - log_factorial(n-k));
}


double hypergeometric(int x,int N,int w,int d)
{
  //natural logarithm of the probability
  if(d > 0)return exp(log_binomial(w,x) + log_binomial(N-w,d-x) - log_binomial(N,d));
  else return 0.0;
}

/* Simple API for FASTA file reading
 * for Bio5495/BME537 Computational Molecular Biology
 * SRE, Sun Sep  8 05:35:11 2002 [AA2721, transatlantic]
 * CVS $Id$
 */

FASTAFILE *
OpenFASTA(const char *seqfile)
{
  FASTAFILE *ffp;

  ffp = (FASTAFILE *)malloc(sizeof(FASTAFILE));
  ffp->fp = fopen(seqfile, "r");              /* Assume seqfile exists & readable!   */
  if (ffp->fp == NULL) 
  { 
    free(ffp); 
    return NULL; 
  } 
  if ((fgets(ffp->buffer, FASTA_MAXLINE, ffp->fp)) == NULL)
    { free(ffp); return NULL; }
  return ffp;
}

int
ReadFASTA(FASTAFILE *ffp, char **ret_seq, char **ret_name, /*char **ret_gene,*/ int *ret_L)
{
  char *s;
  char *name;
  char *seq;
  
  //TODO it would be nice to get the gene name but not all the databases contain it and it is located in 
  //different positions
  //char *gene;
  
  int   n;
  int   nalloc;
  
  /* Peek at the lookahead buffer; see if it appears to be a valid FASTA descline.
   */
  if (ffp->buffer[0] != '>') return 0;    

  /* Parse out the name: the first non-whitespace token after the >
   */
  s  = strtok(ffp->buffer+1, " \t\n");
  name = (char *)malloc(sizeof(char) * (strlen(s)+1));
  strcpy(name, s);

  /* Everything else 'til the next descline is the sequence.
   * Note the idiom for dynamic reallocation of seq as we
   * read more characters, so we don't have to assume a maximum
   * sequence length.
   */
  seq = (char *)malloc(sizeof(char) * 128);     /* allocate seq in blocks of 128 residues */
  nalloc = 128;
  n = 0;
  while (fgets(ffp->buffer, FASTA_MAXLINE, ffp->fp))
    {
      if (ffp->buffer[0] == '>') break;	/* a-ha, we've reached the next descline */

      for (s = ffp->buffer; *s != '\0'; s++)
	{
	  if (! isalpha(*s)) continue;  /* accept any alphabetic character */

	  seq[n] = *s;                  /* store the character, bump length n */
	  n++;
	  if (nalloc == n)	        /* are we out of room in seq? if so, expand */
	    {			        /* (remember, need space for the final '\0')*/
	      nalloc += 128;
	      seq = (char *)realloc(seq, sizeof(char) * nalloc);
	    }
	}
    }
  seq[n] = '\0';

  *ret_name = name;
  *ret_seq  = seq;
  *ret_L    = n;
  /**ret_gene = gene;*/
  
  return 1;
}      

void
CloseFASTA(FASTAFILE *ffp)
{
  fclose(ffp->fp);
  free(ffp);
}

/******************************************************************************************************************/

ProteinFDRestimator::ProteinFDRestimator(unsigned int __minpeplength, unsigned int __minmaxx,
				       unsigned int __maxmass,std::string __decoy_prefix, double __missed_cleavages, 
				       unsigned __nbins, double __targetDecoyRatio, bool __binequalDeepth, unsigned __maxSeqlength)
				       :minpeplength(__minpeplength),minmass(__minmaxx),maxmass(__maxmass),
				       decoy_prefix(__decoy_prefix),missed_cleavages(__missed_cleavages),
				       nbins(__nbins),targetDecoyRatio(__targetDecoyRatio),
				       binequalDeepth(__binequalDeepth),maxSeqlength(__maxSeqlength)
{
  initMassMap();
}

ProteinFDRestimator::~ProteinFDRestimator()
{
  FreeAll(binnedProteins);
  FreeAll(groupedProteins);
  FreeAll(lenghts);
}



void ProteinFDRestimator::parseDataBase(const char* seqfile,const char* seqfileDecoy)
{
  FASTAFILE *ffp;
  FASTAFILE *ffpDecoy;
  char *seq;
  char *seqDecoy;
  char *name;
  char *nameDecoy;
  int   L;
  int   Ldecoy;
  
  std::map<std::string,std::string> targetProteins;
  std::map<std::string,std::string> decoyProteins;
  
  //NOTE I do not parse the gene names because there is not standard definiton for it
  try
  {
    ffp = OpenFASTA(seqfile);
    if(ffp != NULL)
    {
      while (ReadFASTA(ffp, &seq, &name, &L))
	{
	  targetProteins.insert(std::make_pair<std::string,std::string>(name,seq));
	  free(seq);
	  free(name);
	}
    
      CloseFASTA(ffp);
    }
    else
    {
      std::cerr <<  "Error reading target database : " << seqfile <<  std::endl;
      exit(-1);
    }
    
    ffpDecoy = OpenFASTA(seqfileDecoy);
    
    if(ffpDecoy != NULL)
    {
      while(ReadFASTA(ffpDecoy, &seqDecoy, &nameDecoy, &Ldecoy))
	{
	  decoyProteins.insert(std::make_pair<std::string,std::string>(nameDecoy,seqDecoy));
	  free(seqDecoy);
	  free(nameDecoy);
	}
      
      CloseFASTA(ffpDecoy);
    }
    else
    {
      std::cerr <<  "Error reading decoy database : " << seqfileDecoy <<  std::endl;
      exit(-1);
    }
  }
  catch(const std::exception &e)
  {
    e.what();
  }
  
  if(decoyProteins.size() == 0 || targetProteins.size() == 0)
  {
    std::cerr <<  "Error parsing the databases, one of the databases given is empty\n" << std::endl;
    exit(-1);
  }
  else if(VERB > 2)
  {  
    std::cerr << "\nRead " << targetProteins.size() << " target proteins and " << decoyProteins.size() << " decoys proteins\n" << std::endl;
  }
  
  if(targetProteins.size() != decoyProteins.size())
  {
    targetDecoyRatio = (double)targetProteins.size() / (double)decoyProteins.size();
  }
  else
  {  
    targetDecoyRatio = 1.0;
  }
  
  correctIdenticalSequences(targetProteins,decoyProteins);
}

void ProteinFDRestimator::parseDataBase(const char* seqfile)
{
  FASTAFILE *ffp;
  char *seq;
  char *name;
  int   L;
  
  std::map<std::string,std::string> targetProteins;
  std::map<std::string,std::string> decoyProteins;
  
  //NOTE I do not parse the gene names because there is not standard definiton for it
  try
  {
    ffp = OpenFASTA(seqfile);
    
    if(ffp != NULL)
    {
      while (ReadFASTA(ffp, &seq, &name, &L))
	{
	  
	  std::string name2(name);
	  if(name2.find(decoy_prefix) != std::string::npos)
	    decoyProteins.insert(std::make_pair<std::string,std::string>(name2,seq));
	  else
	    targetProteins.insert(std::make_pair<std::string,std::string>(name2,seq));
      
	  free(seq);
	  free(name);
	}
	
      CloseFASTA(ffp);
    }
    else
    {
      std::cerr <<  "Error reading combined database : " << seqfile <<  std::endl;
      exit(-1);
    }
  
    
  }catch(const std::exception &e)
  {
    e.what();
  }
  
  if(decoyProteins.size() == 0 || targetProteins.size() == 0)
  {
    std::cerr << "Error parsing the database, the database given does not contain either decoys or targets\n";
    exit(-1);
  }
  else if(VERB > 2)
  {  
    std::cerr << "\nRead " << targetProteins.size() << " target proteins and " << decoyProteins.size() << " decoys proteins\n" << std::endl;
  }
  
  if(targetProteins.size() != decoyProteins.size())
  {
    targetDecoyRatio = (double)targetProteins.size() / (double)decoyProteins.size();
  }
  else
  {  
    targetDecoyRatio = 1.0;
  }
  
  correctIdenticalSequences(targetProteins,decoyProteins);
}

void ProteinFDRestimator::correctIdenticalSequences(const std::map<std::string,std::string> &targetProteins,
						    const std::map<std::string,std::string> &decoyProteins)
{
  std::map<std::string,std::string>::const_iterator it,it2;
  
  groupedProteins.clear();
  lenghts.clear();
  it = targetProteins.begin();
  it2 = decoyProteins.begin();
  
  std::set<std::string> previouSeqs;
  double length = 0.0;
  
  for(it = targetProteins.begin(); it != targetProteins.end(); it++)
  {
    std::string targetSeq = (*it).second;
    std::string targetName = (*it).first;
    if(previouSeqs.find(targetSeq) != previouSeqs.end())
    {
      length = 0.0;
    }
    else
    {
      length = calculateProtLength(targetSeq);
      previouSeqs.insert(targetSeq);
    }
    groupedProteins.insert(std::make_pair<double,std::string>(length,targetName));
    lenghts.push_back(length);
  }
  
  for(it2 = decoyProteins.begin(); it2 != decoyProteins.end(); it2++)
  {
    std::string decoySeq = (*it2).second;
    std::string decoyName = (*it2).first;
    if(previouSeqs.find(decoySeq) != previouSeqs.end())
    {
      length = 0.0;
    }
    else
    {
      length = calculateProtLength(decoySeq);
      previouSeqs.insert(decoySeq);
    }
    groupedProteins.insert(std::make_pair<double,std::string>(length,decoyName));
    lenghts.push_back(length);
  }
  
  //FreeAll(previouSeqs);
  return;
}

void ProteinFDRestimator::groupProteinsGene()
{
  /**not implemented yet**/
}

void ProteinFDRestimator::estimateFDRthread(unsigned i,const std::set<std::string> &target, const std::set<std::string> &decoy) 
{ 
  unsigned  numberTP = countProteins(i,target);
  unsigned  numberFP = countProteins(i,decoy);
  unsigned  N = getBinProteins(i);
  double fp = estimatePi0HG(N,numberTP,targetDecoyRatio*numberFP);
  if(VERB > 2)
      std::cerr << "\nEstimating FDR for bin " << i << " in thread " << boost::this_thread::get_id() << " with " << numberFP << " Decoy proteins, " 
	  << numberTP << " Target proteins, and " << N << " Total Proteins in the bin " << " with exp fp " << fp << std::endl;
  fptol += fp;
} 

double ProteinFDRestimator::estimateFDR(const std::set<std::string> &__target, const std::set<std::string> &__decoy)
{   
  
    time_t startTime;
    clock_t startClock;
    time(&startTime);
    startClock = clock();
    
    if(binnedProteins.size() > 0)
    {
       FreeAll(binnedProteins);
    } 
    if(binequalDeepth)
    {
      binProteinsEqualDeepth();
    }
    else
    {
      binProteinsEqualWidth();
    }
    
    if(VERB > 2)
      std::cerr << "\nThere are : " << __target.size() << " target proteins and " << __decoy.size() << " decoys proteins that contains high confident PSMs\n" << std::endl;    

    fptol = 0.0;
    boost::thread t[nbins]; 
    
    for(unsigned i = 0; i < nbins; i++)
    {
       t[i] = boost::thread(boost::bind(&ProteinFDRestimator::estimateFDRthread,this,i,__target,__decoy)); 
       t[i].join();
    }
    
    time_t procStart;
    clock_t procStartClock = clock();
    time(&procStart);
    double diff = difftime(procStart, startTime);
    if (VERB > 2) std::cerr << "\nEstimating the protein FDR took : "
      << ((double)(procStartClock - startClock)) / (double)CLOCKS_PER_SEC
      << " cpu seconds or " << diff << " seconds wall time\n" << std::endl;
    
    if(isnan(fptol) || isinf(fptol) || fptol == 0)fptol = -1;
    return fptol ;
}



void ProteinFDRestimator::binProteinsEqualDeepth()
{
  //assuming lengths sorted from less to bigger
  std::sort(lenghts.begin(),lenghts.end());
  unsigned entries = lenghts.size();
  //integer divion and its residue
  unsigned nr_bins = (unsigned)((entries - entries%nbins) / nbins);
  unsigned residues = entries % nbins;
  if(VERB > 2)
    std::cerr << "\nBinning proteins using equal deepth\n" << std::endl;
  while(residues >= nbins && residues != 0)
  {
    nr_bins += (unsigned)((residues - residues%nbins) / nbins);
    residues = residues % nbins; 
  }
  std::vector<double> values;
  for(unsigned i = 0; i < nbins; i++)
  {
    unsigned index = (unsigned)(nr_bins * i);
    double value = lenghts[index];
    values.push_back(value);
    if(VERB > 2)
      std::cerr << "\nValue of bin : " << i << " = " << value << std::endl;
  }
  //there are some elements at the end that are <= nbins that could not be fitted
  if(residues > 0) values.push_back(lenghts.back());
  else nbins--;
  
  std::multimap<double,std::string>::iterator itlow,itup;
  for(unsigned i = 0; i < nbins; i++)
  {
    double lowerbound = values[i];
    double upperbound = values[i+1];
    itlow = groupedProteins.lower_bound(lowerbound);
    itup =  groupedProteins.upper_bound(upperbound);
    std::set<std::string> proteins;
    std::transform(itlow, itup, std::inserter(proteins,proteins.begin()), RetrieveValue());
    binnedProteins.insert(std::make_pair<unsigned,std::set<std::string> >(i,proteins));
  }
  //FreeAll(values);
  return;
}
    
void ProteinFDRestimator::binProteinsEqualWidth()
{
  //assuming lengths sorted from less to bigger
  std::sort(lenghts.begin(),lenghts.end());
  double min = lenghts.front();
  double max = lenghts.back();
  std::vector<double> values;
  int span = abs(max - min);
  double part = span / nbins;
  if(VERB > 2)
    std::cerr << "\nBinning proteins using equal width\n" << std::endl;
  for(unsigned i = 0; i < nbins; i++)
  {
    unsigned index = (unsigned) min + i*part;
    double value = lenghts[index];
    values.push_back(value);
    if(VERB > 2)
      std::cerr << "\nValue of bin : " << i << " = " << value << std::endl;
  }
  values.push_back(max);
  std::multimap<double,std::string>::iterator itlow,itup;
  for(unsigned i = 0; i < nbins; i++)
  {
    double lowerbound = values[i];
    double upperbound = values[i+1];
    itlow = groupedProteins.lower_bound(lowerbound);
    itup =  groupedProteins.upper_bound(upperbound);
    std::set<std::string> proteins;
    std::transform(itlow, itup, std::inserter(proteins,proteins.begin()), RetrieveValue());
    binnedProteins.insert(std::make_pair<unsigned,std::set<std::string> >(i,proteins));
  }
  //FreeAll(values);
  return;
}

double ProteinFDRestimator::estimatePi0HG(unsigned N,unsigned targets,unsigned cf)
{
  std::vector<double> logprob;
  double finalprob = 0;
  double fdr = 0;
  for(unsigned fp = 0; fp <= cf; fp++)
  {
    unsigned tp = targets - fp;
    unsigned w = N - tp;
    double prob = hypergeometric(fp,N,w,cf);
    logprob.push_back(prob);
  }
  //normalization
  double sum = (double)std::accumulate(logprob.rbegin(), logprob.rend(), 0.0);
  std::transform(logprob.begin(), logprob.end(), 
		 logprob.begin(), std::bind2nd(std::divides<double> (),sum));
  
  //exp probability
  for(unsigned i = 0; i < logprob.size(); i++)
    finalprob += logprob[i] * i;

  //FreeAll(logprob);
  if(isnan(finalprob) || isinf(finalprob)) finalprob = 0.0;
  return finalprob;

}

double ProteinFDRestimator::calculatePepMAss(std::string pepsequence,double charge)
{
  double mass  =  0.0;
  if (pepsequence.length () > minpeplength) {
    
    for(unsigned i=0; i<pepsequence.length();i++)
    {
      if(isalpha(pepsequence[i])){
        mass += massMap_[pepsequence[i]];
      }
    }
    
    mass = (mass + massMap_['o'] + (charge * massMap_['h'])); 
  }
  return mass; 
}


unsigned int ProteinFDRestimator::calculateProtLength(std::string protsequence)
{
  size_t length = protsequence.length();
  std::string peptide;
  std::set<std::string> peptides;
  
  if(protsequence[length-1] != '*'){
    protsequence.push_back('*');
    length++;
  }
  
  for(size_t start=0; start<length; start++)
  {
    if((start == 0) || (protsequence[start] == 'K' && protsequence[start+1] != 'P') || (protsequence[start] == 'R' && protsequence[start+1] != 'P')) 
    {
      int numMisCleavages = 0;  
      for(size_t end=start+1;( (end<length) && (((int)(end-start)) < maxSeqlength) && (numMisCleavages <= missed_cleavages) );end++)
      {
        if((protsequence[end] == 'K') || (protsequence[end] == 'R') || (protsequence[end] == '*'))
	{
          if( (end > length) || (end < length && protsequence[end+1] != 'P')  )
	  {
	    peptide = protsequence.substr(start,((int)(end-start))+2);
	    double  mass = calculatePepMAss(peptide);
	    if((mass > minmass) && (mass< maxmass))
	    {
	      peptides.insert(peptide);
	    }
	    numMisCleavages++;
	  }
        }
      } 
    }
  }
  
  unsigned size = peptide.size();
  //FreeAll(peptide);
  return size;
}


unsigned int ProteinFDRestimator::countProteins(unsigned int bin,const std::set<std::string> &proteins)
{
  std::set<std::string> proteinsBins = binnedProteins[bin];
  unsigned count = 0;
  for(std::set<std::string>::const_iterator it = proteins.begin(); it != proteins.end(); it++)
  {
    std::set<std::string>::iterator itfound = proteinsBins.find(*it);
    if(itfound != proteinsBins.end())
    {
      count++;
      proteinsBins.erase(itfound);
    }
  }
  //FreeAll(proteinsBins);
  return count;
}


unsigned int ProteinFDRestimator::getBinProteins(unsigned int bin)
{
  return binnedProteins[bin].size();
}


unsigned int ProteinFDRestimator::getNumberBins()
{
  return nbins;
}

void ProteinFDRestimator::setDecoyPrefix(std::string prefix)
{
  decoy_prefix = prefix;
}


void ProteinFDRestimator::initMassMap(bool useAvgMass){
 

  if (useAvgMass) /*avg masses*/
    {
      massMap_['h']=  1.00794;  
      massMap_['o']= 15.9994;   

      massMap_['G']= 57.05192;
      massMap_['A']= 71.07880;
      massMap_['S']= 87.07820;
      massMap_['P']= 97.11668;
      massMap_['V']= 99.13256;
      massMap_['T']=101.10508;
      massMap_['C']=103.13880;
      massMap_['L']=113.15944;
      massMap_['I']=113.15944;
      massMap_['X']=113.15944;
      massMap_['N']=114.10384;
      massMap_['O']=114.14720;
      massMap_['B']=114.59622;
      massMap_['D']=115.08860;
      massMap_['Q']=128.13072;
      massMap_['K']=128.17408;
      massMap_['Z']=128.62310;
      massMap_['E']=129.11548;
      massMap_['M']=131.19256;
      massMap_['H']=137.14108;
      massMap_['F']=147.17656;
      massMap_['R']=156.18748;
      massMap_['Y']=163.17596;
      massMap_['W']=186.21320;
    }
  else /* monoisotopic masses */
    {
      massMap_['h']=  1.0078250;
      massMap_['o']= 15.9949146;

      massMap_['A']= 71.0371136;
      massMap_['C']=103.0091854;
      massMap_['D']=115.0269428;
      massMap_['E']=129.0425928;
      massMap_['F']=147.0684136;
      massMap_['G']= 57.0214636;
      massMap_['H']=137.0589116;
      massMap_['I']=113.0840636;
      massMap_['K']=128.0949626;
      massMap_['L']=113.0840636;
      massMap_['M']=131.0404854;
      massMap_['N']=114.0429272;
      massMap_['P']= 97.0527636;
      massMap_['Q']=128.0585772;
      massMap_['R']=156.1011106;
      massMap_['S']= 87.0320282;
      massMap_['T']=101.0476782;
      massMap_['U']=149.90419;
      massMap_['V']= 99.0684136;
      massMap_['W']=186.07931;
      massMap_['Y']=163.06333;
    }

}