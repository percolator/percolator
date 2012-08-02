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

/******************************************************************************************************************/

ProteinFDRestimator::ProteinFDRestimator(std::string __decoy_prefix,unsigned __nbins, 
					   double __targetDecoyRatio, bool __binequalDeepth)
				          :decoy_prefix(__decoy_prefix),nbins(__nbins),
				          targetDecoyRatio(__targetDecoyRatio),binequalDeepth(__binequalDeepth)
{
 
}

ProteinFDRestimator::~ProteinFDRestimator()
{
  FreeAll(binnedProteins);
  FreeAll(groupedProteins);
  FreeAll(lenghts);
}


void ProteinFDRestimator::correctIdenticalSequences(const std::map<std::string,std::pair<std::string,double> > &targetProteins,
						       const std::map<std::string,std::pair<std::string,double> > &decoyProteins)
{
  std::map<std::string,std::pair<std::string,double> >::const_iterator it,it2;
  
  groupedProteins.clear();
  lenghts.clear();
  it = targetProteins.begin();
  it2 = decoyProteins.begin();
  unsigned num_corrected = 0;
  std::set<std::string> previouSeqs;
  double length = 0.0;
  
  for(it = targetProteins.begin(); it != targetProteins.end(); it++)
  {
    std::string targetSeq = (*it).second.first;
    std::string targetName = (*it).first;
    if(previouSeqs.find(targetSeq) != previouSeqs.end())
    {
      length = 0.0;
      num_corrected++;
    }
    else
    {
      length = (*it).second.second;
      previouSeqs.insert(targetSeq);
    }
    groupedProteins.insert(std::make_pair<double,std::string>(length,targetName));
    lenghts.push_back(length);
  }
  
  for(it2 = decoyProteins.begin(); it2 != decoyProteins.end(); it2++)
  {
    std::string decoySeq = (*it2).second.first;
    std::string decoyName = (*it2).first;
    if(previouSeqs.find(decoySeq) != previouSeqs.end())
    {
      length = 0.0;
      num_corrected++;
    }
    else
    {
      length = (*it2).second.second;
      previouSeqs.insert(decoySeq);
    }
    groupedProteins.insert(std::make_pair<double,std::string>(length,decoyName));
    lenghts.push_back(length);
  }
  
  if(VERB > 2)
  {
    std::cerr << "There have been " << num_corrected << " of identical sequences corrected to ''" << std::endl;
  }
  return;
}

void ProteinFDRestimator::groupProteinsGene()
{
  /**not implemented yet**/
}

//NOTE this is a version of the code to work with threads
/*void ProteinFDRestimator::estimateFDRthread(unsigned i,const std::set<std::string> &target, const std::set<std::string> &decoy) 
{ 
  unsigned  numberTP = countProteins(i,target);
  unsigned  numberFP = countProteins(i,decoy);
  unsigned  N = getBinProteins(i);
  double fp = estimatePi0HG(N,numberTP,targetDecoyRatio*numberFP);
  if(VERB > 2)
      std::cerr << "\nEstimating FDR for bin " << i << " in thread " << boost::this_thread::get_id() << " with " << numberFP << " Decoy proteins, " 
	  << numberTP << " Target proteins, and " << N << " Total Proteins in the bin " << " with exp fp " << fp << std::endl;
  fptol += fp;
} */

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
    {
      std::cerr << "\nThere are : " << __target.size() << " target proteins and " << __decoy.size() 
      << " decoys proteins that contains high confident PSMs\n" << std::endl;    
    }

    //NOTE this is a version of the code to work with threads
    /*fptol = 0.0;
    boost::thread t[nbins]; 
    
    for(unsigned i = 0; i < nbins; i++)
    {
       t[i] = boost::thread(boost::bind(&ProteinFDRestimator::estimateFDRthread,this,i,__target,__decoy)); 
       t[i].join();
    }*/
    
    double fptol = 0.0;
    for(unsigned i = 0; i < nbins; i++)
    {
      unsigned numberTP = countProteins(i,__target);
      unsigned numberFP = countProteins(i,__decoy);
      unsigned N = getBinProteins(i);
      double fp = estimatePi0HG(N,numberTP,targetDecoyRatio*numberFP);
      
      if(VERB > 2)
      {
	  std::cerr << "\nEstimating FDR for bin " << i << " with " << numberFP << " Decoy proteins, "
         << numberTP << " Target proteins, and " << N << " Total Proteins in the bin " << " with exp fp " << fp << std::endl;
      }

      fptol += fp;
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
  
  /*while(residues >= nbins && residues != 0)
  {
    nr_bins += (unsigned)((residues - residues%nbins) / nbins);
    residues = residues % nbins; 
  }*/
  
  std::vector<double> values;
  for(unsigned i = 0; i <= nbins; i++)
  {
    unsigned index = (unsigned)(nr_bins * i);
    double value = lenghts[index];
    values.push_back(value);
    if(VERB > 2)
      std::cerr << "\nValue of bin : " << i << " with index " << index << " is " << value << std::endl;
  }
  
  //there are some elements at the end that are <= nbins that could not be fitted
  //FIXME I think it fails if it enters here in some scenarios
  
  if(residues > 0)
  {
    values.back() = lenghts.back();
    if(VERB > 2)
      std::cerr << "\nValue of last bin is fixed to : " << values.back() << std::endl;
  }

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
      std::cerr << "\nValue of bin : " << i << " with index " << index << " is " << value << std::endl;
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

  if(isnan(finalprob) || isinf(finalprob)) finalprob = 0.0;
  return finalprob;

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

bool ProteinFDRestimator::getEqualDeppthBinning()
{
  return binequalDeepth;
}

std::string ProteinFDRestimator::getDecoyPrefix()
{
  return decoy_prefix;
}


double ProteinFDRestimator::getTargetDecoyRatio()
{
  return targetDecoyRatio;
}

void ProteinFDRestimator::setDecoyPrefix(std::string prefix)
{
  decoy_prefix = prefix;
}

void ProteinFDRestimator::setEqualDeepthBinning(bool __equal_deepth)
{
  binequalDeepth = __equal_deepth;
}

void ProteinFDRestimator::setNumberBins(unsigned int __nbins)
{
  nbins = __nbins;
}

void ProteinFDRestimator::setTargetDecoyRatio(double __ratio)
{
  targetDecoyRatio = __ratio;
}

