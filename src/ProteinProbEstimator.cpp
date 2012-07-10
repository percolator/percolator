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

#include <iostream>
#include <fstream>
#include "ProteinProbEstimator.h"

/** Helper functions **/

template<class T> void bootstrap(const vector<T>& in, vector<T>& out,
                                 size_t max_size = 1000) {
  out.clear();
  double n = in.size();
  size_t num_draw = min(in.size(), max_size);
  for (size_t ix = 0; ix < num_draw; ++ix) {
    size_t draw = (size_t)((double)rand() / ((double)RAND_MAX + (double)1) * n);
    out.push_back(in[draw]);
  }
  // sort in desending order
  std::sort(out.begin(), out.end());
}

double antiderivativeAt(double m, double b, double xVal)
{
  return (m*xVal*xVal/2.0 + b*xVal);
}

double squareAntiderivativeAt(double m, double b, double xVal)
{
  double u = m*m;
  double v = 2*m*b;
  double t = b*b;
  return (u*xVal*xVal*xVal/3.0 + v*xVal*xVal/2.0 + t*xVal);
}

double area(double x1, double y1, double x2, double y2, double max_x)
{
  double m = (y2-y1)/(x2-x1);
  double b = y1-m*x1;
  double area =  antiderivativeAt(m, b, min(max_x, x2) ) - antiderivativeAt(m, b, x1);
  if(isnan(area)) return 0.0;
  else return area;
}

double areaSq(double x1, double y1, double x2, double y2, double threshold) {
  double m = (y2-y1)/(x2-x1);
  double b = y1-m*x1;
  double area = squareAntiderivativeAt(m, b, min(threshold, x2) ) - squareAntiderivativeAt(m, b, x1);
  if(isnan(area)) return 0.0;
  else return area;
}

    
ProteinProbEstimator::ProteinProbEstimator(double alpha_par, double beta_par, double gamma_par ,bool __tiesAsOneProtein
			 ,bool __usePi0, bool __outputEmpirQVal, bool __groupProteins, bool __noseparate, bool __noprune, 
			  bool __dogridSearch, unsigned __depth,std::string __decoyPattern, bool __mayufdr,bool __outputDecoys, 
			  bool __tabDelimitedOut, std::string __proteinFN, bool __reduceTree) 
{
  peptideScores = 0;
  proteinGraph = 0;
  gamma = gamma_par;
  alpha = alpha_par;
  beta = beta_par;
  numberDecoyProteins = 0;
  numberTargetProteins = 0;
  pi0 = 1.0;
  tiesAsOneProtein = __tiesAsOneProtein;
  usePi0 = __usePi0;
  outputEmpirQVal = __outputEmpirQVal;
  groupProteins = __groupProteins;
  noseparate = __noseparate;
  noprune = __noprune;
  dogridSearch = __dogridSearch;
  depth = __depth;
  decoyPattern = __decoyPattern;
  mayufdr = __mayufdr;
  fdr = 1.0;
  outputDecoys = __outputDecoys;
  tabDelimitedOut = __tabDelimitedOut;
  proteinFN = __proteinFN;
  rocN = default_rocN; 
  reduceTree = __reduceTree;
}

ProteinProbEstimator::~ProteinProbEstimator()
{
  
  if(proteinGraph)
    delete proteinGraph;
  
  FreeAll(qvalues);
  FreeAll(qvaluesEmp);
  FreeAll(pvalues);
  
  if(mayufdr && fastReader)
    delete fastReader;
  
  for(std::multimap<double,std::vector<std::string> >::iterator it = pepProteins.begin();
      it != pepProteins.end(); it++)
      {
	FreeAll(it->second);
      }
      
  for(std::map<const std::string,Protein*>::iterator it = proteins.begin(); 
      it != proteins.end(); it++)
      {
	if(it->second)
	  delete it->second;
      }
}


bool ProteinProbEstimator::initialize(Scores* fullset){
  
  peptideScores = fullset;
  setTargetandDecoysNames();
  
  double peptidePrior_local = peptidePrior;
  if(computePriors)
  {
    peptidePrior_local = estimatePriors();
  }
  if(VERB > 1 && computePriors)
  {
    std::cerr << "The estimated peptide level prior probability is : " << peptidePrior_local << std::endl;
  }
  
  proteinGraph = new GroupPowerBigraph (alpha,beta,gamma,groupProteins,noseparate,noprune);
  proteinGraph->setMaxAllowedConfigurations(max_allow_configurations);
  proteinGraph->setPeptidePrior(peptidePrior_local);
  
  if(reduceTree && dogridSearch)
  {
    //NOTE lets create a smaller tree to estimate the parameters faster
    if(VERB > 1)
    {
      std::cerr << "Reducing tree of proteins to increase speed.." << std::endl;
    }
    proteinGraph->setProteinThreshold(reduced_proteinThreshold);
    proteinGraph->setPsmThreshold(reduced_psmThreshold);
    proteinGraph->setPeptideThreshold(reduced_peptideThreshold);
    proteinGraph->setGroupProteins(true);
    proteinGraph->setSeparateProteins(false);
    proteinGraph->setPruneProteins(false);
  }
  else
  {
    proteinGraph->setProteinThreshold(proteinThreshold);
    proteinGraph->setPsmThreshold(psmThreshold);
    proteinGraph->setPeptideThreshold(peptideThreshold);
  }
  
  proteinGraph->read(peptideScores);
}

void ProteinProbEstimator::run(){
  
  time_t startTime;
  clock_t startClock;
  time(&startTime);
  startClock = clock();

  if(mayufdr)
  { 
    /** MAYUS method for estimation of Protein FDR **/
    
    std::cerr << "\nEstimating Protein FDR ... " << std::endl;
    
    fastReader = new ProteinFDRestimator();
    fastReader->setDecoyPrefix(decoyPattern);
    fastReader->setTargetDecoyRatio(target_decoy_ratio);
    fastReader->setEqualDeepthBinning(binning_equal_deepth);
    fastReader->setNumberBins(number_bins);

    fastReader->correctIdenticalSequences(targetProteins,decoyProteins);
    //These guys are the number of target and decoys proteins but from the subset of PSM with FDR < threshold
    std::set<std::string> numberTP;
    std::set<std::string> numberFP;
    getTPandPFfromPeptides(psmThresholdMayu,numberTP,numberFP);
    
    double fptol = fastReader->estimateFDR(numberTP,numberFP);
    
    if(fptol == -1)
    {
      fdr = 1.0;
      std::cerr << "\nThere was an error estimating the Protein FDR..\n" << std::endl;
    }
    else
    {	
      fdr = (fptol/(double)numberTP.size());
      
      if(fdr <= 0 || fdr >= 1.0) fdr = 1.0;
      
      if(VERB > 1)
	std::cerr << "\nEstimated Protein FDR at ( " << psmThresholdMayu << ") PSM FDR is : " << fdr << " with " 
	<< fptol << " expected number of false positives proteins\n" << std::endl;
    }
    
    /** MAYUS method for estimationg of Protein FDR**/
  }
  
  if(dogridSearch) 
  {
    if(VERB > 1) 
    {
      std::cerr << "\nThe parameters for the model will be estimated by grid search.\n" << std::endl;
    }
    
    gridSearch(alpha,gamma,beta);
    
    time_t procStart;
    clock_t procStartClock = clock();
    time(&procStart);
    double diff = difftime(procStart, startTime);
    if (VERB > 1) cerr << "\nEstimating the parameters took : "
      << ((double)(procStartClock - startClock)) / (double)CLOCKS_PER_SEC
      << " cpu seconds or " << diff << " seconds wall time\n" << endl;
  }

  if(VERB > 1) 
  {
      cerr << "\nThe following parameters have been chosen:\n";
      cerr << "gamma = " << gamma << endl;
      cerr << "alpha = " << alpha << endl;
      cerr << "beta  = " << beta << endl;
      cerr << "\nProtein level probabilities will now be estimated\n";
  }


  if(dogridSearch && reduceTree)
  {
    //NOTE lets create the tree again with all the members
    proteinGraph->setProteinThreshold(proteinThreshold);
    proteinGraph->setPsmThreshold(psmThreshold);
    proteinGraph->setPeptideThreshold(peptideThreshold);
    proteinGraph->setGroupProteins(groupProteins);
    proteinGraph->setSeparateProteins(noseparate);
    proteinGraph->setPruneProteins(noprune);
    proteinGraph->read(peptideScores);
  }
  
  proteinGraph->setAlphaBetaGamma(alpha,beta,gamma);
  proteinGraph->getProteinProbs();
  pepProteins.clear();
  proteinGraph->getProteinProbsPercolator(pepProteins);
  estimateQValues();
  
  if(usePi0 && !mayufdr && ProteinProbEstimator::getOutputEmpirQval())
  {
    estimatePValues();
    pi0 = estimatePi0();
    if(pi0 <= 0.0 || pi0 > 1.0) pi0 = *qvalues.rbegin();
  }
  else
  {
    pi0 = fdr;
  }
  
  estimateQValuesEmp();
  updateProteinProbabilities();

  if(VERB > 1)
  {
    std::cerr << "\nThe number of Proteins idenfified at q-value = 0.01 is : " << getQvaluesBelowLevel(0.01) << std::endl;
  }
  
  time_t procStart;
  clock_t procStartClock = clock();
  time(&procStart);
  double diff = difftime(procStart, startTime);
  if (VERB > 1) cerr << "Estimating Protein Probabilities took : "
    << ((double)(procStartClock - startClock)) / (double)CLOCKS_PER_SEC
    << " cpu seconds or " << diff << " seconds wall time" << endl;
  

   if(tabDelimitedOut && !proteinFN.empty())
   {
      ofstream proteinOut(proteinFN.data(), ios::out);
      print(proteinOut);
      proteinOut.close();
   }
   else
   {
      print(std::cout);
   }
  
}

double ProteinProbEstimator::estimatePriors()
{
  /* Compute a priori probabilities of peptide presence */
  /* it is the mean of the probabilities, maybe one prior for each charge */
  double prior_peptide = 0.0;
  double prior_peptide2 = 0.0;
  unsigned total_peptides = 0;
  double prior = 0.0;
  for (vector<ScoreHolder>::iterator psm = peptideScores->begin(); psm!= peptideScores->end(); ++psm) 
  {
    if(!psm->isDecoy())
    {
      unsigned size = psm->pPSM->proteinIds.size();
      double prior = prior_protein * size;
      double tmp_prior = prior;
      // for each protein
      for(set<string>::iterator protIt = psm->pPSM->proteinIds.begin(); protIt != psm->pPSM->proteinIds.end(); protIt++)
      {
	unsigned index = std::distance(psm->pPSM->proteinIds.begin(), protIt);
	tmp_prior = (tmp_prior * prior_protein * (size - index)) / (index + 1);
	prior +=  pow(-1,index) * tmp_prior;
      }
      /* update computed prior */
      prior_peptide += prior;
      prior_peptide2 += psm->pPSM->pep;
      ++total_peptides;
    }
  }
  //prior = prior_peptide2 / (double)total_peptides;
  prior = prior_peptide / (double)total_peptides;
  if(prior > 0.99) prior = 0.99;
  if(prior < 0.01) prior = 0.01;
  return prior;
}



void ProteinProbEstimator::print(ostream& myout) {
  
  std::vector<std::pair<std::string,Protein*> > myvec(proteins.begin(), proteins.end());
  std::sort(myvec.begin(), myvec.end(), IntCmpProb());
  
  myout
      << "ProteinId\tq-value\tposterior_error_prob\tpeptideIds"
      << std::endl;
      
  for (std::vector<std::pair<std::string,Protein*> > ::const_iterator myP = myvec.begin(); 
	 myP != myvec.end(); myP++) 
  {
    if( (!outputDecoys && !myP->second->getIsDecoy()) || (outputDecoys))
    {
      myout << myP->second->getName() << "\t" << myP->second->getQ() << "\t" << myP->second->getPEP() << "\t";
      std::vector<Protein::Peptide*> peptides = myP->second->getPeptides();
      for(std::vector<Protein::Peptide*>::const_iterator peptIt = peptides.begin(); peptIt != peptides.end(); peptIt++)
      {
	 if((*peptIt)->name != "")
	 {
	    myout << (*peptIt)->name << "  ";
	 }
      }
      myout << std::endl;
    }
  }
}

void ProteinProbEstimator::estimatePValues()
{
  // assuming combined sorted in best hit first order
  std::vector<std::pair<double , bool> > combined;
  for (std::multimap<double,std::vector<std::string> >::const_iterator it = pepProteins.begin();
       it != pepProteins.end(); it++)
  {
     double prob = it->first;
     std::vector<std::string> proteinList = it->second;
     for(std::vector<std::string>::const_iterator itP = proteinList.begin();
	itP != proteinList.end(); itP++)
      {
	std::string proteinName = *itP;
	bool isdecoy = proteins[proteinName]->getIsDecoy();
	combined.push_back(std::make_pair<double,bool>(prob,isdecoy));
      }
  }
  pvalues.clear();
  std::vector<pair<double, bool> >::const_iterator myPair = combined.begin();
  size_t nDecoys = 0, posSame = 0, negSame = 0;
  double prevScore = -4711.4711; // number that hopefully never turn up first in sequence
  while (myPair != combined.end()) {
    if (myPair->first != prevScore) {
      for (size_t ix = 0; ix < posSame; ++ix) {
        pvalues.push_back((double)nDecoys + (((double)negSame)
            / (double)(posSame + 1)) * (ix + 1));
      }
      nDecoys += negSame;
      negSame = 0;
      posSame = 0;
      prevScore = myPair->first;
    }
    if (myPair->second) {
      ++negSame;
    } else {
      ++posSame;
    }
    ++myPair;
  }
  std::transform(pvalues.begin(), pvalues.end(), pvalues.begin(), 
		 std::bind2nd(std::divides<double> (),(double)nDecoys));
}

void ProteinProbEstimator::getTPandPFfromPeptides(double threshold, std::set<std::string> &numberTP, 
						  std::set<std::string> &numberFP)
{
  /* The original paper of Mayu describes a protein as :
   * FP = if only if all its peptides with q <= threshold are decoy
   * TP = at least one of its peptides with q <= threshold is target
   * However, the way mayu estimates it on the program is like this :
   * FP = any protein that contains a decoy psm with q <= threshold
   * TP = any protein that contains a target psm with q <= threshold
   * They do not consider protein containing both decoy and target psms
   * Also, mayu estimates q as the empirical (target-decoy) q value.
   * Percolator estimates q as the empirical (target-decoy) q value and adjusted by pi0
   * Mayu extract the list of TP and FP proteins from PSM level whereas percolator
   * extract the list of TP and FP proteins from peptide level, this avoids redundancy and
   * gives a better calibration since peptide level q values are re-adjusted in percolator.
   * This creates sometimes a difference in the number of TP and FP proteins which causes a difference
   * in the estimated protein FDR
   */
  for (std::map<std::string,Protein*>::const_iterator it = proteins.begin();
       it != proteins.end(); it++)
  {
     unsigned num_target_confident = 0;
     unsigned num_decoy_confident = 0;
     std::string protname = it->first;
     std::vector<Protein::Peptide*> peptides = it->second->getPeptides();
     for(std::vector<Protein::Peptide*>::const_iterator itP = peptides.begin();
	itP != peptides.end(); itP++)
      {
	Protein::Peptide *p = *itP;
	if(p->q <= threshold && p->isdecoy)
	  ++num_decoy_confident;
	if(p->q <= threshold && !p->isdecoy)
	//else if(p->q <= threshold && !p->isdecoy)
	  ++num_target_confident;
	
	
      }
      
      //if(num_decoy_confident > 0 && num_target_confident == 0)
      if(num_decoy_confident > 0)
      {
	numberFP.insert(protname);
      }
      if(num_target_confident > 0)
      //else if(num_target_confident > 0)
      {
	numberTP.insert(protname);
      }
  }
  return;
}


double ProteinProbEstimator::estimatePi0(const unsigned int numBoot) 
{
  std::vector<double> pBoot, lambdas, pi0s, mse;
  std::vector<double>::iterator start;
  int numLambda = 100;
  double maxLambda = 0.5;
  size_t n = pvalues.size();
  // Calculate pi0 for different values for lambda
  // N.B. numLambda and maxLambda are global variables.
  for (unsigned int ix = 0; ix <= numLambda; ++ix) 
  {
    double lambda = ((ix + 1) / (double)numLambda) * maxLambda;
    // Find the index of the first element in p that is < lambda.
    // N.B. Assumes p is sorted in ascending order.
    start = lower_bound(pvalues.begin(), pvalues.end(), lambda);
    // Calculates the difference in index between start and end
    double Wl = (double)distance(start, pvalues.end());
    double pi0 = Wl / n / (1 - lambda);
    if (pi0 > 0.0) 
    {
      lambdas.push_back(lambda);
      pi0s.push_back(pi0);
    }
  }
  if(pi0s.size()==0)
  {
    cerr << "Error in the input data: too good separation between target "
        << "and decoy Proteins.\nImpossible to estimate pi0. Taking the highest estimated q value as pi0.\n";
    return -1;
  }
  double minPi0 = *min_element(pi0s.begin(), pi0s.end());
  // Initialize the vector mse with zeroes.
  fill_n(back_inserter(mse), pi0s.size(), 0.0);
  // Examine which lambda level that is most stable under bootstrap
  for (unsigned int boot = 0; boot < numBoot; ++boot) 
  {
    // Create an array of bootstrapped p-values, and sort in ascending order.
    bootstrap<double> (pvalues, pBoot);
    n = pBoot.size();
    for (unsigned int ix = 0; ix < lambdas.size(); ++ix) 
    {
      start = lower_bound(pBoot.begin(), pBoot.end(), lambdas[ix]);
      double Wl = (double)distance(start, pBoot.end());
      double pi0Boot = Wl / n / (1 - lambdas[ix]);
      // Estimated mean-squared error.
      mse[ix] += (pi0Boot - minPi0) * (pi0Boot - minPi0);
    }
  }
  // Which index did the iterator get?
  unsigned int minIx = distance(mse.begin(), min_element(mse.begin(),mse.end()));
  double pi0 = max(min(pi0s[minIx], 1.0), 0.0);
  return pi0;
}

unsigned ProteinProbEstimator::getQvaluesBelowLevel(double level)
{   
    unsigned nP = 0;
    for (std::map<const std::string,Protein*>::const_iterator myP = proteins.begin(); 
	 myP != proteins.end(); ++myP) {
	 if(myP->second->getQ() <= level && !myP->second->getIsDecoy()) nP++;
    }
    return nP;
}

unsigned ProteinProbEstimator::getQvaluesBelowLevelDecoy(double level)
{   
    unsigned nP = 0;
    for (std::map<const std::string,Protein*>::const_iterator myP = proteins.begin(); 
	 myP != proteins.end(); ++myP) {
	 if(myP->second->getQ() <= level && myP->second->getIsDecoy()) nP++;
    }
    return nP;
}


void ProteinProbEstimator::estimateQValues()
{
  unsigned nP = 0;
  double sum = 0.0;
  double qvalue = 0.0;
  qvalues.clear();

  for (std::multimap<double,std::vector<std::string> >::const_iterator it = pepProteins.begin(); 
       it != pepProteins.end(); it++) 
  {
    
    if(tiesAsOneProtein)
    {
      int ntargets = countTargets(it->second);
      //NOTE in case I want to count and use target and decoys proteins while estimateing qvalue from PEP
      if(countDecoyQvalue)
      {
	int ndecoys = it->second.size() - ntargets;
	sum += (double)(it->first * (ntargets + ndecoys));
	nP += (ntargets + ndecoys);
      }
      else
      {
	sum += (double)(it->first * ntargets);
	nP += ntargets;
      }
      
      qvalue = (sum / (double)nP);
      if(isnan(qvalue) || isinf(qvalue) || qvalue > 1.0) qvalue = 1.0;
      qvalues.push_back(qvalue);
    }    
    else
    {
      std::vector<std::string> proteins = it->second;
      for(std::vector<std::string>::const_iterator it2 = proteins.begin(); it2 != proteins.end(); it2++)
      {
	std::string protein = *it2;
	if(!countDecoyQvalue)
	{
	  if(isTarget(protein))
	  {
	    sum += it->first;
	    nP++;
	  }
	}
	else
	{
	  //NOTE in case I want to count and use target and decoys proteins while estimateing qvalue from PEP
	  sum += it->first;
	  nP++;
	}
	qvalue = (sum / (double)nP);
	if(isnan(qvalue) || isinf(qvalue) || qvalue > 1.0) qvalue = 1.0;
	qvalues.push_back(qvalue);
      }
    }

  }
  std::partial_sum(qvalues.rbegin(),qvalues.rend(),qvalues.rbegin(),myminfunc);
}

void ProteinProbEstimator::estimateQValuesEmp()
{
    // assuming combined sorted in decending order
  unsigned nDecoys = 0;
  unsigned numTarget = 0;
  unsigned nTargets = 0;
  double qvalue = 0.0;
  unsigned numDecoy = 0;
  pvalues.clear();
  qvaluesEmp.clear();
  double TargetDecoyRatio = (double)numberTargetProteins / (double)numberDecoyProteins;
 
  for (std::multimap<double,std::vector<std::string> >::const_iterator it = pepProteins.begin(); 
       it != pepProteins.end(); it++) 
  {

    if(tiesAsOneProtein)
    {
      numTarget = countTargets(it->second);
      numDecoy = it->second.size() - numTarget;
      
      nDecoys += numDecoy;
      nTargets += numTarget;
      
      if(nTargets) qvalue = (double)(nDecoys * pi0 * TargetDecoyRatio) / (double)nTargets;
      if(isnan(qvalue) || isinf(qvalue) || qvalue > 1.0) qvalue = 1.0;
      
      qvaluesEmp.push_back(qvalue);
      
      if(numDecoy > 0)
        pvalues.push_back((nDecoys)/(double)(numberDecoyProteins));
      else 
        pvalues.push_back((nDecoys+(double)1)/(numberDecoyProteins+(double)1));
    }
    else
    {
      std::vector<std::string> proteins = it->second;
      for(std::vector<std::string>::const_iterator it2 = proteins.begin(); it2 != proteins.end(); it2++)
      {
	 std::string protein = *it2;
	 if(isDecoy(protein))
	 {  
	   nDecoys++;
	   pvalues.push_back((nDecoys)/(double)(numberDecoyProteins));
	 }
	 else
	 {
	   nTargets++;
	   pvalues.push_back((nDecoys+(double)1)/(numberDecoyProteins+(double)1));
	 }
	 
	 if(nTargets) qvalue = (double)(nDecoys * pi0 * TargetDecoyRatio) / (double)nTargets;
	 if(isnan(qvalue) || isinf(qvalue) || qvalue > 1.0) qvalue = 1.0;
	 qvaluesEmp.push_back(qvalue);
      }
    }
  }
  std::partial_sum(qvaluesEmp.rbegin(), qvaluesEmp.rend(), qvaluesEmp.rbegin(), myminfunc);
}

void ProteinProbEstimator::updateProteinProbabilities()
{
  std::vector<double> peps;
  std::vector<std::vector<std::string> > proteinNames;
  std::transform(pepProteins.begin(), pepProteins.end(), std::back_inserter(peps), RetrieveKey());
  std::transform(pepProteins.begin(), pepProteins.end(), std::back_inserter(proteinNames), RetrieveValue());
  unsigned qindex = 0;
  for (unsigned i = 0; i < peps.size(); i++) 
  {
    double pep = peps[i];
    std::vector<std::string> proteinlist = proteinNames[i];
    for(unsigned j = 0; j < proteinlist.size(); j++)
    { 
      std::string proteinName = proteinlist[j];
      if(tiesAsOneProtein)
      {
	proteins[proteinName]->setPEP(pep);
	proteins[proteinName]->setQ(qvalues[i]);
	proteins[proteinName]->setQemp(qvaluesEmp[i]);
	proteins[proteinName]->setP(pvalues[i]);
      }
      else
      {	
	proteins[proteinName]->setPEP(pep);
	proteins[proteinName]->setQ(qvalues[qindex]);
	proteins[proteinName]->setQemp(qvaluesEmp[qindex]);
	proteins[proteinName]->setP(pvalues[qindex]);
      }
      qindex++;
    }
  }

}



std::map<const std::string,Protein*> ProteinProbEstimator::getProteins()
{
  return this->proteins;
}


void ProteinProbEstimator::setTargetandDecoysNames()
{
  for (vector<ScoreHolder>::iterator psm = peptideScores->begin(); psm!= peptideScores->end(); ++psm) 
  {    
    // for each protein
    for(set<string>::iterator protIt = psm->pPSM->proteinIds.begin(); protIt != psm->pPSM->proteinIds.end(); protIt++)
    {
      Protein::Peptide *peptide = new Protein::Peptide(psm->pPSM->getPeptideSequence(),psm->isDecoy(),
							psm->pPSM->pep,psm->pPSM->q,psm->pPSM->p);
      if(proteins.find(*protIt) == proteins.end())
      {
	Protein *newprotein = new Protein(*protIt,0.0,0.0,0.0,0.0,psm->isDecoy(),peptide);
	proteins.insert(std::make_pair<std::string,Protein*>(*protIt,newprotein));
	
	if(psm->isDecoy())
	{
	  falsePosSet.insert(*protIt);
	}
	else
	{
	  truePosSet.insert(*protIt);
	}
      }
      else
      {
	proteins[*protIt]->setPeptide(peptide);
      }
    }
  }  
  numberDecoyProteins = falsePosSet.size();
  numberTargetProteins = truePosSet.size();
}


void ProteinProbEstimator::addProteinDb(const percolatorInNs::protein& protein)
{
  if(protein.isDecoy())
    decoyProteins.insert(std::make_pair<std::string,std::pair<std::string,double> >
    (protein.name(),std::make_pair<std::string,double>(protein.sequence(),protein.length())));
  else
    targetProteins.insert(std::make_pair<std::string,std::pair<std::string,double> >
    (protein.name(),std::make_pair<std::string,double>(protein.sequence(),protein.length())));
}

void ProteinProbEstimator::gridSearch(double __alpha,double __gamma,double __beta)
{
 
  double gamma_best, alpha_best, beta_best;
  gamma_best = alpha_best = beta_best = -1.0;
  double best_objective = -100000000;
  std::vector<std::vector<std::string> > names;
  std::vector<double> probs,empq,estq;
  std::vector<unsigned> numberFP,numberTP;  
  std::vector<double> gamma_search,beta_search,alpha_search;

  switch(depth)
  {
    case 0:
      gamma_search = boost::assign::list_of(0.1)(0.25)(0.5)(0.75);
      beta_search = boost::assign::list_of(0.0)(0.01)(0.015)(0.025)(0.035)(0.05)(0.1);
      alpha_search = boost::assign::list_of(0.01)(0.04)(0.09)(0.16)(0.25)(0.36)(0.5);
      break;
    
    case 1:
      gamma_search = boost::assign::list_of(0.1)(0.25)(0.5);
      beta_search = boost::assign::list_of(0.0)(0.01)(0.15)(0.025)(0.035)(0.05);
      alpha_search = boost::assign::list_of(0.01)(0.04)(0.09)(0.16)(0.25)(0.36);
      break;
      
    case 2:
      gamma_search = boost::assign::list_of(0.1)(0.5);
      beta_search = boost::assign::list_of(0.0)(0.01)(0.15)(0.030)(0.05);
      alpha_search = boost::assign::list_of(0.01)(0.04)(0.16)(0.25)(0.36);
      break;
    
    default:
      gamma_search = boost::assign::list_of(0.5);
      beta_search = boost::assign::list_of(0.0)(0.01)(0.15)(0.030)(0.05);
      alpha_search = boost::assign::list_of(0.01)(0.04)(0.16)(0.25)(0.36);
  }

  if(__alpha != -1)
    alpha_search = boost::assign::list_of(__alpha);
  if(__beta != -1)
    beta_search = boost::assign::list_of(__beta);
  if(__gamma != -1)
    gamma_search = boost::assign::list_of(__gamma);
  
  //NOTE paralellize it for gamma, build copy constructor for fido, be careful with shared variables (mutex)
  
  for (unsigned int i = 0; i < gamma_search.size(); i++)
  {
    for (unsigned int j = 0; j < alpha_search.size(); j++)
    {
      for (unsigned int k = 0; k < beta_search.size(); k++)
      {
	double gamma_local = gamma_search[i];
	double alpha_local = alpha_search[j];
	double beta_local = beta_search[k];
	
	if(VERB > 2)
	{
	  std::cerr << "Grid searching Alpha= " << alpha_local << " Beta= " << beta_local << " Gamma= " 
	  << gamma_local << std::endl;
	}
	
	proteinGraph->setAlphaBetaGamma(alpha_local, beta_local, gamma_local);
	proteinGraph->getProteinProbs();
	proteinGraph->getProteinProbsAndNames(names,probs);
	getEstimated_and_Empirical_FDR(names,probs,empq,estq);
	getROC(names,numberFP,numberTP);
	
	double rocR = getROC_N(numberFP, numberTP, rocN);
	double fdr_mse = getFDR_divergence(estq, empq, threshold);
	double current_objective = (lambda * rocR) - fabs(((1-lambda) * (fdr_mse)));
	
	if (current_objective > best_objective)
	{
	  best_objective = current_objective;
	  gamma_best = gamma_local;
	  alpha_best = alpha_local;
	  beta_best = beta_local;
	}
	
	if(VERB > 2)
	{
	  std::cerr << "Roc " << rocN <<" , MSE and objective function value " << " : " << rocR << " " 
	  << fabs(fdr_mse) << " " << current_objective << std::endl;
	}
      }
    }
  }
  
  alpha = alpha_best;
  beta = beta_best;
  gamma = gamma_best;
}


/**
 * output protein level probabilites results in xml format
 */
void ProteinProbEstimator::writeOutputToXML(string xmlOutputFN){
  
  
  std::vector<std::pair<std::string,Protein*> > myvec(proteins.begin(), proteins.end());
  std::sort(myvec.begin(), myvec.end(), IntCmpProb());

  ofstream os;
  os.open(xmlOutputFN.data(), ios::app);
  // append PROTEINs tag
  os << "  <proteins>" << endl;
  for (std::vector<std::pair<std::string,Protein*> > ::const_iterator myP = myvec.begin(); 
	 myP != myvec.end(); myP++) {
     
        if( (!outputDecoys && !myP->second->getIsDecoy()) || (outputDecoys))
	{

	  os << "    <protein p:protein_id=\"" << myP->second->getName() << "\"";
  
	  if (outputDecoys) {
	    if(myP->second->getIsDecoy()) os << " p:decoy=\"true\"";
	    else  os << " p:decoy=\"false\"";
	  }
	  os << ">" << endl;
	  os << "      <pep>" << myP->second->getPEP() << "</pep>" << endl;
	  if(ProteinProbEstimator::getOutputEmpirQval())
	    os << "      <q_value_emp>" << myP->second->getQemp() << "</q_value_emp>\n";
	  os << "      <q_value>" << myP->second->getQ() << "</q_value>\n";
	  if(ProteinProbEstimator::getOutputEmpirQval())
	    os << "      <p_value>" << myP->second->getP() << "</p_value>\n";
	  std::vector<Protein::Peptide*> peptides = myP->second->getPeptides();
	  for(std::vector<Protein::Peptide*>::const_iterator peptIt = peptides.begin(); 
	      peptIt != peptides.end(); peptIt++)
	  {
	    if((*peptIt)->name != "")
	    {
	      os << "      <peptide_seq seq=\"" << (*peptIt)->name << "\"/>"<<endl;
	    }
	    
	  }
	  os << "    </protein>" << endl;
	}
  }
    
  os << "  </proteins>" << endl << endl;
  os.close();
  
  FreeAll(myvec);
}


string ProteinProbEstimator::printCopyright(){
  ostringstream oss;
  oss << "Copyright (c) 2008-9 University of Washington. All rights reserved.\n"
      << "Written by Oliver R. Serang (orserang@u.washington.edu) in the\n"
      << "Department of Genome Sciences at the University of Washington.\n" << std::endl;
  return oss.str();
}


double ProteinProbEstimator::getROC_N(const std::vector<unsigned> &fpArray, const std::vector<unsigned> &tpArray, int N)
{
  double rocNvalue = 0.0;

  if ( fpArray.back() < N )
    {
      std::cerr << "There are not enough false positives; needed " << N << " and was only given " 
      << fpArray.back() << std::endl << std::endl;
      exit(1);
    }

  for (int k=0;( (k<fpArray.size()-1) && (fpArray[k] < N)); k++)
  {
      if ( fpArray[k] != fpArray[k+1] )
	{
	  double currentArea = area(fpArray[k], tpArray[k], fpArray[k+1], tpArray[k+1], N);
	  rocNvalue += currentArea;
	}
  }
  return rocNvalue / (N * tpArray.back());
}

void ProteinProbEstimator::getEstimated_and_Empirical_FDR(const std::vector<std::vector<string> > &names,
							       const std::vector<double> &probabilities,
							       std::vector<double> &empq,
							       std::vector<double> &estq) 
{
  empq.clear();
  estq.clear();
  double fpCount = 0.0, tpCount = 0.0;
  double totalFDR = 0.0, estFDR = 0.0, empFDR = 0.0;
  double TargetDecoyRatio = (double)numberTargetProteins / (double)numberDecoyProteins;
  double previousEmpQ = 0.0;
  double previousEstQ = 0.0;
  if(updateRocN) rocN = 0;
  //NOTE no need to store more q values since they will not be taken into account while estimating MSE FDR divergence
  for (int k=0; (k<names.size() && estFDR <= threshold); k++)
    {
      double prob = probabilities[k];

      if(tiesAsOneProtein)
      {
	unsigned tpChange = countTargets(names[k]);
	unsigned fpChange = names[k].size() - tpChange;
	fpCount += (double)fpChange;
	tpCount += (double)tpChange;
	
	if(countDecoyQvalue)
	{
	  //NOTE in case I want to count target and decoys while estimateing qvalue from PEP
	  totalFDR += (prob) * (double)(tpChange + fpChange);
	  estFDR = totalFDR / (tpCount + fpCount);
	}
	else
	{
	  totalFDR += (prob) * (double)(tpChange);
	  estFDR = totalFDR / (tpCount);  
	}

	if(tpCount) empFDR = (fpCount * pi0 * TargetDecoyRatio) / tpCount; 
	
	if(empFDR > 1.0 || isnan(empFDR) || isinf(empFDR)) empFDR = 1.0;
	if(estFDR > 1.0 || isnan(estFDR) || isinf(estFDR)) estFDR = 1.0;
	    
	if(estFDR < previousEstQ) estFDR = previousEstQ;
	else previousEstQ = estFDR;
	    
	if(empFDR < previousEmpQ) empFDR = previousEmpQ;
	else previousEmpQ = empFDR;
	
	if(estFDR <= thresholdRoc && updateRocN)
	{ 
	  rocN = (unsigned)std::max(rocN,(unsigned)std::max(50,std::min((int)fpCount,500)));
	}
	estq.push_back(estFDR);
	empq.push_back(empFDR);
	
      }
      else
      {
	for(unsigned i=0; i<names[k].size(); i++)
	{
	    std::string protein = names[k][i];
	    if(isDecoy(protein))
	    {
	      fpCount++;
	    }
	    else
	    {
	      tpCount++;
	      if(!countDecoyQvalue)
	      {
		totalFDR += (prob);
		estFDR = totalFDR / (tpCount);
	      }
	    }
	    if(countDecoyQvalue)
	    {
	      totalFDR += (prob);
	      estFDR = totalFDR / (tpCount + fpCount);
	    }
	    
	    if(tpCount) empFDR = (fpCount * pi0 * TargetDecoyRatio) / tpCount; 
	    if(empFDR > 1.0 || isnan(empFDR) || isinf(empFDR)) empFDR = 1.0;
	    if(estFDR > 1.0 || isnan(estFDR) || isinf(estFDR)) estFDR = 1.0;
	    
	    if(estFDR < previousEstQ) estFDR = previousEstQ;
	    else previousEstQ = estFDR;
	    
	    if(empFDR < previousEmpQ) empFDR = previousEmpQ;
	    else previousEmpQ = empFDR;
	    
	    if(estFDR <= thresholdRoc && updateRocN)
	    {
	      rocN = (unsigned)std::max(rocN,(unsigned)std::max(50,std::min((int)fpCount,500)));
	    }
	    estq.push_back(estFDR);
	    empq.push_back(empFDR);
	    
	 }
      }
	
    }
   
  return;
}


double ProteinProbEstimator::getFDR_divergence(const std::vector<double> &estFDR, const std::vector<double> &empFDR, double THRESH)
{
  Vector diff = Vector(estFDR) - Vector(empFDR);
  double tot = 0.0;
  int k;
  for( k=0; k<diff.size()-1; k++)
  {
      if(conservative)
	tot += area(estFDR[k], diff[k], estFDR[k+1], diff[k+1], estFDR[k+1]);
      else
	tot += areaSq(estFDR[k], diff[k], estFDR[k+1], diff[k+1], estFDR[k+1]);
  }

  double xRange = min(THRESH, estFDR[k]) - estFDR[0];

  if ( isinf(tot) || tot == 0.0)
    return tot;
  else
    return tot / xRange;
}


void ProteinProbEstimator::getROC(const std::vector<std::vector<string> > &names,
				  std::vector<unsigned> &numberFP,std::vector<unsigned> &numberTP)
{
  numberFP.clear();
  numberTP.clear();
  unsigned fpCount, tpCount;
  fpCount = tpCount = 0;
  //NOTE no need to store more fp since they will not be taken into account while estimating the ROC curve divergence
  for (int k=0; (k<names.size() && fpCount <= rocN); k++)
  {
    
      unsigned tpChange = countTargets(names[k]);
      unsigned fpChange = names[k].size() - tpChange;
      
      //NOTE possible alternative is to only sum up when the new prob is different that the previous one
      fpCount += fpChange;
      tpCount += tpChange;
      
      numberFP.push_back( fpCount );
      numberTP.push_back( tpCount );

  }

  numberFP.push_back( fpCount );
  numberTP.push_back( tpCount );	  
  numberFP.push_back( falsePosSet.size() );
  numberTP.push_back( truePosSet.size() );
  
  return;
}

unsigned ProteinProbEstimator::countTargets(const std::vector<std::string> &proteinList)
{
  //NOTE this is probably faster
  //return std::count_if(proteinList.begin(),proteinList.end(),isTarget);
 
  unsigned count = 0;
  for(std::vector<std::string>::const_iterator it = proteinList.begin();
      it != proteinList.end(); it++)
  {
    if(useDecoyPrefix)
    {
      if((*it).find(decoyPattern) == std::string::npos)
      {
	count++;
      }
    }
    else
    {
      if(truePosSet.count(*it) != 0)
      {
	count++;
      }
    }
  }
  return count;
}

unsigned ProteinProbEstimator::countDecoys(const std::vector<std::string> &proteinList)
{
  //NOTE this is probably faster
  //return std::count_if(proteinList.begin(),proteinList.end(),isDecoy);
  
  unsigned count = 0;
  for(std::vector<std::string>::const_iterator it = proteinList.begin();
      it != proteinList.end(); it++)
  {
    if(useDecoyPrefix)
    {
      if((*it).find(decoyPattern) != std::string::npos)
      {
	count++;
      }
    }
    else
    {
      if(falsePosSet.count(*it) != 0)
      {
	count++;
      }
    }
  }
  return count;
}

bool ProteinProbEstimator::isDecoy(const std::string& proteinName) 
{
  //NOTE faster with decoyPrefix but I assume the label that identifies decoys is in decoyPattern
   bool isdecoy;
   if(!useDecoyPrefix)
     isdecoy = falsePosSet.count(proteinName) != 0;
   else
     isdecoy = proteinName.find(decoyPattern) != std::string::npos;
   
   return isdecoy;
}
    
bool ProteinProbEstimator::isTarget(const std::string& proteinName) 
{  
  //NOTE faster with decoyPrefix but I assume the label that identifies decoys is in decoyPattern
   bool istarget;
   if(!useDecoyPrefix)
     istarget = truePosSet.count(proteinName) != 0;
   else
     istarget = proteinName.find(decoyPattern) == std::string::npos;
   
   return istarget;
}

void ProteinProbEstimator::setOutputEmpirQval(bool __outputEmpirQVal)
{
  outputEmpirQVal = __outputEmpirQVal;
}

void ProteinProbEstimator::setTiesAsOneProtein(bool __tiesAsOneProtein)
{
  tiesAsOneProtein = __tiesAsOneProtein;
}

void ProteinProbEstimator::setUsePio(bool __usePi0)
{
  usePi0 = __usePi0;
}

void ProteinProbEstimator::setGroupProteins(bool __groupProteins)
{
  groupProteins = __groupProteins;
}

void ProteinProbEstimator::setPruneProteins(bool __noprune)
{
  noprune = __noprune;
}

void ProteinProbEstimator::setSeparateProteins(bool __noseparate)
{
  noseparate = __noseparate;
}

bool ProteinProbEstimator::getOutputEmpirQval()
{
  return outputEmpirQVal;
}

void ProteinProbEstimator::setDepth(unsigned int __depth)
{
  depth = __depth;
}

void ProteinProbEstimator::setGridSearch(bool __dogridSearch)
{
  dogridSearch = __dogridSearch;
}

void ProteinProbEstimator::setMayusFDR(bool __mayufdr)
{
  mayufdr = __mayufdr;
}

void ProteinProbEstimator::setFDR(double __fdr)
{
  fdr = __fdr;
}

void ProteinProbEstimator::setOutputDecoys(bool __outputDecoys)
{
  outputDecoys = __outputDecoys;
}

void ProteinProbEstimator::setProteinFN(std::string __proteinFN)
{
  proteinFN = __proteinFN;
}

void ProteinProbEstimator::setTabDelimitedOutput(bool __tabDelimitedOut)
{
  tabDelimitedOut = __tabDelimitedOut;
}

bool ProteinProbEstimator::getTiesAsOneProtein()
{
  return tiesAsOneProtein;
}

bool ProteinProbEstimator::getUsePio()
{
  return usePi0;
}

double ProteinProbEstimator::getPi0()
{
  return pi0;
}

bool ProteinProbEstimator::getGroupProteins()
{
  return groupProteins;
}


bool ProteinProbEstimator::getPruneProteins()
{
  return noprune;
}

bool ProteinProbEstimator::getSeparateProteins()
{
  return noseparate;
}

double ProteinProbEstimator::getAlpha()
{
 return alpha;
}

double ProteinProbEstimator::getBeta()
{
  return beta;
}

double ProteinProbEstimator::getGamma()
{
  return gamma;
}


string ProteinProbEstimator::getDecoyPatter()
{
  return decoyPattern;
}

bool ProteinProbEstimator::getDepth()
{
  return depth;
}


bool ProteinProbEstimator::getMayuFdr()
{
  return mayufdr;
}

double ProteinProbEstimator::getFDR()
{
  return fdr; 
}

bool ProteinProbEstimator::getGridSearch()
{
  return dogridSearch;
}

bool ProteinProbEstimator::getOutputDecoys()
{
  return outputDecoys;
}

string ProteinProbEstimator::getProteinFN()
{
  return proteinFN;
}

bool ProteinProbEstimator::getTabDelimitedOutput()
{
  return tabDelimitedOut;
}