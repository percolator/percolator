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

/**
 * output protein PEPs
 */
const bool ProteinProbEstimator::outputPEPs = true;


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
  sort(out.begin(), out.end());
}

double antiderivativeAt(double m, double b, double xVal)
{
  return m*xVal*xVal/2.0 + b*xVal;
}

double area(double x1, double y1, double x2, double y2, double max_x)
{
  double m = (y2-y1)/(x2-x1);
  double b = y1-m*x1;

  double area =  antiderivativeAt(m, b, min(max_x, x2) ) - antiderivativeAt(m, b, x1);
  
  if(isnan(area)) return 0.0;
  else return area;
}

double squareAntiderivativeAt(double m, double b, double xVal)
{
  // turn into ux^2+vx+t
  double u = m*m;
  double v = 2*m*b;
  double t = b*b;

  return u*xVal*xVal*xVal/3.0 + v*xVal*xVal/2.0 + t*xVal;
}

ProteinProbEstimator::ProteinProbEstimator(double alpha_par, double beta_par, double gamma_par ,bool __tiesAsOneProtein
			 ,bool __usePi0, bool __outputEmpirQVal, bool __groupProteins, bool __noseparate, bool __noprune, bool __dogridSearch) {
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
}

ProteinProbEstimator::~ProteinProbEstimator(){
  delete proteinGraph;
}


bool ProteinProbEstimator::initialize(Scores* fullset){
  // percolator's peptide level statistics
  peptideScores = fullset;
  setTargetandDecoysNames();
}

void ProteinProbEstimator::run(){
  
  srand(time(NULL)); cout.precision(8); cerr.precision(8);
  // by default, a grid search is executed to estimate the values of the
  // parameters gamma alpha and beta
  if(dogridSearch) {
    if(VERB > 1) {
      cerr << "The parameters for the model will be estimated by grid search."
          << endl;
    }
    gridSearch();
  }
  else
  {
    if(alpha == -1)
      alpha = default_alpha;
    if(beta == -1)
      beta = default_beta;
    if(gamma == -1)
      gamma = default_gamma;
  }
  if(VERB > 1) {
      cerr << "\nThe following parameters have been chosen;\n";
      cerr << "gamma = " << gamma << endl;
      cerr << "alpha = " << alpha << endl;
      cerr << "beta  = " << beta << endl;
      cerr << "Protein level probabilities will now be calculated\n";
  }

  delete proteinGraph;
  proteinGraph = new GroupPowerBigraph (peptideScores,alpha,beta,gamma,groupProteins,noseparate,noprune);
  proteinGraph->getProteinProbs();
  pepProteins.clear();
  pepProteins = proteinGraph->getProteinProbsPercolator();
  
  estimateQValues();
  
  if(usePi0)
  {
    pi0 = *qvalues.rbegin();
    //estimatePValues();
    //pi0 = estimatePi0();
  }
  estimateQValuesEmp();
  updateProteinProbabilities();
  
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
	bool isdecoy = proteins[proteinName].getIsDecoy();
	combined.push_back(std::make_pair<double,bool>(prob,isdecoy));
      }
  }
  pvalues.clear();
  vector<pair<double, bool> >::const_iterator myPair = combined.begin();
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
  transform(pvalues.begin(), pvalues.end(), pvalues.begin(), bind2nd(divides<double> (),
                                                   (double)nDecoys));
}


double ProteinProbEstimator::estimatePi0(const unsigned int numBoot) 
{
  vector<double> pBoot, lambdas, pi0s, mse;
  vector<double>::iterator start;
  int numLambda = 100;
  double maxLambda = 0.5;
  size_t n = pvalues.size();
  // Calculate pi0 for different values for lambda
  // N.B. numLambda and maxLambda are global variables.
  for (unsigned int ix = 0; ix <= numLambda; ++ix) {
    double lambda = ((ix + 1) / (double)numLambda) * maxLambda;
    // Find the index of the first element in p that is < lambda.
    // N.B. Assumes p is sorted in ascending order.
    start = lower_bound(pvalues.begin(), pvalues.end(), lambda);
    // Calculates the difference in index between start and end
    double Wl = (double)distance(start, pvalues.end());
    double pi0 = Wl / n / (1 - lambda);
    if (pi0 > 0.0) {
      lambdas.push_back(lambda);
      pi0s.push_back(pi0);
    }
  }
  if(pi0s.size()==0){
    cerr << "Error in the input data: too good separation between target "
        << "and decoy PSMs.\nImpossible to estimate pi0. Terminating.\n";
    exit(0);
  }
  double minPi0 = *min_element(pi0s.begin(), pi0s.end());
  // Initialize the vector mse with zeroes.
  fill_n(back_inserter(mse), pi0s.size(), 0.0);
  // Examine which lambda level that is most stable under bootstrap
  for (unsigned int boot = 0; boot < numBoot; ++boot) {
    // Create an array of bootstrapped p-values, and sort in ascending order.
    bootstrap<double> (pvalues, pBoot);
    n = pBoot.size();
    for (unsigned int ix = 0; ix < lambdas.size(); ++ix) {
      start = lower_bound(pBoot.begin(), pBoot.end(), lambdas[ix]);
      double Wl = (double)distance(start, pBoot.end());
      double pi0Boot = Wl / n / (1 - lambdas[ix]);
      // Estimated mean-squared error.
      mse[ix] += (pi0Boot - minPi0) * (pi0Boot - minPi0);
    }
  }
  // Which index did the iterator get?
  unsigned int minIx = distance(mse.begin(), min_element(mse.begin(),
                                                         mse.end()));
  double pi0 = max(min(pi0s[minIx], 1.0), 0.0);
  return pi0;
}

unsigned ProteinProbEstimator::getQvaluesBelowLevel(double level)
{   
    unsigned nP = 0;
    for (std::map<const std::string,Protein>::const_iterator myP = proteins.begin(); 
	 myP != proteins.end(); ++myP) {
	 if(myP->second.getQ() < level && !myP->second.getIsDecoy()) nP++;
    }
    return nP;
}

unsigned ProteinProbEstimator::getQvaluesBelowLevelDecoy(double level)
{   
    unsigned nP = 0;
    for (std::map<const std::string,Protein>::const_iterator myP = proteins.begin(); 
	 myP != proteins.end(); ++myP) {
	 if(myP->second.getQ() < level && myP->second.getIsDecoy()) nP++;
    }
    return nP;
}


void ProteinProbEstimator::estimateQValues()
{

  int nP = 0;
  double sum = 0.0;
  double qvalue = 0.0;
  qvalues.clear();
  // assuming probabilites sorted in decending order
  for (std::multimap<double,std::vector<std::string> >::const_iterator it = pepProteins.begin(); 
       it != pepProteins.end(); it++) 
  {
    int ntargets = countTargets(it->second);
    sum += (double)(it->first * ntargets);
    nP += ntargets;
    qvalue = (sum / (double)nP);
    qvalues.push_back(qvalue);
  }
  std::partial_sum(qvalues.rbegin(),qvalues.rend(),qvalues.rbegin(),mymin);
}

void ProteinProbEstimator::estimateQValuesEmp()
{
    // assuming combined sorted in decending order
  double nDecoys = 0;
  double nTargets = 0;
  double qvalue = 0.0;
  double numDecoy = 0;
  pvalues.clear();
  qvaluesEmp.clear();
  for (std::multimap<double,std::vector<std::string> >::const_iterator it = pepProteins.begin(); 
       it != pepProteins.end(); it++) 
  {
    nTargets += countTargets(it->second);
    numDecoy = countDecoys(it->second);
    nDecoys += numDecoy;
    qvalue = (double)nDecoys / (double)nTargets;
    if(qvalue > 1.0) qvalue = 1.0;
    qvaluesEmp.push_back(qvalue);
    if(numDecoy > 0)
      pvalues.push_back((nDecoys)/(double)(numberDecoyProteins));
    else 
      pvalues.push_back((nDecoys+(double)1)/(numberDecoyProteins+(double)1));
  }

  double factor = pi0 * ((double)nTargets / (double)nDecoys);
  std::transform(qvaluesEmp.begin(), qvaluesEmp.end(), 
		 qvaluesEmp.begin(), bind2nd(multiplies<double> (),factor));
  std::partial_sum(qvaluesEmp.rbegin(), qvaluesEmp.rend(), qvaluesEmp.rbegin(), mymin);
}

void ProteinProbEstimator::updateProteinProbabilities()
{
  std::vector<double> peps;
  std::vector<std::vector<std::string> > proteinNames;
  transform(pepProteins.begin(), pepProteins.end(), back_inserter(peps), RetrieveKey());
  transform(pepProteins.begin(), pepProteins.end(), back_inserter(proteinNames), RetrieveValue());

  for (unsigned i = 0; i < peps.size(); i++) 
  {
    double pep = peps[i];
    std::vector<std::string> proteinlist = proteinNames[i];
    for(unsigned j = 0; j < proteinlist.size(); j++)
    {
      std::string proteinName = proteinlist[j];
      proteins[proteinName].setPEP(pep);
      proteins[proteinName].setQ(qvalues[i]);
      proteins[proteinName].setQemp(qvaluesEmp[i]);
      proteins[proteinName].setP(pvalues[i]);
    }
  }
}



std::map<const std::string,Protein> ProteinProbEstimator::getProteins()
{
  return this->proteins;
}


void ProteinProbEstimator::setTargetandDecoysNames()
{
  vector<ScoreHolder>::iterator psm = peptideScores->begin();
  
  for (; psm!= peptideScores->end(); ++psm) {
    
    set<string>::iterator protIt = psm->pPSM->proteinIds.begin();
    // for each protein
    for(; protIt != psm->pPSM->proteinIds.end(); protIt++){

      if(proteins.find(*protIt) == proteins.end())
      {
	Protein newprotein(*protIt,0.0,0.0,0.0,0.0,psm->isDecoy());
	newprotein.setPeptides(std::vector<std::string>(1,psm->pPSM->getPeptideSequence()));
	proteins.insert(std::make_pair<std::string,Protein>(*protIt,newprotein));
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
	proteins[*protIt].setPeptide(psm->pPSM->getPeptideSequence());
      }
    }
  }
  numberDecoyProteins = falsePosSet.size();
  numberTargetProteins = truePosSet.size();
}

void ProteinProbEstimator::gridSearch()
{

  double gamma_temp, alpha_temp, beta_temp;
  gamma_temp = alpha_temp = beta_temp = 0.01;

  GroupPowerBigraph gpb( peptideScores, default_gamma, default_alpha, default_beta,groupProteins,noprune,noseparate);

  double gamma_best, alpha_best, beta_best;
	 gamma_best = alpha_best = beta_best = -1.0;
  double best_objective = -100000000;

  double threshold = 0.01;
  double threshold2 = 0.05;
  double threshold3 = 0.1;
  int rocN = 75;
  
  //TODO make the range of the grid search and the N parametizable or according to data size
  
//   double gamma_search[] = {0.1, 0.25, 0.5, 075};
//   double beta_search[] = {0.0, 0.01, 0.15, 0.020, 0.025, 0.05};
//   double alpha_search[] = {0.01, 0.04, 0.09, 0.16, 0.25, 0.36};
  
//   double gamma_search[] = {0.1,0.25, 0.5, 075, 0.9};
//   double beta_search[] = {0.0, 0.01, 0.15, 0.025,0.35,0.05,0.1};
//   double alpha_search[] = {0.01, 0.04,0.09, 0.16, 0.25, 0.36,0.5};
  
  double gamma_search[] = {0.1, 0.5, 075};
  double beta_search[] = {0.0, 0.01, 0.15, 0.025, 0.05};
  double alpha_search[] = {0.01, 0.04, 0.16, 0.25, 0.36};
  
  for (unsigned int i=0; i<sizeof(gamma_search)/sizeof(double); i++)
  {
    for (unsigned int j=0; j<sizeof(alpha_search)/sizeof(double); j++)
    {
      for (unsigned int k=0; k<sizeof(beta_search)/sizeof(double); k++)
      {
	gamma = gamma_search[i];
	alpha = alpha_search[j];
	beta = beta_search[k];
	//std::cout << "Grid searching : " << alpha << " " << beta << " " << gamma << std::endl;
	gpb.setAlphaBetaGamma(alpha, beta, gamma);
	gpb.getProteinProbs();
	pair< vector< vector< string > >, std::vector< double > > NameProbs;
	NameProbs = gpb.getProteinProbsAndNames();
	//pepProteins.clear();
	//pepProteins = gpb.getProteinProbsPercolator();
	std::vector<double> prot_probs = NameProbs.second;
	std::vector<std::vector<std::string> > prot_names = NameProbs.first;
	//transform(pepProteins.begin(), pepProteins.end(), back_inserter(prot_probs), RetrieveKey());
	//transform(pepProteins.begin(), pepProteins.end(), back_inserter(prot_names), RetrieveValue());
	std::pair<std::vector<double>,std::vector<double> > EstEmp = getEstimated_and_Empirical_FDR(prot_names,prot_probs);
	pair<std::vector<int>, std::vector<int> > roc = getROC(prot_names);
	double rocR = getROC_N(roc.first, roc.second, rocN);
	double fdr_mse = getFDR_divergence(EstEmp.first, EstEmp.second, threshold3);
	double lambda = 0.15;
	double current_objective = lambda * rocR - (1-lambda) * fdr_mse;
	if (current_objective > best_objective)
	{
	  best_objective = current_objective;
	  gamma_best = gamma;
	  alpha_best = alpha;
	  beta_best = beta;
	}
	//cerr << gamma << " " << alpha << " " << beta << " : " << roc50 << " " << fdr_mse << " " << current_objective << endl;
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
  
  
  std::vector<std::pair<std::string,Protein> > myvec(proteins.begin(), proteins.end());
  std::sort(myvec.begin(), myvec.end(), IntCmpProb());

  ofstream os;
  os.open(xmlOutputFN.data(), ios::app);
  // append PROTEINs tag
  os << "  <proteins>" << endl;
  for (std::vector<std::pair<std::string,Protein> > ::const_iterator myP = myvec.begin(); 
	 myP != myvec.end(); myP++) {

        os << "    <protein p:protein_id=\"" << myP->second.getName() << "\"";
  
        if (Scores::isOutXmlDecoys()) {
          if(myP->second.getIsDecoy()) os << " p:decoy=\"true\"";
          else  os << " p:decoy=\"false\"";
        }
        os << ">" << endl;
        if(ProteinProbEstimator::outputPEPs)
          os << "      <pep>" << myP->second.getPEP() << "</pep>" << endl;
        if(ProteinProbEstimator::getOutputEmpirQval())
          os << "      <q_value_emp>" << myP->second.getQemp() << "</q_value_emp>\n";
        os << "      <q_value>" << myP->second.getQ() << "</q_value>\n";
	os << "      <p_value>" << myP->second.getP() << "</p_value>\n";
	std::vector<std::string> peptides = myP->second.getPeptides();
	for(std::vector<std::string>::const_iterator peptIt = peptides.begin(); peptIt != peptides.end(); peptIt++)
	{
	  if(!peptIt->empty())
	  {
	    os << "      <peptide_seq seq=\"" << peptIt->c_str() << "\"/>"<<endl;
	  }
	    
	}
        os << "    </protein>" << endl;
  }
    
  os << "  </proteins>" << endl << endl;
  os.close();
}


string ProteinProbEstimator::printCopyright(){
  ostringstream oss;
  oss << "Copyright (c) 2008-9 University of Washington. All rights reserved.\n"
      << "Written by Oliver R. Serang (orserang@u.washington.edu) in the\n"
      << "Department of Genome Sciences at the University of Washington.\n";
  return oss.str();
}


double ProteinProbEstimator::getROC_N(const std::vector<int> & fpArray, const std::vector<int> & tpArray, int N)
{
  double rocN = 0.0;

  if ( fpArray.back() < N )
    {
      cerr << "There are not enough false positives; needed " << N << " and was only given " << fpArray.back() << endl << endl;
      exit(1);
    }

  for (int k=0; k<fpArray.size()-1; k++)
    {
      // find segments where the fp value changes
	  
      if ( fpArray[k] >= N )
	break;

      if ( fpArray[k] != fpArray[k+1] )
	{
	  // this line segment is a function
	      
	  double currentArea = area(fpArray[k], tpArray[k], fpArray[k+1], tpArray[k+1], N);
	  rocN += currentArea;
	}
    }
  return rocN / (N * tpArray.back());
}

pair<std::vector<double>, std::vector<double> > ProteinProbEstimator::getEstimated_and_Empirical_FDR(std::vector<std::vector<string> > names, 
								   std::vector<double> probabilities)
{
  std::vector<double> estFDR_array, empFDR_array;
  double fpCount = 0.0, tpCount = 0.0;
  double totalFDR = 0.0, estFDR = 0.0, empFDR = 0.0;
  for (int k=0; k<names.size(); k++)
    {
      double prob = probabilities[k];
      int fpChange = countDecoys(names[k]);
      int tpChange = countTargets(names[k]);
      
      if ( tpChange > 0 && tiesAsOneProtein)
	tpChange = 1;
      if ( fpChange > 0 && tiesAsOneProtein )
	fpChange = 1;

      fpCount += (double)fpChange;
      tpCount += (double)tpChange;
      
      totalFDR += (prob) * (double)(tpChange);
      estFDR = totalFDR / (tpCount);
      empFDR = fpCount / tpCount; 
      if(empFDR > 1.0) empFDR = 1.0;
      if(estFDR > 1.0) estFDR = 1.0;
      
      estFDR_array.push_back(estFDR);
      //empFDR_array.push_back(totalFDR);
      empFDR_array.push_back(empFDR);

    }
  std::partial_sum(estFDR_array.rbegin(),estFDR_array.rend(),estFDR_array.rbegin(),mymin);
  double factor = pi0 * ((double)tpCount / (double)fpCount);
  std::transform(empFDR_array.begin(), empFDR_array.end(), 
		 empFDR_array.begin(), bind2nd(multiplies<double> (),factor));
  std::partial_sum(empFDR_array.rbegin(),empFDR_array.rend(),empFDR_array.rbegin(),mymin);
  return pair<std::vector<double>, std::vector<double> >(estFDR_array, empFDR_array);
}

std::vector<double> diffVector(const std::vector<double> &a, 
			       const std::vector<double> &b)
{
  std::vector<double> result = std::vector<double>(a.size());

  for (int k=0; k<result.size(); k++)
    {
      result[k] = a[k] - b[k];
    }
  
  return result;
}

double ProteinProbEstimator::getFDR_divergence(const std::vector<double> estFDR, const std::vector<double> empFDR, double THRESH)
{
  std::vector<double> diff = diffVector(estFDR,empFDR);

  double tot = 0.0;

  int k;
  for (k=0; k<diff.size()-1; k++)
    {
      // stop if no part of the estFDR is < threshold
      if ( estFDR[k] >= THRESH )
	{
	  if ( k == 0 )
	    tot = 1.0 / 0.0;

	  break;
	}

      tot += area(estFDR[k], diff[k], estFDR[k+1], diff[k+1], estFDR[k+1]);
    }

  double xRange = min(THRESH, estFDR[k]) - estFDR[0];

  if ( isinf(tot) )
    return tot;

  return tot / xRange;
}

pair<std::vector<int>, std::vector<int> > ProteinProbEstimator::getROC(std::vector<std::vector<string> > names)
{
  std::vector<int> fps, tps;
  int fpCount, tpCount;
  fpCount = tpCount = 0;
  
  for (int k=0; k<names.size(); k++)
    {
      int fpChange = countDecoys(names[k]);
      int tpChange = countTargets(names[k]);
      
      if ( tpChange > 0 && tiesAsOneProtein )
	tpChange = 1;
      if ( fpChange > 0 && tiesAsOneProtein )
	fpChange = 1;
      
      fpCount += fpChange;
      tpCount += tpChange;
      
      fps.push_back( fpCount );
      tps.push_back( tpCount );
	
    }

  fps.push_back( fpCount );
  tps.push_back( tpCount );	  
  fps.push_back( falsePosSet.size() );
  tps.push_back( truePosSet.size() );
  
  return pair<std::vector<int>, std::vector<int> >(fps, tps);
}


void ProteinProbEstimator::setOutputEmpirQval(bool outputEmpirQVal)
{
  this->outputEmpirQVal = outputEmpirQVal;
}

void ProteinProbEstimator::setTiesAsOneProtein(bool tiesAsOneProtein)
{
  this->tiesAsOneProtein = tiesAsOneProtein;
}

void ProteinProbEstimator::setUsePio(bool usePi0)
{
  this->usePi0 = usePi0;
}

void ProteinProbEstimator::setGroupProteins(bool groupProteins)
{
  this->groupProteins = groupProteins;
}

void ProteinProbEstimator::setPruneProteins(bool noprune)
{
  this->noprune = noprune;
}

void ProteinProbEstimator::setSeparateProteins(bool noseparate)
{
  this->noseparate = noseparate;
}

bool ProteinProbEstimator::getOutputEmpirQval()
{
  return this->outputEmpirQVal;
}

bool ProteinProbEstimator::getTiesAsOneProtein()
{
  return this->tiesAsOneProtein;
}

bool ProteinProbEstimator::getUsePio()
{
  return this->usePi0;
}

double ProteinProbEstimator::getPi0()
{
  return this->pi0;
}

bool ProteinProbEstimator::getGroupProteins()
{
  return this->groupProteins;
}


bool ProteinProbEstimator::getPruneProteins()
{
  return this->noprune;
}

bool ProteinProbEstimator::getSeparateProteins()
{
  return this->noseparate;
}


int ProteinProbEstimator::countTargets(std::vector<std::string> proteinList)
{
  int count = 0;
  for(int i = 0; i < proteinList.size(); i++)
  {
    if(truePosSet.count(proteinList[i]) > 0)
    {
      count++;
    }
  }
  if(tiesAsOneProtein && count > 0)
    return 1;
  else
    return count;
}

int ProteinProbEstimator::countDecoys(std::vector<std::string> proteinList)
{
  int count = 0;
  for(int i = 0; i < proteinList.size(); i++)
  {
    if(falsePosSet.count(proteinList[i]) > 0)
    {
      count++;
    }
  }
  
  if(tiesAsOneProtein && count > 0)
    return 1;
  else
    return count;
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


