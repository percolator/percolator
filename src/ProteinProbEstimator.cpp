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

double antiderivativeAt(double m, double b, double xVal)
{
  return m*xVal*xVal/2.0 + b*xVal;
}

double area(double x1, double y1, double x2, double y2, double max_x)
{
  double m = (y2-y1)/(x2-x1);
  double b = y1-m*x1;

  return antiderivativeAt(m, b, min(max_x, x2) ) - antiderivativeAt(m, b, x1);
}

double squareAntiderivativeAt(double m, double b, double xVal)
{
  // turn into ux^2+vx+t
  double u = m*m;
  double v = 2*m*b;
  double t = b*b;

  //  cout << "\t\tsquareAntiderivativeAt " << xVal << " is " << u*xVal*xVal*xVal/3.0 + v*xVal*xVal/2.0 + t*xVal << endl;
  //  cout << "u, v, t = " << u << " " << v << " " << t << endl;

  return u*xVal*xVal*xVal/3.0 + v*xVal*xVal/2.0 + t*xVal;
}

ProteinProbEstimator::ProteinProbEstimator(double alpha_par, double beta_par, double gamma_par ,bool tiesAsOneProtein
			 ,bool usePi0, bool outputEmpirQVal, bool groupProteins, bool noseparate, bool noprune) {
  peptideScores = 0;
  proteinGraph = 0;
  gamma = gamma_par;
  alpha = alpha_par;
  beta = beta_par;
  numberDecoyProteins = 0;
  numberTargetProteins = 0;
  pi0 = 0.0;
  this->tiesAsOneProtein = tiesAsOneProtein;
  this->usePi0 = usePi0;
  this->outputEmpirQVal = outputEmpirQVal;
  this->groupProteins = groupProteins;
  this->noseparate = noseparate;
  this->noprune = noprune;
}

ProteinProbEstimator::~ProteinProbEstimator(){
  delete proteinGraph;
}


bool ProteinProbEstimator::initialize(Scores* fullset){
  // percolator's peptide level statistics
  peptideScores = fullset;
  setTargetandDecoysNames();
}

void ProteinProbEstimator::run(bool startGridSearch){
  
  srand(time(NULL)); cout.precision(8); cerr.precision(8);
  // by default, a grid search is executed to estimate the values of the
  // parameters gamma alpha and beta
  if(startGridSearch) {
    if(VERB > 1) {
      cerr << "The parameters for the model will be estimated by grid search."
          << endl;
    }
    gridSearch();
    if(VERB > 1) {
      cerr << "\nThe following parameters have been chosen;\n";
      cerr << "gamma = " << gamma << endl;
      cerr << "alpha = " << alpha << endl;
      cerr << "beta  = " << beta << endl;
      cerr << "Protein level probabilities will now be calculated\n";
    }
  }
  else
  {
    alpha = default_alpha;
    beta = default_beta;
    gamma = default_gamma;
  }

  //GroupPowerBigraph::LOG_MAX_ALLOWED_CONFIGURATIONS = ;
  
  delete proteinGraph;
  //proteinGraph = new GroupPowerBigraph (peptideScores,alpha,beta,gamma,groupProteins,noprune,noseparate);
  proteinGraph = new GroupPowerBigraph (peptideScores,alpha,beta,gamma);
  
//   std::cerr << "Size of peptide Scores : " << peptideScores->size() << std::endl;
  proteinGraph->getProteinProbs();
  estimateQValues();
  if(usePi0)
    /** call get Pvalues and calculate Pi0 as qvality **/
    //TODO include the ties as one protein functionality and the use Pi0 as well
    estimateQValuesEmp(estQ.rbegin()->second);
  else
    estimateQValuesEmp(1.0);
  updateProteinProbabilities();
//   std::cerr << "Size of protein Prob Graph : " << proteinGraph->probabilityR.size() << std::endl;
}


unsigned int ProteinProbEstimator::getQvaluesBelowLevel(double level)
{   
    unsigned nP = 0;
    for (std::map<std::string,Protein>::const_iterator myP = proteins.begin(); 
	 myP != proteins.end(); ++myP) {
	 if(myP->second.getQ() <= level) nP++;
    }
    return nP;
}

void ProteinProbEstimator::estimateQValues()
{
  int nP = 1;
  double sum = 0.0;
  double qvalue = 0.0;
  // assuming combined sorted in decending order
  Array<double> sorted = proteinGraph->getDescendingProteinsAndWeights().second;
  Array< Array<string> > protein_ids = proteinGraph->getDescendingProteinsAndWeights().first;
  for (int i=0; i<sorted.size(); i++,nP++)
  {
    double pep = 1.0 - sorted[i];
    if(pep < 0.0 || !isnan(pep)) pep = 0.0;
    if(pep > 1.0 || isinf(pep))  pep = 1.0;
    sum += pep;
    qvalue = (sum / (double)nP);
    for (int j=0; j<protein_ids[i].size(); j++)
      estQ.insert(std::pair<std::string,double>(std::string(protein_ids[i][j]),qvalue));
  }
  partial_sum(estQ);
}

void ProteinProbEstimator::estimateQValuesEmp(double pi0)
{
    // assuming combined sorted in decending order
  double nDecoys = 0;
  double nTargets = 0;
  Array<double> pps = proteinGraph->getDescendingProteinsAndWeights().second;
  Array< Array<string> > protein_ids = proteinGraph->getDescendingProteinsAndWeights().first;
  for (int i=0; i<pps.size(); i++) {
    double pep = 1.0 - pps[i];
    if(pep < 0.0 || !isnan(pep)) pep = 0.0;
    for (int j=0; j<protein_ids[i].size(); j++) {
      if (proteins[protein_ids[i][j]].getIsDecoy()) {
	++nDecoys;
      } else {
	++nTargets;
	empQ.insert(std::make_pair<std::string,double>(std::string(protein_ids[i][j]),(((double)nDecoys) / (double)nTargets)));
	
      }
    }
  }
  
  double factor = pi0 * ((double)nTargets / (double)nDecoys);
  std::for_each(empQ.begin(), empQ.end(), mul_x(factor));
  partial_sum(empQ);

}

void ProteinProbEstimator::updateProteinProbabilities()
{
  Array<double> pps = proteinGraph->getDescendingProteinsAndWeights().second;
  Array< Array<string> > protein_ids = proteinGraph->getDescendingProteinsAndWeights().first;
  for (int i=0; i<pps.size(); i++) {
    double pep = 1.0 - pps[i];
    for (int j=0; j<protein_ids[i].size(); j++) {
      proteins[protein_ids[i][j]].setPEP(pep);
      proteins[protein_ids[i][j]].setQ(estQ[protein_ids[i][j]]);
      proteins[protein_ids[i][j]].setQemp(empQ[protein_ids[i][j]]);
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
      // CHECK PROTEIN EXITS 
      Protein newprotein(*protIt,0.0,0.0,0.0,0.0,psm->label != -1 ? false : true);
      newprotein.setPeptide(*psm);
      proteins.insert(std::make_pair<std::string,Protein>(*protIt,newprotein));
   
      if(std::string(*protIt).find("random") != string::npos) //is decoy
      {
	//std::cout << "Reading Decoy protein : " << *protIt << std::endl;
	falsePosSet.insert(*protIt);
	numberDecoyProteins++;
      }
      else
      {
	//std::cout << "Reading Target protein : " << *protIt << std::endl;
	truePosSet.insert(*protIt);
	numberTargetProteins++;

      }
	
    }
  }
}

void ProteinProbEstimator::gridSearch()
{

  double gamma_temp, alpha_temp, beta_temp;
  gamma_temp = alpha_temp = beta_temp = 0.01;

  //GroupPowerBigraph gpb( peptideScores, default_gamma, default_alpha, default_beta,groupProteins,noprune,noseparate);
  GroupPowerBigraph gpb( peptideScores, default_gamma, default_alpha, default_beta);
  double gamma_best, alpha_best, beta_best;
	 gamma_best = alpha_best = beta_best = -1.0;
  double best_objective = -100000000;

  double gamma_search[] = {0.1, 0.5, 0.9};
  double alpha_search[] = {0.01, 0.04, 0.09, 0.16, 0.25, 0.36};
  double beta_search[] = {0.0, 0.01, 0.025, 0.05};
  
  for (unsigned int i=0; i<sizeof(gamma_search)/sizeof(double); i++)
  {
    for (unsigned int j=0; j<sizeof(alpha_search)/sizeof(double); j++)
    {
      for (unsigned int k=0; k<sizeof(beta_search)/sizeof(double); k++)
      {
	gamma = gamma_search[i];
	alpha = alpha_search[j];
	beta = beta_search[k];
	std::cout << "Grid searching : " << alpha << " " << beta << " " << gamma << std::endl;
	gpb.setAlphaBetaGamma(alpha, beta, gamma);
	gpb.getProteinProbs();
	pair<Array<Array<string> >, Array<double> > prot_groups_and_probs = gpb.getDescendingProteinsAndWeights();
	Array<Array<string> > prot_names = prot_groups_and_probs.first;
	Array<double> prot_probs = prot_groups_and_probs.second;
	pair<Array<int>, Array<int> > roc = getROC(prot_names, prot_probs, falsePosSet, truePosSet);
	double roc50 = getROC_N(roc.first, roc.second, 50);
	pair<Array<double>, Array<double> > est_and_emp_fdr = getEstimated_and_Empirical_FDR(prot_names, prot_probs, falsePosSet, truePosSet);
	double fdr_mse = getFDR_divergence(est_and_emp_fdr.first, est_and_emp_fdr.second, 0.10);
	double lambda = 0.15;
	double current_objective = lambda * roc50 - (1-lambda) * fdr_mse;
	if (current_objective > best_objective)
	{
	  best_objective = current_objective;
	  gamma_best = gamma;
	  alpha_best = alpha;
	  beta_best = beta;
	}
	cerr << gamma << " " << alpha << " " << beta << " : " << roc50 << " " << fdr_mse << " " << current_objective << endl;
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
  
  ofstream os;
  os.open(xmlOutputFN.data(), ios::app);
  // append PROTEINs tag
  os << "  <proteins>" << endl;
  for (std::map<std::string,Protein>::const_iterator myP = proteins.begin(); 
	 myP != proteins.end(); ++myP) {

        os << "    <protein p:protein_id=\"" << myP->second.getName() << "\"";
  
        if (Scores::isOutXmlDecoys()) {
          if(myP->second.getIsDecoy()) os << " p:decoy=\"true\"";
          else  os << " p:decoy=\"false\"";
        }
        os << ">" << endl;
        if(ProteinProbEstimator::outputPEPs)
          os << "      <pep>" << myP->second.getPEP() << "</pep>" << endl;
        if(ProteinProbEstimator::getOutputEmpirQval())
          os << "      <q_value>" << myP->second.getQemp() << "</q_value>\n";
        else
          os << "      <q_value>" << myP->second.getQ() << "</q_value>\n";
	for(std::vector<ScoreHolderPeptide>::const_iterator peptIt = myP->second.getPeptide().begin(); 
	    peptIt != myP->second.getPeptide().end(); peptIt++){
	  string pept = (*peptIt).pPSM->getPeptideSequence();
	  os << " 	<peptide_seq seq=\"" << pept << "\"/>"<<endl;
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

int ProteinProbEstimator::matchCount( const set<string> & positiveNames, const Array<string> & atThreshold )
{
  int count = 0;
  
  for (int k=0; k<atThreshold.size(); k++)
    {
      if ( positiveNames.count( atThreshold[k] ) > 0 )
	count++;
    }

  return count;
}

double ProteinProbEstimator::getROC_N(const Array<int> & fpArray, const Array<int> & tpArray, int N)
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

pair<Array<double>, Array<double> > ProteinProbEstimator::getEstimated_and_Empirical_FDR(Array<Array<string> > names, 
								   Array<double> probabilities, 
								   const set<string> & falsePosSet, 
								   const set<string> & truePosSet)
{
  Array<double> estFDR_array, empFDR_array;
  estFDR_array.add(0);
  empFDR_array.add(0);
  
  int fpCount, tpCount;
  fpCount = tpCount = 0;
  double totalFDR = 0.0, estFDR = 0.0;

  bool scheduledUpdate = false;
  
  double lastProb = -1.0;
  for (int k=0; k<names.size(); k++)
    {
      double prob = probabilities[k];
      int fpChange = matchCount(falsePosSet, names[k]);
      int tpChange = matchCount(truePosSet, names[k]);
      
      // for different style of grading, counting groups as a
      // single protein and throwing away any groups that include
      // TPs and FPs
      if ( tpChange > 0 && fpChange > 0 )
	tpChange = fpChange = 0;
      
      if ( tpChange > 0 )
	tpChange = 1;
      if ( fpChange > 0 )
	fpChange = 1;
      
      if ( prob != lastProb && lastProb != -1 )
	{
	  scheduledUpdate = true;
	}

      if ( scheduledUpdate )
	{
	  scheduledUpdate = false;

	  totalFDR += (1-prob) * (fpChange + tpChange);
	  estFDR = totalFDR / (fpCount + tpCount);
	  
	  estFDR_array.add(estFDR);
	  empFDR_array.add(totalFDR);
	}

      fpCount += fpChange;
      tpCount += tpChange;

      lastProb = prob;
    }
  
  return pair<Array<double>, Array<double> >(estFDR_array, empFDR_array);
}

double ProteinProbEstimator::getFDR_divergence(const Array<double> estFDR, const Array<double> empFDR, double THRESH)
{
  Vector diff = Vector(estFDR) - Vector(empFDR);

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

pair<Array<int>, Array<int> > ProteinProbEstimator::getROC(Array<Array<string> > names, Array<double> probabilities, 
							  const set<string> & falsePosSet, const set<string> & truePosSet)
{
  Array<int> fps, tps;
  fps.add(0);
  tps.add(0);
  
  int fpCount, tpCount;
  fpCount = tpCount = 0;

  bool scheduledUpdate = false;
  
  double lastProb = -1.0;
  for (int k=0; k<names.size(); k++)
    {
      double prob = probabilities[k];
      int fpChange = matchCount(falsePosSet, names[k]);
      int tpChange = matchCount(truePosSet, names[k]);
      
      // for different style of grading, counting groups as a
      // single protein and throwing away any groups that include
      // TPs and FPs
      if ( tpChange > 0 && fpChange > 0 )
	tpChange = fpChange = 0;
      
      if ( tpChange > 0 )
	tpChange = 1;
      if ( fpChange > 0 )
	fpChange = 1;
      
      if ( prob != lastProb && lastProb != -1 )
	{
	  scheduledUpdate = true;
	}

      if ( scheduledUpdate )
	{
	  fps.add( fpCount );
	  tps.add( tpCount );
	  scheduledUpdate = false;

	  //	  cout << fpCount << " " << tpCount << endl;

	  //	  totalFDR += (1-prob) * (fpChange + tpChange);
	  //	  estFDR = totalFDR / (fpCount + tpCount);
	}

      fpCount += fpChange;
      tpCount += tpChange;

      lastProb = prob;
    }

  fps.add( fpCount );
  tps.add( tpCount );	  
  
  fps.add( falsePosSet.size() );
  tps.add( truePosSet.size() );
  
  return pair<Array<int>, Array<int> >(fps, tps);
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





