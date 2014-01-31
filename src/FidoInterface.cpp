/*******************************************************************************
 Copyright 2006-2012 Lukas Käll <lukas.kall@scilifelab.se>

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

#include "FidoInterface.h"

const double FidoInterface::psmThreshold = 0.0;
const double FidoInterface::reduced_psmThreshold = 0.2;
const double FidoInterface::peptideThreshold = 0.001;
const double FidoInterface::reduced_peptideThreshold = 0.2;
const double FidoInterface::proteinThreshold = 0.01; 
const double FidoInterface::reduced_proteinThreshold = 0.2;
const double FidoInterface::peptidePrior = 0.1; 
const double FidoInterface::max_allow_configurations = 18;
const double FidoInterface::lambda = 0.15;

double trapezoid_area(double x1, double x2, double y1, double y2)
{
  double base = abs(x1 - x2);
  double height_avg = abs((y1 + y2) / 2);
  return base * height_avg;
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

double area(double x1, double y1, double x2, double y2)
{
  double m = (y2-y1)/(x2-x1);
  double b = y1-m*x1;
  double area =  antiderivativeAt(m, b, x2) - antiderivativeAt(m, b, x1);
  return area;
}

double areaSq(double x1, double y1, double x2, double y2) 
{
  double m = (y2-y1)/(x2-x1);
  double b = y1-m*x1;
  double area = squareAntiderivativeAt(m, b, x2) - squareAntiderivativeAt(m, b, x1);
  return area;
}

double Round(double dbVal, int nPlaces /* = 0 */)
{
    const double dbShift = pow(10.0, nPlaces);
    return  floor(dbVal * dbShift + 0.5) / dbShift; 
}
    
// get the number of decimal places
int GetDecimalPlaces(double dbVal)
{
    static const int MAX_DP = 10;
    static const double THRES = pow(0.1, MAX_DP);
    if (dbVal == 0.0)
        return 0;
    int nDecimal = 0;
    while (dbVal - floor(dbVal) > THRES && nDecimal < MAX_DP)
    {
        dbVal *= 10.0;
        nDecimal++;
    }
    return nDecimal;
}


FidoInterface::FidoInterface(double __alpha,double __beta,double __gamma,bool __nogroupProteins, 
			      bool __noseparate, bool __noprune, unsigned __depth,bool __reduceTree, 
			      bool __truncate, double mse_threshold,bool tiesAsOneProtein, bool usePi0, 
			      bool outputEmpirQVal, std::string decoyPattern,bool __trivialGrouping)
			      :ProteinProbEstimator(tiesAsOneProtein,usePi0,outputEmpirQVal,decoyPattern)
{
  alpha = __alpha;
  beta = __beta;
  gamma = __gamma;
  trivialGrouping = __trivialGrouping;
  nogroupProteins = __nogroupProteins;
  noseparate = __noseparate;
  noprune = __noprune;
  depth = __depth;
  reduceTree = __reduceTree;
  truncate = __truncate;
  threshold = mse_threshold;
  dogridSearch = false;
}

FidoInterface::~FidoInterface()
{  if(proteinGraph)
  {
    delete proteinGraph;
  }
  proteinGraph = 0;

}

void FidoInterface::run()
{
  dogridSearch = !(alpha != -1 && beta != -1 && gamma != -1);
  
  double peptidePrior_local = peptidePrior;
  if(computePriors)
  {
    peptidePrior_local = estimatePriors();
  
    if(VERB > 1)
    {
      std::cerr << "The estimated peptide level prior probability is : " << peptidePrior_local << std::endl;
    }
  }
  
  double local_protein_threshold = proteinThreshold;
  if(truncate) local_protein_threshold = 0.0;
  
  proteinGraph = new GroupPowerBigraph (alpha,beta,gamma,nogroupProteins,noseparate,noprune,trivialGrouping);
  proteinGraph->setMaxAllowedConfigurations(max_allow_configurations);
  proteinGraph->setPeptidePrior(peptidePrior_local);
  
  if(reduceTree && dogridSearch)
  {
    //NOTE lets create a smaller tree to estimate the parameters faster
    if(VERB > 1)
    {
      std::cerr << "Reducing the tree of proteins to increase the speed of the grid search.." << std::endl;
    }
    proteinGraph->setProteinThreshold(reduced_proteinThreshold);
    proteinGraph->setPsmThreshold(reduced_psmThreshold);
    proteinGraph->setPeptideThreshold(reduced_peptideThreshold);
    proteinGraph->setGroupProteins(false);
    proteinGraph->setSeparateProteins(false);
    proteinGraph->setPruneProteins(false);
    proteinGraph->setTrivialGrouping(true);
    proteinGraph->setMultipleLabeledPeptides(false);
  }
  else
  {
    proteinGraph->setProteinThreshold(local_protein_threshold);
    proteinGraph->setPsmThreshold(psmThreshold);
    proteinGraph->setPeptideThreshold(peptideThreshold);
    proteinGraph->setMultipleLabeledPeptides(allow_multiple_labeled_peptides);
  }
 
}

void FidoInterface::computeProbabilities()
{
  
  proteinGraph->read(peptideScores);
  
  time_t startTime;
  clock_t startClock;
  time(&startTime);
  startClock = clock();
  
  if(mayufdr)
  {
    computeFDR();
  }
  
  if(dogridSearch) 
  {
    if(VERB > 1) 
    {
      std::cerr << "The parameters for the model will be estimated by grid search.\n" << std::endl;
    }
    
    if(optimize)
      gridSearchOptimize(); 
    else
      gridSearch();
    
    time_t procStart;
    clock_t procStartClock = clock();
    time(&procStart);
    double diff = difftime(procStart, startTime);
    if (VERB > 1) cerr << "Estimating the parameters took : "
      << ((double)(procStartClock - startClock)) / (double)CLOCKS_PER_SEC
      << " cpu seconds or " << diff << " seconds wall time" << endl;
  }

  if(VERB > 1) 
  {
      cerr << "The following parameters have been chosen:\n";
      std::cerr.precision(10);
      cerr << "gamma = " << gamma << endl;
      cerr << "alpha = " << alpha << endl;
      cerr << "beta  = " << beta << endl;
      std::cerr.unsetf(std::ios::floatfield);
      cerr << "\nProtein level probabilities will now be estimated";
  }


  if(dogridSearch && reduceTree)
  {
    //NOTE lets create the tree again with all the members
    double local_protein_threshold = proteinThreshold;
    if(truncate) local_protein_threshold = 0.0;
    proteinGraph->setProteinThreshold(local_protein_threshold);
    proteinGraph->setPsmThreshold(psmThreshold);
    proteinGraph->setPeptideThreshold(peptideThreshold);
    proteinGraph->setGroupProteins(nogroupProteins);
    proteinGraph->setSeparateProteins(noseparate);
    proteinGraph->setPruneProteins(noprune);
    proteinGraph->setTrivialGrouping(trivialGrouping);
    proteinGraph->setMultipleLabeledPeptides(allow_multiple_labeled_peptides);
    proteinGraph->read(peptideScores);
  }
  
  proteinGraph->setAlphaBetaGamma(alpha,beta,gamma);
  proteinGraph->getProteinProbs();
  pepProteins.clear();
  proteinGraph->getProteinProbsPercolator(pepProteins);
}

//NOTE almost entirely duplicated of computeProbabilities, it could be refactored
void FidoInterface::computeProbabilitiesFromFile(ifstream &fin)
{
  
  proteinGraph->read(fin);
  
  time_t startTime;
  clock_t startClock;
  time(&startTime);
  startClock = clock();
  
  if(mayufdr)
  {
    computeFDR();
  }
  
  if(dogridSearch) 
  {
    if(VERB > 1) 
    {
      std::cerr << "The parameters for the model will be estimated by grid search.\n" << std::endl;
    }
    
    if(optimize)
      gridSearchOptimize(); 
    else
      gridSearch();
    
    time_t procStart;
    clock_t procStartClock = clock();
    time(&procStart);
    double diff = difftime(procStart, startTime);
    if (VERB > 1) cerr << "Estimating the parameters took : "
      << ((double)(procStartClock - startClock)) / (double)CLOCKS_PER_SEC
      << " cpu seconds or " << diff << " seconds wall time" << endl;
  }

  if(VERB > 1) 
  {
      cerr << "The following parameters have been chosen:\n";
      std::cerr.precision(10);
      cerr << "gamma = " << gamma << endl;
      cerr << "alpha = " << alpha << endl;
      cerr << "beta  = " << beta << endl;
      std::cerr.unsetf(std::ios::floatfield);
      cerr << "\nProtein level probabilities will now be estimated";
  }


  if(dogridSearch && reduceTree)
  {
    //NOTE lets create the tree again with all the members
    double local_protein_threshold = proteinThreshold;
    if(truncate) local_protein_threshold = 0.0;
    proteinGraph->setProteinThreshold(local_protein_threshold);
    proteinGraph->setPsmThreshold(psmThreshold);
    proteinGraph->setPeptideThreshold(peptideThreshold);
    proteinGraph->setGroupProteins(nogroupProteins);
    proteinGraph->setSeparateProteins(noseparate);
    proteinGraph->setPruneProteins(noprune);
    proteinGraph->setMultipleLabeledPeptides(allow_multiple_labeled_peptides);
    proteinGraph->read(fin);
  }
  
  proteinGraph->setAlphaBetaGamma(alpha,beta,gamma);
  proteinGraph->getProteinProbs();
  pepProteins.clear();
  proteinGraph->getProteinProbsPercolator(pepProteins);
}

void FidoInterface::gridSearch()
{
  double gamma_best, alpha_best, beta_best;
  gamma_best = alpha_best = beta_best = -1.0;
  double best_objective = -100000000;
  std::vector<std::vector<std::string> > names;
  std::vector<double> probs,empq,estq;
  std::vector<long double> gamma_search,beta_search,alpha_search;

  double roc, mse,current_objective;
  
  switch(depth)
  {
    case 0:    
      gamma_search = boost::assign::list_of(0.5);
      beta_search = boost::assign::list_of(0.001);
      alpha_search = boost::assign::list_of(0.008)(0.032)(0.128);
      break;
    
    case 1:
      gamma_search = boost::assign::list_of(0.1)(0.5)(0.9);
      beta_search = boost::assign::list_of(0.001);
      for (double k = 0.002; k <= 0.4; k*=4)
      {
       alpha_search.push_back(k);
      }
      break;
      
    case 2:
      gamma_search = boost::assign::list_of(0.1)(0.3)(0.5)(0.75)(0.9);
      beta_search = boost::assign::list_of(0.001);
      for (double k = 0.001; k <= 0.4; k*=2)
      {
       alpha_search.push_back(k);
      }
      break;
    
    default:
      gamma_search = boost::assign::list_of(0.5);
      beta_search = boost::assign::list_of(0.001);
      alpha_search = boost::assign::list_of(0.008)(0.032)(0.128);
  }

  if(alpha != -1)
    alpha_search = boost::assign::list_of(alpha);
  if(beta != -1)
    beta_search = boost::assign::list_of(beta);
  if(gamma != -1)
    gamma_search = boost::assign::list_of(gamma);
  
  for (unsigned int i = 0; i < gamma_search.size(); i++)
  {
    double gamma_local = gamma_search[i];
    
    for (unsigned int j = 0; j < alpha_search.size(); j++)
    {
      double alpha_local = alpha_search[j];
      
      for (unsigned int k = 0; k < beta_search.size(); k++)
      {

	double beta_local = beta_search[k];
	
	proteinGraph->setAlphaBetaGamma(alpha_local, beta_local, gamma_local);
	proteinGraph->getProteinProbs();
	proteinGraph->getProteinProbsAndNames(names,probs);
	getEstimated_and_Empirical_FDR(names,probs,empq,estq);
	getROC_AUC(names,probs,roc);
	getFDR_MSE(estq,empq,mse);
	
	current_objective = (lambda * roc) - fabs(((1-lambda) * (mse)));
	
	if(VERB > 2)
	{
	  std::cerr.precision(10);
	  std::cerr << "Grid searching Alpha= "  << alpha_local << " Beta= " << beta_local << " Gamma= "  << gamma_local << std::endl;
	  std::cerr.unsetf(std::ios::floatfield);
	  std::cerr << "The ROC AUC estimated values is : " << roc <<  std::endl;
	  std::cerr << "The MSE FDR estimated values is : " <<  mse << std::endl;
	  std::cerr << "Objective function with second roc and mse is : " << current_objective << std::endl;
	}  
	if (current_objective > best_objective)
	{
	  best_objective = current_objective;
	  gamma_best = gamma_local;
	  alpha_best = alpha_local;
	  beta_best = beta_local;
	}

      }
    }
  }
  
  alpha = alpha_best;
  beta = beta_best;
  gamma = gamma_best;
}

void FidoInterface::gridSearchOptimize()
{
 
  if(VERB > 1)
  {
    std::cerr << "Running super grid search..." << std::endl;
  }
  
  double gamma_best, alpha_best, beta_best;
  gamma_best = alpha_best = beta_best = -1.0;
  double best_objective = -100000000;
  std::vector<std::vector<std::string> > names;
  std::vector<double> probs,empq,estq; 
  double roc,mse,current_objective;
  
  double alpha_step = 0.05;
  double beta_step = 0.05;
  double gamma_step = 0.05;
  
  double beta_init = 0.00001;
  double alpha_init = 0.001;
  double gamma_init = 0.1;
  
  double gamma_limit = 0.5;
  double beta_limit = 0.05;
  double alpha_limit = 0.5;
  
  if(alpha != -1)
  {
    alpha_init = alpha_limit = alpha;
  }
  
  if(beta != -1)
  {
    beta_init = beta_limit = beta;
  }
  
  if(gamma != -1)
  {
    gamma_init = gamma_limit = gamma;
  }
  
  //NOTE very annoying the residue error of the floats that get acummulated in every iteration
  
  for (double i = gamma_init; i <= gamma_limit; i+=gamma_step)
  { 
    double gamma_local = i;
    
    for (double j = log10(beta_init); j <= Round(log10(beta_limit),2); j+=beta_step)
    {
      double original = pow(10,j);
      double beta_local = original - beta_init;
      if(beta_local > 0.0) beta_local = original;
      
      for (double k = log10(alpha_init); k <= Round(log10(alpha_limit),2); k+=alpha_step)
      {
       
	double alpha_local = pow(10,k);
	
	proteinGraph->setAlphaBetaGamma(alpha_local, beta_local, gamma_local);
	proteinGraph->getProteinProbs();
	proteinGraph->getProteinProbsAndNames(names,probs);
	getEstimated_and_Empirical_FDR(names,probs,empq,estq);
	getROC_AUC(names,probs,roc);
	getFDR_MSE(estq,empq,mse);
	
	current_objective = (lambda * roc) - fabs(((1-lambda) * (mse)));
	
	if(VERB > 2)
	{
	  std::cerr.precision(10);
	  std::cerr << "Grid searching Alpha= "  << alpha_local << " Beta= " << beta_local << " Gamma= "  << gamma_local << std::endl;
	  std::cerr.unsetf(std::ios::floatfield);
	  std::cerr << "The ROC AUC estimated values is : " << roc <<  std::endl;
	  std::cerr << "The MSE FDR estimated values is : " <<  mse << std::endl;
	  std::cerr << "Objective function with second roc and mse is : " << current_objective << std::endl;
	  
	}    
	if (current_objective > best_objective)
	{
	  best_objective = current_objective;
	  gamma_best = gamma_local;
	  alpha_best = alpha_local;
	  beta_best = beta_local;
	}
      }
    }
  }
  
  alpha = alpha_best;
  beta = beta_best;
  gamma = gamma_best;

}


void FidoInterface::getROC_AUC(const std::vector<std::vector<string> > &names,
					  const std::vector<double> &probabilities, double &auc)
{
  /* Estimate ROC auc1 area as : (So - no(no + 1) / 2) / (no*n1)
   * where no = number of target
   * where n1 = number of decoy
   * where So = SUM ri
   * where ri is the rank of i target in the ranked list of target and decoys
   */
  
  /* Estimate ROC auc2 area as : sum trapezoid area of each segment (integral of absolute value)
   * A_segment(i) = abs(X1-Xo) * abs((y1 + y2 ) / 2)
   * Where yo = number TP at segment i
   * Where y1 = number TP at segment i + 1
   * Where Xo = number FP at segment i
   * Where X1 = number FP at segment i + 1
   * Total Area = Total Area / total_TP * total_FP
   */
  
  /* Estimate ROC auc3 area as : sum trapezoid area with antiderivatives of each segment (absolute value of the integral)
   * A_segment(i) = ((yo - m*Xo)*X1 + m/2 * X1^2) - ((yo - m*Xo)*Xo - m/2 * X2^2))
   * Where yo = number TP at segment i
   * Where y1 = number TP at segment i + 1
   * Where Xo = number FP at segment i
   * Where X2 = number FP at segment i + 1
   * Where m = (y1 - y0) / (X1 - X0)
   * Total Area = abs(Total Area / total_TP * total_FP)
   */
  
  std::vector<bool> ranked_list; // true if is decoy
  std::vector<unsigned> tpArray,fpArray;
  
  unsigned prev_tp,prev_fp,tp,fp;
  prev_tp = prev_fp = tp = fp = 0;
  double prev_prob = -1;
  auc = 0.0;
  
  //assuming names and probabilities same size
  for (unsigned k=0; k < names.size() && fp <= rocN; k++)
  {
      double prob = probabilities[k];
      unsigned tpChange = countTargets(names[k]);
      unsigned fpChange = names[k].size() - tpChange;
      //if ties activated count groups as 1 protein
      if(tiesAsOneProtein)
      {
	if(tpChange) tpChange = 1;
	if(fpChange) fpChange = 1;
      }

      tp += tpChange;
      fp += fpChange;
      //should only do it when fp changes and either of them is != 0
      if(prev_prob != -1 && fp != 0 && tp != 0 && fp != prev_fp)
      {
	double trapezoid = trapezoid_area(fp,prev_fp,tp,prev_tp);
	prev_fp = fp;
	prev_tp = tp;
	auc += trapezoid;
      }   
      
      prev_prob = prob;
  }

  unsigned normalizer = (tp * fp);
  
  if(normalizer)
  {
    auc /= normalizer;
  }
  else
  {
    auc = 0.0;
  }
  
  return;
}


void FidoInterface::getEstimated_and_Empirical_FDR(const std::vector<std::vector<string> > &names,
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
  
  if(updateRocN) rocN = 50;
  
  //NOTE no need to store more q values since they will not be taken into account while estimating MSE FDR divergence
  for (unsigned int k=0; (k<names.size() && (estFDR <= threshold)); k++)
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
	  totalFDR += (prob) * (double)(tpChange + fpChange);
	  estFDR = totalFDR / (tpCount + fpCount);
	}
	else
	{
	  totalFDR += (prob) * (double)(tpChange);
	  estFDR = totalFDR / (tpCount);  
	}

	if(tpCount) empFDR = (fpCount * pi0 * TargetDecoyRatio) / tpCount; 
	
	if(empFDR > 1.0 || std::isnan(empFDR) || std::isinf(empFDR)) empFDR = 1.0;
	if(estFDR > 1.0 || std::isnan(estFDR) || std::isinf(estFDR)) estFDR = 1.0;
	    
	if(estFDR < previousEstQ) estFDR = previousEstQ;
	else previousEstQ = estFDR;
	    
	if(empFDR < previousEmpQ) empFDR = previousEmpQ;
	else previousEmpQ = empFDR;
	
	if(updateRocN)
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
	    
	    bool isdecoy = isDecoy(protein);
	    
	    if(isdecoy)
	    {
	      fpCount++;
	    }
	    else
	    {
	      tpCount++;
	    }
	    
	    if(countDecoyQvalue)
	    {
	      totalFDR += (prob);
	      estFDR = totalFDR / (tpCount + fpCount);
	    }
	    else if(tpCount)
	    {
	      if(!((bool)isdecoy)) totalFDR += (prob);
	      estFDR = totalFDR / (tpCount);
	    }
	    
	    if(tpCount) empFDR = (fpCount * pi0 * TargetDecoyRatio) / tpCount; 
	    
	    if(empFDR > 1.0 || std::isnan(empFDR) || std::isinf(empFDR)) empFDR = 1.0;
	    if(estFDR > 1.0 || std::isnan(estFDR) || std::isinf(estFDR)) estFDR = 1.0;
	    
	    if(estFDR < previousEstQ) estFDR = previousEstQ;
	    else previousEstQ = estFDR;
	    
	    if(empFDR < previousEmpQ) empFDR = previousEmpQ;
	    else previousEmpQ = empFDR;
	    
	    if(updateRocN)
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


void FidoInterface::getFDR_MSE(const std::vector<double> &estFDR, 
				 const std::vector<double> &empFDR,double &mse)
{
  /* Estimate MSE mse1 as : 1/N multiply by the SUM from k=1 to N of (estFDR(k) - empFDR(k))^2 */
  
  /* Estimate MSE mse2 area as : sum trapezoid area of each segment  (integral of the absolute value)
   * A_segment(i) = abs(X1-Xo) * abs((y1 + y2 ) / 2)
   * Where yo = estimated FDR at segment i
   * Where y1 = estimated FDR at segment i + 1
   * Where Xo = empirical FDR at segment i
   * Where X1 = empirical FDR at segment i + 1
   * Total Area = Total Area / range of X
   */
  
  /* Estimate MSE mse3 area as : sum trapezoid area with antiderivatives of each segment (absolute value of the integral)
   * A_segment(i) = ((yo - m*Xo)*X1 + m/2 * X1^2) - ((yo - m*Xo)*Xo - m/2 * X2^2))
   * Where yo = estimated FDR at segment i
   * Where y1 = estimated FDR at segment i + 1
   * Where Xo = empirical FDR at segment i
   * Where X2 = empirical FDR at segment i + 1
   * Where m = (y1 - y0) / (X1 - X0)
   * Total Area = abs(Total Area / range of X)
   */
  
   /* Estimate MSE mse4 area as : sum trapezoid squared area with antiderivatives of each segment 
   * A_segment(i) = ((yo - m*Xo)*X1 + m/2 * X1^2) - ((yo - m*Xo)*Xo - m/2 * X2^2))
   * Where yo = estimated FDR at segment i
   * Where y1 = estimated FDR at segment i + 1
   * Where Xo = empirical FDR at segment i
   * Where X2 = empirical FDR at segment i + 1
   * Where m = (y1 - y0) / (X1 - X0)
   * Total Area = Total Area / range of X
   */
  
  if(   (*min_element(estFDR.begin(),estFDR.end()) >= threshold) 
     || (estFDR.size() != empFDR.size()) 
     || (estFDR.empty() || empFDR.empty())
     || (((*max_element(estFDR.begin(),estFDR.end()) <= 0.0) 
     && (*max_element(empFDR.begin(),empFDR.end()) <= 0.0)))
    )
  {
    //no elements into the confidence interval or vectors empty 
    //or differnt size or all zeroes 
    mse = threshold;
    //mse1 = mse2 = mse3 = mse4 = 1.0;
    return;
  }
  mse = 0.0;
  //mse1 = mse2 = mse3 = mse4 = 0.0;
  double x1,x2,y1,y2;

  for(unsigned k = 0; k < estFDR.size()-1; k++) 
  {    
    if(estFDR[k] <= threshold && empFDR[k] <= threshold)
    {
      //empFDR and estFDR below threshoold, y2,x2 are the diff of them
      x1 = estFDR[k];
      x2 = estFDR[k+1];
      y1 = x1 - empFDR[k];
      y2 = x2 - empFDR[k+1];
    }
    else 
    {
      //empFDR is above threshold, penalize the area positive
      x1 = estFDR[k];
      x2 = estFDR[k+1];
      y1 = x1;
      y2 = x2;
    }
    
    if( x1 != x2 && x2 != 0 && y2 != 0 ) //if there is an area
    {
      x2 = min(x2,threshold); //in case x2 is above threshold
      //mse2 += trapezoid_area(x1,x2,y1,y2);
      //mse3 += abs(area(x1, y1, x2, y2));
      mse += areaSq(x1, y1, x2, y2);
    }
    
    //mse1 += pow(y1,2);
  }

  //mse1 += pow(y2,2); //last element of diff between vectors
  
  double normalizer1 = abs(std::min(estFDR.back(),threshold) - estFDR.front()); //normalize by x axis range (threshold on top always)
  //double normalizer2 = (double)estFDR.size(); //normalize by the number of elements
  
  //mse1 /= normalizer2;
  //mse2 /= normalizer1;
  //mse3 /= normalizer1;
  mse /= normalizer1;
  
  return;
}

string FidoInterface::printCopyright(){
  ostringstream oss;
  oss << "Copyright (c) 2008-9 University of Washington. All rights reserved.\n"
      << "Written by Oliver R. Serang (orserang@u.washington.edu) in the\n"
      << "Department of Genome Sciences at the University of Washington.\n" << std::endl;
  return oss.str();
}
