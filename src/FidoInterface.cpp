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

const double FidoInterface::kPsmThreshold = 0.0;
const double FidoInterface::kPeptideThreshold = 0.001;
const double FidoInterface::kPeptidePrior = 0.1; 
const double FidoInterface::LOG_MAX_ALLOWED_CONFIGURATIONS = 18;
const double FidoInterface::kObjectiveLambda = 0.15;

double trapezoid_area(double x1, double x2, double y1, double y2) {
  double base = abs(x1 - x2);
  double height_avg = abs((y1 + y2) / 2);
  return base * height_avg;
}

double antiderivativeAt(double m, double b, double xVal) {
  return (m*xVal*xVal/2.0 + b*xVal);
}

double squareAntiderivativeAt(double m, double b, double xVal) {
  double u = m*m;
  double v = 2*m*b;
  double t = b*b;
  return (u*xVal*xVal*xVal/3.0 + v*xVal*xVal/2.0 + t*xVal);
}

double area(double x1, double y1, double x2, double y2) {
  double m = (y2-y1)/(x2-x1);
  double b = y1-m*x1;
  double area =  antiderivativeAt(m, b, x2) - antiderivativeAt(m, b, x1);
  return area;
}

double areaSq(double x1, double y1, double x2, double y2) {
  double m = (y2-y1)/(x2-x1);
  double b = y1-m*x1;
  double area = squareAntiderivativeAt(m, b, x2) - squareAntiderivativeAt(m, b, x1);
  return area;
}

double Round(double dbVal, int nPlaces /* = 0 */) {
  const double dbShift = pow(10.0, nPlaces);
  return  floor(dbVal * dbShift + 0.5) / dbShift; 
}
    
// get the number of decimal places
int GetDecimalPlaces(double dbVal) {
  static const int MAX_DP = 10;
  static const double THRES = pow(0.1, MAX_DP);
  if (dbVal == 0.0)
    return 0;
  int nDecimal = 0;
  while (dbVal - floor(dbVal) > THRES && nDecimal < MAX_DP) {
    dbVal *= 10.0;
    nDecimal++;
  }
  return nDecimal;
}

FidoInterface::FidoInterface(double alpha, double beta, double gamma, 
    bool noPartitioning, bool noClustering, bool noPruning, 
    unsigned gridSearchDepth, double gridSearchThreshold, 
    double proteinThreshold, double mseThreshold, 
    bool usePi0, bool outputEmpirQVal, 
    std::string decoyPattern, bool trivialGrouping) :
  ProteinProbEstimator(trivialGrouping, usePi0, outputEmpirQVal, decoyPattern), 
  alpha_(alpha), beta_(beta), gamma_(gamma),
  noPartitioning_(noPartitioning), noClustering_(noClustering),
  noPruning_(noPruning), proteinThreshold_(proteinThreshold), 
  gridSearchDepth_(gridSearchDepth), 
  gridSearchThreshold_(gridSearchThreshold), mseThreshold_(mseThreshold),
  doGridSearch_(false), rocN_(kDefaultRocN) {}
      
FidoInterface::~FidoInterface() {  
  if (proteinGraph_) {
    delete proteinGraph_;
  }
  proteinGraph_ = 0;
}

void FidoInterface::run() {
  doGridSearch_ = !(alpha_ != -1 && beta_ != -1 && gamma_ != -1);
  
  double localPeptidePrior = kPeptidePrior;
  if (kComputePriors) {
    localPeptidePrior = estimatePriors();
    if (VERB > 1) {
      std::cerr << "The estimated peptide level prior probability is : " << localPeptidePrior << std::endl;
    }
  }
  
  proteinGraph_ = new GroupPowerBigraph(alpha_, beta_, gamma_, noClustering_, noPartitioning_, noPruning_, trivialGrouping_);
  proteinGraph_->setMaxAllowedConfigurations(LOG_MAX_ALLOWED_CONFIGURATIONS);
  proteinGraph_->setPeptidePrior(localPeptidePrior);
  
  if (gridSearchThreshold_ > 0.0 && doGridSearch_) {
    //NOTE lets create a smaller tree to estimate the parameters faster
    if (VERB > 1) {
      std::cerr << "Reducing the tree of proteins to increase the speed of the grid search.." << std::endl;
    }
    proteinGraph_->setProteinThreshold(gridSearchThreshold_);
    proteinGraph_->setPsmThreshold(gridSearchThreshold_);
    proteinGraph_->setPeptideThreshold(gridSearchThreshold_);
    proteinGraph_->setNoClustering(false); // i.e. do clustering
    proteinGraph_->setNoPartitioning(false); // i.e. do partitioning
    proteinGraph_->setNoPruning(false); // i.e. do pruning
    proteinGraph_->setTrivialGrouping(true);
    proteinGraph_->setMultipleLabeledPeptides(false);
  } else {
    proteinGraph_->setProteinThreshold(proteinThreshold_);
    proteinGraph_->setPsmThreshold(kPsmThreshold);
    proteinGraph_->setPeptideThreshold(kPeptideThreshold);
    proteinGraph_->setMultipleLabeledPeptides(kAddPeptideDecoyLabel);
  }
 
}

void FidoInterface::updateTargetDecoySizes() {
  std::vector<std::vector<std::string> > proteinNames;
  proteinGraph_->getProteinNames(proteinNames);
  
  numberTargetProteins_ = 0;
  numberDecoyProteins_ = 0;
  for (unsigned int k = 0; k < proteinNames.size(); ++k) {
    unsigned tpChange = countTargets(proteinNames[k]);
    unsigned fpChange = proteinNames[k].size() - tpChange;
    if (trivialGrouping_) {
      if (tpChange > 0) numberTargetProteins_ += 1;
      if (fpChange > 0) numberDecoyProteins_ += 1;
    }
  }
}

void FidoInterface::computeProbabilities(const std::string& fname) {
  ifstream fin;
  if (fname.size() > 0) {
    fin.open(fname.c_str());
    proteinGraph_->read(fin);
  } else {
    proteinGraph_->read(peptideScores_);
  }
  
  if (trivialGrouping_) updateTargetDecoySizes();
  
  time_t startTime;
  clock_t startClock;
  time(&startTime);
  startClock = clock();
  
  if (mayufdr) {
    computeFDR();
  }
  
  if (doGridSearch_) {
    if (VERB > 1) {
      std::cerr << "The parameters for the model will be estimated by grid search.\n" << std::endl;
    }
    
    if (kOptimizeParams)
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

  if (VERB > 1) {
    cerr << "The following parameters have been chosen:\n";
    std::cerr.precision(10);
    cerr << "alpha = " << alpha_ << endl;
    cerr << "beta  = " << beta_ << endl;
    cerr << "gamma = " << gamma_ << endl;
    std::cerr.unsetf(std::ios::floatfield);
    cerr << "\nProtein level probabilities will now be estimated\n";
  }


  if (gridSearchThreshold_ > 0.0 && doGridSearch_) {
    //NOTE reset the tree after grid searching
    proteinGraph_->setProteinThreshold(proteinThreshold_);
    proteinGraph_->setPsmThreshold(kPsmThreshold);
    proteinGraph_->setPeptideThreshold(kPeptideThreshold);
    proteinGraph_->setNoClustering(noClustering_);
    proteinGraph_->setNoPartitioning(noPartitioning_);
    proteinGraph_->setNoPruning(noPruning_);
    proteinGraph_->setTrivialGrouping(trivialGrouping_);
    proteinGraph_->setMultipleLabeledPeptides(kAddPeptideDecoyLabel);
    
    if (fname.size() > 0) {
      proteinGraph_->read(fin);
    } else {
      proteinGraph_->read(peptideScores_);
    }
    if (trivialGrouping_) updateTargetDecoySizes();
  }
  
  proteinGraph_->setAlphaBetaGamma(alpha_, beta_, gamma_);
  proteinGraph_->getProteinProbs();
  pepProteinMap_.clear();
  proteinGraph_->getProteinProbsPercolator(pepProteinMap_);
}

void FidoInterface::gridSearch() {
  std::vector<double> gamma_search, beta_search, alpha_search;
  
  switch(gridSearchDepth_) {
    case 0:    
      gamma_search.push_back(0.5);
      
      beta_search.push_back(0.001);
      
      alpha_search.push_back(0.008);
      alpha_search.push_back(0.032);
      alpha_search.push_back(0.128);
      break;
    
    case 1:
      gamma_search.push_back(0.1);
      gamma_search.push_back(0.5);
      gamma_search.push_back(0.9);
      
      beta_search.push_back(0.001);
      for (double k = 0.002; k <= 0.4; k*=4) {
       alpha_search.push_back(k);
      }
      break;
      
    case 2:
      gamma_search.push_back(0.1);
      gamma_search.push_back(0.3);
      gamma_search.push_back(0.5);
      gamma_search.push_back(0.75);
      gamma_search.push_back(0.9);
      
      beta_search.push_back(0.001);
      for (double k = 0.001; k <= 0.4; k*=2) {
       alpha_search.push_back(k);
      }
      break;
    
    default:
      gamma_search.push_back(0.5);
      
      beta_search.push_back(0.001);
      
      alpha_search.push_back(0.008);
      alpha_search.push_back(0.032);
      alpha_search.push_back(0.128);
      break;
  }

  if (alpha_ != -1) alpha_search.push_back(alpha_);
  if ( beta_ != -1) beta_search.push_back(beta_);
  if (gamma_ != -1) gamma_search.push_back(gamma_);
  
  gridSearch(alpha_search, beta_search, gamma_search);
}

void FidoInterface::gridSearchOptimize() {
  if (VERB > 1) {
    std::cerr << "Running super grid search..." << std::endl;
  }
  
  double alpha_step = 0.05;
  double beta_step = 0.05;
  double gamma_step = 0.05;
  
  double beta_init = 0.00001;
  double alpha_init = 0.001;
  double gamma_init = 0.1;
  
  double gamma_limit = 0.5;
  double beta_limit = 0.05;
  double alpha_limit = 0.5;
  
  if (alpha_ != -1) alpha_init = alpha_limit = alpha_;
  if ( beta_ != -1) beta_init = beta_limit = beta_;
  if (gamma_ != -1) gamma_init = gamma_limit = gamma_;
  
  std::vector<double> gamma_search, beta_search, alpha_search;
  
  //NOTE very annoying the residue error of the floats that get acummulated in every iteration
  for (int i = 0; gamma_init + gamma_step*i <= gamma_limit; ++i) { 
    gamma_search.push_back(gamma_init + gamma_step*i);
  }
  
  double lbi = log10(beta_init);
  for (int i = 0; lbi + i*beta_step <= Round(log10(beta_limit),2); ++i) {
    double j = lbi + i*beta_step;
    double original = pow(10,j);
    double beta_local = original - beta_init;
    if (beta_local > 0.0) beta_local = original;
    beta_search.push_back(beta_local);
  }
  
  double lai = log10(alpha_init);
  for (int i = 0; lai + i*alpha_step <= Round(log10(alpha_limit),2); ++i) {
    double k = lai + i*alpha_step;
    alpha_search.push_back(pow(10,k));
  }
  
  gridSearch(alpha_search, beta_search, gamma_search);
}

void FidoInterface::gridSearch(std::vector<double>& alpha_search, 
    std::vector<double>& beta_search, 
    std::vector<double>& gamma_search) {
  double gamma_best = -1.0, alpha_best = -1.0, beta_best = -1.0;
  double best_objective = -100000000;
  double current_objective;
  
  for (unsigned int i = 0; i < gamma_search.size(); i++) {
    double gamma_local = gamma_search[i];
    for (unsigned int j = 0; j < alpha_search.size(); j++) {
      double alpha_local = alpha_search[j];
      for (unsigned int k = 0; k < beta_search.size(); k++) {
        double beta_local = beta_search[k];
        
        current_objective = calcObjective(alpha_local, beta_local, gamma_local);
        if (current_objective > best_objective) {
          best_objective = current_objective;
          gamma_best = gamma_local;
          alpha_best = alpha_local;
          beta_best = beta_local;
        }
      }
    }
  }
  alpha_ = alpha_best;
  beta_ = beta_best;
  gamma_ = gamma_best;
}

double FidoInterface::calcObjective(double alpha, double beta, double gamma) {
  std::vector<std::vector<std::string> > names;
  std::vector<double> probs, empq, estq; 
  double roc ,mse, objective;
  
  proteinGraph_->setAlphaBetaGamma(alpha, beta, gamma);
  proteinGraph_->getProteinProbs();
  proteinGraph_->getProteinProbsAndNames(names, probs);
  
  getEstimated_and_Empirical_FDR(names, probs, empq, estq);
  getROC_AUC(names, probs, roc);
  getFDR_MSE(estq, empq, mse);
  
  objective = (kObjectiveLambda * roc) - fabs((1-kObjectiveLambda) * mse);
  
  if (VERB > 2) {
    std::cerr.precision(10);
    std::cerr << "Grid searching Alpha= "  << alpha << 
                 " Beta= " << beta << 
                 " Gamma= "  << gamma << std::endl;
    std::cerr.unsetf(std::ios::floatfield);
    std::cerr << "The ROC AUC estimated values is : " << roc << std::endl;
    std::cerr << "The MSE FDR estimated values is : " << mse << std::endl;
    std::cerr << "Objective function with second roc and mse is : " << 
                 objective << std::endl;
  }
  return objective;
}

void FidoInterface::getROC_AUC(const std::vector<std::vector<string> > &names,
            const std::vector<double> &probabilities, double &auc) {
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
  
  // assuming names and probabilities same size; rocN_ set by getEstimated_and_Empirical_FDR()
  for (unsigned k = 0; k < names.size() && fp <= rocN_; k++) {
    double prob = probabilities[k];
    unsigned tpChange = countTargets(names[k]);
    unsigned fpChange = names[k].size() - tpChange;
    //if ties activated count groups as 1 protein
    if (trivialGrouping_) {
      if (tpChange > 0) tpChange = 1;
      if (fpChange > 0) fpChange = 1;
    }
    tp += tpChange;
    fp += fpChange;
    //should only do it when fp changes and either of them is != 0
    if (prev_prob != -1 && fp != 0 && tp != 0 && fp != prev_fp) {
      double trapezoid = trapezoid_area(fp,prev_fp,tp,prev_tp);
      prev_fp = fp;
      prev_tp = tp;
      auc += trapezoid;
    }    
    prev_prob = prob;
  }

  unsigned normalizer = (tp * fp);
  
  if (normalizer > 0) {
    auc /= normalizer;
  } else {
    auc = 0.0;
  }
  
  return;
}


void FidoInterface::getEstimated_and_Empirical_FDR(
    const std::vector<std::vector<string> >& proteinNames,
    const std::vector<double>& probabilities,
    std::vector<double>& empq, 
    std::vector<double>& estq) {
  empq.clear();
  estq.clear();
  
  double targetDecoyRatio = static_cast<double>(numberTargetProteins_) / numberDecoyProteins_;
  FDRCalculator fdrCalculator(usePi0_, targetDecoyRatio, pi0_, countDecoyQvalue_);
  
  //NOTE no need to store more q values since they will not be taken into account while estimating MSE FDR divergence
  for (unsigned int k = 0; (k < proteinNames.size() && 
        (fdrCalculator.getPreviousEstQ() <= mseThreshold_)); k++) {
    double prob = probabilities[k];
    unsigned tpChange = countTargets(proteinNames[k]);
    unsigned fpChange = proteinNames[k].size() - tpChange;
    if (trivialGrouping_) {
      if (tpChange > 0) tpChange = 1;
      if (fpChange > 0) fpChange = 1;
    }
    fdrCalculator.calcFDRs(fpChange, tpChange, prob, estq, empq);
  }
  if (kUpdateRocN) rocN_ = fdrCalculator.getRocN();
}


void FidoInterface::getFDR_MSE(const std::vector<double> &estFDR, 
         const std::vector<double> &empFDR,double &mse) {
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
  
  if(   (*min_element(estFDR.begin(),estFDR.end()) >= mseThreshold_) 
     || (estFDR.size() != empFDR.size()) 
     || (estFDR.empty() || empFDR.empty())
     || (((*max_element(estFDR.begin(),estFDR.end()) <= 0.0) 
     && (*max_element(empFDR.begin(),empFDR.end()) <= 0.0))) ) {
    //no elements into the confidence interval or vectors empty 
    //or different size or all zeroes 
    mse = mseThreshold_;
    //mse1 = mse2 = mse3 = mse4 = 1.0;
    return;
  }
  mse = 0.0;
  //mse1 = mse2 = mse3 = mse4 = 0.0;
  double x1,x2,y1,y2;

  for (unsigned k = 0; k < estFDR.size()-1; k++) {
    if (estFDR[k] <= mseThreshold_ && empFDR[k] <= mseThreshold_) {
      //empFDR and estFDR below threshold, y2,x2 are the diff of them
      x1 = estFDR[k];
      x2 = estFDR[k+1];
      y1 = x1 - empFDR[k];
      y2 = x2 - empFDR[k+1];
    } else {
      //empFDR is above mseThreshold_, penalize the area positive
      x1 = estFDR[k];
      x2 = estFDR[k+1];
      y1 = x1;
      y2 = x2;
    }
    
    
    if ( x1 != x2 && x2 != 0 && y2 != 0 ) { //if there is an area
      x2 = min(x2,mseThreshold_); //in case x2 is above mseThreshold_
      //mse2 += trapezoid_area(x1,x2,y1,y2);
      //mse3 += abs(area(x1, y1, x2, y2));
      mse += areaSq(x1, y1, x2, y2);
    }
    
    //mse1 += pow(y1,2);
  }

  //mse1 += pow(y2,2); //last element of diff between vectors
  
  double normalizer1 = abs(std::min(estFDR.back(),mseThreshold_) - estFDR.front()); //normalize by x axis range (mseThreshold_ on top always)
  //double normalizer2 = (double)estFDR.size(); //normalize by the number of elements
  
  //mse1 /= normalizer2;
  //mse2 /= normalizer1;
  //mse3 /= normalizer1;
  mse /= normalizer1;
  
  return;
}

std::ostream& FidoInterface::printParametersXML(std::ostream &os) {
  os << "    <alpha>" << alpha_ <<"</alpha>" << endl;
  os << "    <beta>"  << beta_ <<"</beta>" << endl;
  os << "    <gamma>" << gamma_ <<"</gamma>" << endl;
  return os;
}

string FidoInterface::printCopyright() {
  ostringstream oss;
  oss << "Copyright (c) 2008-9 University of Washington. All rights reserved.\n"
      << "Written by Oliver R. Serang (orserang@u.washington.edu) in the\n"
      << "Department of Genome Sciences at the University of Washington.\n" << std::endl;
  return oss.str();
}

void FDRCalculator::calcFDRs(double fpChange, double tpChange, double prob,
               std::vector<double>& empq, std::vector<double>& estq) {
  double estFDR = 0.0, empFDR = 0.0;
  fpCount_ += fpChange;
  tpCount_ += tpChange;
  
  if (countDecoyQvalue_) {
    totalFDR_ += prob * (tpChange + fpChange);
    estFDR = totalFDR_ / (tpCount_ + fpCount_);
  } else if (tpCount_ > 0) {
    totalFDR_ += prob * tpChange;
    estFDR = totalFDR_ / tpCount_;
  }

  if (tpCount_ > 0) {
    if (usePi0_) {
      empFDR = (fpCount_ * pi0_ * targetDecoyRatio_) / tpCount_;
    } else {
      empFDR = fpCount_ / tpCount_;
    }
  }
  
  if (empFDR > 1.0 || std::isnan(empFDR) || std::isinf(empFDR)) empFDR = 1.0;
  if (estFDR > 1.0 || std::isnan(estFDR) || std::isinf(estFDR)) estFDR = 1.0;
      
  if (estFDR < previousEstQ_) estFDR = previousEstQ_;
  else previousEstQ_ = estFDR;
      
  if (empFDR < previousEmpQ_) empFDR = previousEmpQ_;
  else previousEmpQ_ = empFDR;
  
  estq.push_back(estFDR);
  empq.push_back(empFDR);
  
  rocN_ = std::max(rocN_, (unsigned)std::max(50,std::min((int)fpCount_,500)));
}
