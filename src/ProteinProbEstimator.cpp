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

#include "ProteinProbEstimator.h"

const double ProteinProbEstimator::target_decoy_ratio = 1.0;
const double ProteinProbEstimator::psmThresholdMayu = 0.90;
const double ProteinProbEstimator::prior_protein = 0.5;
bool ProteinProbEstimator::calcProteinLevelProb = false;
/** Helper functions **/

template<class T> 
void bootstrap(const vector<T>& in, vector<T>& out, size_t max_size = 1000) {
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


ProteinProbEstimator::ProteinProbEstimator(bool trivialGrouping, double pi0, 
					     bool outputEmpirQVal, std::string decoyPattern) : 
	  trivialGrouping_(trivialGrouping), pi0_(pi0), peptideScores_(NULL), 
	  numberDecoyProteins_(0u), numberTargetProteins_(0u), usePi0_(true),
	  outputEmpirQVal_(outputEmpirQVal), decoyPattern_(decoyPattern), fdr_(1.0) {}

ProteinProbEstimator::~ProteinProbEstimator() {
  FreeAll(qvalues);
  FreeAll(qvaluesEmp);
  FreeAll(pvalues);
  
  if (mayufdr && fastReader) {
    delete fastReader;
  }
  fastReader = 0;
  
  for(std::multimap<double,std::vector<std::string> >::iterator it = pepProteinMap_.begin();
        it != pepProteinMap_.end(); it++) {
	  FreeAll(it->second);
  }
      
  for(std::map<const std::string,Protein*>::iterator it = proteins.begin(); 
        it != proteins.end(); it++) {
	  if(it->second)
	    delete it->second;
  }
}

bool ProteinProbEstimator::initialize(Scores* fullset) {
  peptideScores_ = fullset;
  setTargetandDecoysNames();
  return true;
}

void ProteinProbEstimator::computeFDR() {
  if (VERB > 1)
    std::cerr << "Estimating Protein FDR ... " << std::endl;
    
  fastReader = new ProteinFDRestimator();
  fastReader->setDecoyPrefix(decoyPattern_);
  fastReader->setTargetDecoyRatio(target_decoy_ratio);
  fastReader->setEqualDeepthBinning(binning_equal_deepth);
  fastReader->setNumberBins(number_bins);
  fastReader->correctIdenticalSequences(targetProteins,decoyProteins);
  //These guys are the number of target and decoys proteins but from the subset of PSM with FDR < threshold
  std::set<std::string> numberTP;
  std::set<std::string> numberFP;
  getTPandPFfromPeptides(psmThresholdMayu,numberTP,numberFP);
    
  double fptol = fastReader->estimateFDR(numberTP,numberFP);
    
  if (fptol == -1) {
    fdr_ = 1.0;
    
    if(VERB > 1)
      std::cerr << "There was an error estimating the Protein FDR..\n" << std::endl;
  } else {	
    fdr_ = (fptol/(double)numberTP.size());
    if(fdr_ <= 0 || fdr_ >= 1.0) fdr_ = 1.0;      
    
    if(VERB > 1) {
      std::cerr << "Estimated Protein FDR at ( " << psmThresholdMayu << ") PSM FDR is : " 
      << fdr_ << " with " << fptol << " expected number of false positives proteins\n" << std::endl;
    }
  }
}

void ProteinProbEstimator::computeStatistics() {
  if (usePi0_ && !mayufdr && outputEmpirQVal_ && pvalues.size() == 0) { 
    estimatePValues();
    /*pi0_ = estimatePi0();
    if (pi0_ <= 0.0 || pi0_ > 1.0) pi0_ = *qvalues.rbegin();
  } else {
    pi0_ = fdr_;
    */
  }
  
  /** computing q values **/
  estimateQValues();
  estimateQValuesEmp();
  updateProteinProbabilities();

  if (VERB > 1) {
    if (trivialGrouping_) {
      std::cerr << "Number of protein groups identified at q-value = 0.01: ";
    } else {
      std::cerr << "Number of proteins identified at q-value = 0.01: ";
    }
    std::cerr << getQvaluesBelowLevel(0.01) << std::endl;
  }
}

void ProteinProbEstimator::printOut(const std::string &proteinFN, 
				      const std::string &proteinDecoyFN) {
  if(!proteinFN.empty() || !proteinDecoyFN.empty()) {
    if(!proteinFN.empty()) {
      ofstream proteinOut(proteinFN.data(), ios::out);
      print(proteinOut,false);
      proteinOut.close();	
    }
    if(!proteinDecoyFN.empty()) {
      ofstream proteinOut(proteinDecoyFN.data(), ios::out);
      print(proteinOut,true);
      proteinOut.close();
    }
  } else {
    print(std::cout);
  }
}

double ProteinProbEstimator::estimatePriors() {
  /* Compute a priori probabilities of peptide presence (correctness?) */
  /* prior = the mean of the probabilities, maybe one prior for each charge *
   * prior2 = assuming a peptide is present if only if the protein is present and counting
   * the size of protein and prior protein probabily in the computation
   * prior3 = the ratio of confident peptides among all the peptides 
   * prior4 = the ratio of decoy peptides versus target peptides (TDC) */
  
  double prior_peptide = 0.0;
  double prior_peptide2 = 0.0;
  unsigned confident_peptides = 0;
  unsigned decoy_peptides = 0;
  unsigned total_peptides = 0;
  double prior, prior2, prior3, prior4;
  prior = prior2 = prior3 = prior4 = 0.0;
  for (vector<ScoreHolder>::iterator psm = peptideScores_->begin(); 
         psm!= peptideScores_->end(); ++psm) {
    if(!psm->isDecoy()) {
      unsigned size = psm->pPSM->proteinIds.size();
      double prior = prior_protein * size;
      double tmp_prior = prior;
      // for each protein
      for(set<string>::iterator protIt = psm->pPSM->proteinIds.begin(); 
	          protIt != psm->pPSM->proteinIds.end(); protIt++) {
	      unsigned index = std::distance(psm->pPSM->proteinIds.begin(), protIt);
	      tmp_prior = (tmp_prior * prior_protein * (size - index)) / (index + 1);
	      prior +=  pow(-1.0,(int)index) * tmp_prior;
      }
      /* update computed prior */
      prior_peptide += (1.0-prior);
      if(psm->q <= 0.1) ++confident_peptides;
      prior_peptide2 += (1.0-psm->pep);
      ++total_peptides;
    } else {
      ++decoy_peptides;
    }
  }
  
  prior = prior_peptide2 / (double)total_peptides;
  prior2 = prior_peptide / (double)total_peptides;
  prior3 = confident_peptides / (double)total_peptides;
  prior4 = (total_peptides - decoy_peptides) / (double)total_peptides;
  
  if (prior > 0.99) prior = 0.99;
  if (prior < 0.01) prior = 0.01;
  if (prior2 > 0.99) prior2 = 0.99;
  if (prior2 < 0.01) prior2 = 0.01;
  if (prior3 > 0.99) prior3 = 0.99;
  if (prior3 < 0.01) prior3 = 0.01;
  if (prior4 > 0.99) prior4 = 0.99;
  if (prior4 < 0.01) prior4 = 0.01;
  
  return prior4;
}

void ProteinProbEstimator::getCombinedList(std::vector<std::pair<double , bool> > &combined) {
  for (std::multimap<double,std::vector<std::string> >::const_iterator it = pepProteinMap_.begin();
       it != pepProteinMap_.end(); it++) {
    double prob = it->first;
    std::vector<std::string> proteinList = it->second;
    for(std::vector<std::string>::const_iterator itP = proteinList.begin();
	        itP != proteinList.end(); itP++) {
      std::string proteinName = *itP;
      bool isdecoy = proteins[proteinName]->getIsDecoy();
      combined.push_back(std::make_pair(prob,isdecoy));
    }
  }
}

void ProteinProbEstimator::estimatePValues() {
  // assuming combined sorted in best hit first order
  std::vector<std::pair<double , bool> > combined;
  getCombinedList(combined);
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

void ProteinProbEstimator::getTPandPFfromPeptides(double psm_threshold, 
						  std::set<std::string> &numberTP, 
						  std::set<std::string> &numberFP) {
  /* The original paper of Mayu describes a protein as :
   * FP = if only if all its peptides with q <= threshold are decoy
   * TP = at least one of its peptides with q <= threshold is target
   * However, the way mayu estimates it on the program is like this :
   * FP = any protein that contains a decoy psm with q <= threshold
   * TP = any protein that contains a target psm with q <= threshold
   * They do not consider protein containing both decoy and target psms
   * Also, mayu estimates q as the empirical (target-decoy) q value.
   * Percolator estimates q as the empirical (target-decoy) q value and adjusted by pi0
   * Mayu extracts the list of TP and FP proteins from PSM level whereas percolator
   * extract the list of TP and FP proteins from peptide level, this avoids redundancy and
   * gives a better calibration since peptide level q values are re-adjusted in percolator.
   * This creates sometimes a difference in the number of TP and FP proteins between percolator and Mayus 
   * which causes a slight difference in the estimated protein FDR
   */
  for (std::map<std::string,Protein*>::const_iterator it = proteins.begin();
       it != proteins.end(); it++) {
    unsigned num_target_confident = 0;
    unsigned num_decoy_confident = 0;
    std::string protname = it->first;
    std::vector<Protein::Peptide*> peptides = it->second->getPeptides();
    for(std::vector<Protein::Peptide*>::const_iterator itP = peptides.begin();
          itP != peptides.end(); itP++) {
      Protein::Peptide *p = *itP;
      if(p->q <= psm_threshold && p->isdecoy)
        ++num_decoy_confident;
      if(p->q <= psm_threshold && !p->isdecoy)
        ++num_target_confident;
    }
    if (num_decoy_confident > 0) {
      numberFP.insert(protname);
    }
    if (num_target_confident > 0) {
      numberTP.insert(protname);
    }
  }
}

double ProteinProbEstimator::estimatePi0(const unsigned int numBoot) 
{
  std::vector<double> pBoot, lambdas, pi0s, mse;
  std::vector<double>::iterator start;
  unsigned int numLambda = 100;
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
  if (pi0s.size()==0) {
    cerr << "Error in the input data: too good separation between target "
        << "and decoy Proteins.\nImpossible to estimate pi0. Taking the highest estimated q value as pi0.\n";
    return -1;
  }
  double minPi0 = *min_element(pi0s.begin(), pi0s.end());
  // Initialize the vector mse with zeroes.
  fill_n(back_inserter(mse), pi0s.size(), 0.0);
  // Examine which lambda level that is most stable under bootstrap
  for (unsigned int boot = 0; boot < numBoot; ++boot) {
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

unsigned ProteinProbEstimator::getQvaluesBelowLevel(double level) {
  std::set<int> identifiedGroupIds;
  unsigned nP = 0;
  for (std::map<const std::string,Protein*>::const_iterator myP = proteins.begin(); 
          myP != proteins.end(); ++myP) {
    if (myP->second->getQ() < level && !myP->second->getIsDecoy()) {
      nP++;
      identifiedGroupIds.insert(myP->second->getGroupId());
    }
  }
  
  if (trivialGrouping_) {
    return identifiedGroupIds.size();
  } else {
    return nP;
  }
}

unsigned ProteinProbEstimator::getQvaluesBelowLevelDecoy(double level) { 
  std::set<int> identifiedGroupIds;
  unsigned nP = 0;
  for (std::map<const std::string,Protein*>::const_iterator myP = proteins.begin(); 
          myP != proteins.end(); ++myP) {
    if (myP->second->getQ() < level && myP->second->getIsDecoy()) {
      nP++;
      identifiedGroupIds.insert(myP->second->getGroupId());
    }
  }
  
  if (trivialGrouping_) {
    return identifiedGroupIds.size();
  } else {
    return nP;
  }
}


void ProteinProbEstimator::estimateQValues() {
  unsigned nP = 0;
  double sum = 0.0;
  double qvalue = 0.0;
  qvalues.clear();
  
  for (std::multimap<double,std::vector<std::string> >::const_iterator it = pepProteinMap_.begin(); 
          it != pepProteinMap_.end(); it++) {
    int ntargets = countTargets(it->second);
    int ndecoys = it->second.size() - ntargets;
    if (trivialGrouping_) {
      if (ntargets > 0) ntargets = 1;
      if (ndecoys > 0) ndecoys = 1;
    }
    
    if (ndecoys == 0) {
      //NOTE in case I want to count and use target and decoys proteins while estimating qvalue from PEP
      if (countDecoyQvalue_) {
        sum += (double)(it->first * (ntargets + ndecoys));
        nP += (ntargets + ndecoys);
      } else {
        sum += (double)(it->first * ntargets);
        nP += ntargets;
      }
    }
    qvalue = (sum / (double)nP);
    if (std::isnan(qvalue) || std::isinf(qvalue) || qvalue > 1.0) qvalue = 1.0;
    
    for (size_t i = 0; i < it->second.size(); ++i) {
      qvalues.push_back(qvalue);
    }
  }
  std::partial_sum(qvalues.rbegin(),qvalues.rend(),qvalues.rbegin(),myminfunc);
}

void ProteinProbEstimator::estimateQValuesEmp() {
    // assuming combined sorted in decending order
  unsigned nDecoys = 0;
  unsigned numTarget = 0;
  unsigned nTargets = 0;
  double qvalue = 0.0;
  unsigned numDecoy = 0;
  //pvalues.clear();
  qvaluesEmp.clear();
  double targetDecoyRatio = (double)numberTargetProteins_ / (double)numberDecoyProteins_;
  
  if (VERB > 1) {
    if (trivialGrouping_) {
      std::cerr << "Number of target protein groups: " << numberTargetProteins_ << std::endl;
      std::cerr << "Number of decoy protein groups: " << numberDecoyProteins_ << std::endl;
    } else {
      std::cerr << "Number of target proteins: " << numberTargetProteins_ << std::endl;
      std::cerr << "Number of decoy proteins: " << numberDecoyProteins_ << std::endl;
    }
  }
  
  for (std::multimap<double,std::vector<std::string> >::const_iterator it = pepProteinMap_.begin(); it != pepProteinMap_.end(); it++) {
    numTarget = countTargets(it->second);
    numDecoy = it->second.size() - numTarget;
    if (trivialGrouping_) {
      if (numTarget > 0) numTarget = 1;
      if (numDecoy > 0) numDecoy = 1;
    }
    nDecoys += numDecoy;
    nTargets += numTarget;
    
    if (nTargets) {
      if (usePi0_) {
        qvalue = (double)(nDecoys * pi0_ * targetDecoyRatio) / (double)nTargets;
      } else {
        qvalue = (double)(nDecoys) / (double)nTargets;
      }
    }
    if (std::isnan(qvalue) || std::isinf(qvalue) || qvalue > 1.0) qvalue = 1.0;
    
    for (size_t i = 0; i < it->second.size(); ++i) {
      qvaluesEmp.push_back(qvalue);
      /*
      if (numDecoy > 0) {
        pvalues.push_back((nDecoys)/(double)(numberDecoyProteins_));
      } else {
        pvalues.push_back((nDecoys+(double)1)/(numberDecoyProteins_+(double)1));
      }
      */
    }
  }
  std::partial_sum(qvaluesEmp.rbegin(), qvaluesEmp.rend(), qvaluesEmp.rbegin(), myminfunc);
}

void ProteinProbEstimator::updateProteinProbabilities() {
  std::vector<double> peps; // posterior error probabilities, not peptides
  std::vector<std::vector<std::string> > proteinNames;
  std::transform(pepProteinMap_.begin(), pepProteinMap_.end(), std::back_inserter(peps), RetrieveKey());
  std::transform(pepProteinMap_.begin(), pepProteinMap_.end(), std::back_inserter(proteinNames), RetrieveValue());
  unsigned protIdx = 0;
  for (unsigned pepIdx = 0; pepIdx < peps.size(); pepIdx++) {
    double pep = peps[pepIdx];
    std::vector<std::string> proteinlist = proteinNames[pepIdx];
    for (unsigned j = 0; j < proteinlist.size(); j++) { 
      int protGroupId = protIdx + 1;
      if (trivialGrouping_) protGroupId = pepIdx + 1;
      
      std::string proteinName = proteinlist[j];
      proteins[proteinName]->setPEP(pep);
      proteins[proteinName]->setQ(qvalues[protIdx]);
      proteins[proteinName]->setQemp(qvaluesEmp[protIdx]);
      proteins[proteinName]->setP(pvalues[protIdx]);
      proteins[proteinName]->setGroupId(protGroupId);
      ++protIdx;
    }
  }

}

void ProteinProbEstimator::setTargetandDecoysNames() {
  unsigned int numGroups = 0;
  for (vector<ScoreHolder>::iterator psm = peptideScores_->begin(); psm!= peptideScores_->end(); ++psm) {
    // for each protein
    for (set<string>::iterator protIt = psm->pPSM->proteinIds.begin(); protIt != psm->pPSM->proteinIds.end(); protIt++) {
      Protein::Peptide *peptide = new Protein::Peptide(psm->pPSM->getPeptideSequence(),psm->isDecoy(),psm->p,
							psm->pep,psm->q,psm->p);
      if (proteins.find(*protIt) == proteins.end()) {
	      Protein *newprotein = new Protein(*protIt,0.0,0.0,0.0,0.0,psm->isDecoy(),peptide,++numGroups);
	      proteins.insert(std::make_pair(*protIt,newprotein));
	
	      if (psm->isDecoy()) {
	        falsePosSet.insert(*protIt);
	      } else {
	        truePosSet.insert(*protIt);
	      }
      } else {
      	proteins[*protIt]->setPeptide(peptide);
      }
    }
  }  
  numberDecoyProteins_ = falsePosSet.size();
  numberTargetProteins_ = truePosSet.size();
}

void ProteinProbEstimator::addProteinDb(bool isDecoy, std::string name, std::string sequence, double length) {
  if (isDecoy)
    decoyProteins.insert(std::make_pair(name,std::make_pair(sequence,length)));
  else
    targetProteins.insert(std::make_pair(name,std::make_pair(sequence,length)));
}

unsigned ProteinProbEstimator::countTargets(const std::vector<std::string> &proteinList) {
  unsigned count = 0;
  for (std::vector<std::string>::const_iterator it = proteinList.begin(); it != proteinList.end(); it++) {
    if (useDecoyPrefix) {
      if ((*it).find(decoyPattern_) == std::string::npos) {
      	count++;
      }
    } else {
      if (truePosSet.count(*it) != 0) {
      	count++;
      }
    }
  }
  return count;
}

unsigned ProteinProbEstimator::countDecoys(const std::vector<std::string> &proteinList) {
  unsigned count = 0;
  for(std::vector<std::string>::const_iterator it = proteinList.begin(); it != proteinList.end(); it++) {
    if (useDecoyPrefix) {
      if((*it).find(decoyPattern_) != std::string::npos) {
      	count++;
      }
    } else {
      if(falsePosSet.count(*it) != 0) {
      	count++;
      }
    }
  }
  return count;
}

bool ProteinProbEstimator::isDecoy(const std::string& proteinName) {
  //NOTE faster with decoyPrefix but I assume the label that identifies decoys is in decoyPattern_
   return (bool)(useDecoyPrefix ? proteinName.find(decoyPattern_) 
	    != std::string::npos : falsePosSet.count(proteinName) != 0);
}
    
bool ProteinProbEstimator::isTarget(const std::string& proteinName) {
  //NOTE faster with decoyPrefix but I assume the label that identifies decoys is in decoyPattern_
   return (bool)(useDecoyPrefix ? proteinName.find(decoyPattern_) 
		  == std::string::npos : truePosSet.count(proteinName) != 0);
}


void ProteinProbEstimator::writeOutputToXML(string xmlOutputFN, bool outputDecoys) {
  std::vector<std::pair<std::string,Protein*> > myvec(proteins.begin(), proteins.end());
  std::sort(myvec.begin(), myvec.end(), IntCmpProb());

  ofstream os;
  os.open(xmlOutputFN.data(), ios::app);
  // append PROTEINs tag
  os << "  <proteins>" << endl;
  for (std::vector<std::pair<std::string,Protein*> > ::const_iterator myP = myvec.begin(); 
	      myP != myvec.end(); myP++) {
    if ( (!outputDecoys && !myP->second->getIsDecoy()) || (outputDecoys)) {
      os << "    <protein p:protein_id=\"" << myP->second->getName() << "\"";
	    if (outputDecoys) {
	      if (myP->second->getIsDecoy()) 
	        os << " p:decoy=\"true\"";
	      else  
	        os << " p:decoy=\"false\"";
	    }
	    os << ">" << endl;
	    
	    os << "      <pep>" << scientific << myP->second->getPEP() << "</pep>" << endl;
	  
	    if (outputEmpirQVal_) {
	      os << "      <q_value_emp>" << scientific << myP->second->getQemp() << "</q_value_emp>\n";
	    }
	  
	    os << "      <q_value>" << scientific << myP->second->getQ() << "</q_value>\n";
	    
	    if (outputEmpirQVal_) {
	      os << "      <p_value>" << scientific << myP->second->getP() << "</p_value>\n";
	    }
	  
	    std::vector<Protein::Peptide*> peptides = myP->second->getPeptides();
	    for (std::vector<Protein::Peptide*>::const_iterator peptIt = peptides.begin(); 
	        peptIt != peptides.end(); peptIt++) {
	      if ((*peptIt)->name != "") {
	        os << "      <peptide_seq seq=\"" << (*peptIt)->name << "\"/>"<<endl;
	      }
	    }
	    os << "    </protein>" << endl;
	  }
  }
    
  os << "  </proteins>" << endl << endl;
  os.close();
}

void ProteinProbEstimator::print(ostream& myout, bool decoy) {
  
  std::vector<std::pair<std::string,Protein*> > myvec(proteins.begin(), proteins.end());
  std::sort(myvec.begin(), myvec.end(), IntCmpProb());
  
  myout << "ProteinId\tProteinGroupId\tq-value\tposterior_error_prob\tpeptideIds" << std::endl;
  
  for (std::vector<std::pair<std::string,Protein*> > ::const_iterator myP = myvec.begin(); 
	        myP != myvec.end(); myP++) {
    if( (decoy && myP->second->getIsDecoy()) || (!decoy && !myP->second->getIsDecoy())) {
      myout << myP->second->getName() << "\t" << myP->second->getGroupId() << "\t" 
            << myP->second->getQ() << "\t" << myP->second->getPEP() << "\t";
      std::vector<Protein::Peptide*> peptides = myP->second->getPeptides();
      for(std::vector<Protein::Peptide*>::const_iterator peptIt = peptides.begin(); peptIt != peptides.end(); peptIt++) {
        if((*peptIt)->name != "") {
          myout << (*peptIt)->name << "  ";
        }
      }
      myout << std::endl;
    }
  }
}
