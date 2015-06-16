// Written by Oliver Serang 2009
// see license for more information

#include "GroupPowerBigraph.h"

GroupPowerBigraph::~GroupPowerBigraph() { }

Array<double> GroupPowerBigraph::proteinProbs() {
  Array<double> result;
  for (int k = 0; k < subgraphs_.size(); k++) {
    subgraphs_[k].getProteinProbs(params_);
    result.append( subgraphs_[k].proteinProbabilities() );
  }
  return result;
}

void GroupPowerBigraph::getProteinProbs() {
  probsPresentProteins_ = proteinProbs();
}

void GroupPowerBigraph::getGroupProtNames() {
  for (int k = 0; k < subgraphs_.size(); k++) {
    const BasicGroupBigraph& bgb = subgraphs_[k];
    for (int j = 0; j < bgb.proteinsToPSMs.size(); j++) {
      groupProtNames_.add(bgb.proteinGroupNames()[j]);
    }
  }
}

void GroupPowerBigraph::printProteinWeights() const {
  //NOTE this is printing out PEPs not posterior probabilities
  Array<double> sorted = probsPresentProteins_;
  Array<int> indices = sorted.sort();
  for (int k = 0; k < sorted.size(); k++) {
    cout << 1.0 - sorted[k] << " " << groupProtNames_[ indices[k] ] << endl;
  }
  if (severedProteins_.size() != 0)
    cout << "1.0 " << severedProteins_ << endl;
}

// return a map of PEPs and their respectives proteins
void GroupPowerBigraph::getProteinProbsPercolator(std::multimap<double, std::vector<std::string> > &pepProteins) const {
  Array<double> sorted = probsPresentProteins_;
  Array<int> indices = sorted.sort();
  for (int k = 0; k < sorted.size(); k++) {
    double pep = (1.0 - sorted[k]);
    if (pep < 0.0) pep = 0.0;
    if (pep > 1.0) pep = 1.0;
 
    pepProteins.insert(std::make_pair(pep, groupProtNames_[ indices[k] ].getVector()));
  }
  
  if (severedProteins_.size()!=0) {
    pepProteins.insert(std::make_pair(1.0,severedProteins_.getVector()));
  }
 
  return;
}

void GroupPowerBigraph::getProteinProbsAndNames(std::vector<std::vector<std::string> > &names, std::vector<double> &probs) const {
  names.clear();
  probs.clear();
  
  Array<double> sorted = probsPresentProteins_;
  Array<int> indices = sorted.sort();
  for (int k=0; k<sorted.size(); k++) {
    double pep = (1.0 - sorted[k]);
    if (pep <= 0.0) pep = 0.0;
    if (pep >= 1.0) pep = 1.0;
    names.push_back(groupProtNames_[ indices[k] ].getVector());
    probs.push_back(pep);
  }
  
  if (severedProteins_.size() != 0) {
    names.push_back(severedProteins_.getVector());
    probs.push_back(1.0);
  }
  return;
}


double GroupPowerBigraph::getLogNumberStates() const {
  double total = 0;
  for (int k=0; k<subgraphs_.size(); k++) {
    total = Numerical::logAdd(total, subgraphs_[k].logNumberOfConfigurations());
  }
  // add one because each peptide needs to be estimated once
  return total + 1;
}

Array<BasicBigraph> GroupPowerBigraph::iterativePartitionSubgraphs(BasicBigraph & bb, double newPeptideThreshold) {
  bb.setPeptideThreshold(newPeptideThreshold);
  bb.prune();
  severedProteins_.append( bb.severedProteins );

  Array<BasicBigraph> preResult = bb.partitionSections();
  Array<BasicBigraph> result;

  for (int k = 0; k < preResult.size(); k++) {
    BasicGroupBigraph bgb = BasicGroupBigraph(preResult[k], noClustering_/*,trivialGrouping_*/);
    double logNumConfig = bgb.logNumberOfConfigurations();
    if ( logNumConfig > LOG_MAX_ALLOWED_CONFIGURATIONS && 
         log2(bgb.PSMsToProteins.size())+log2(bgb.getOriginalN()[0].size+1) <= LOG_MAX_ALLOWED_CONFIGURATIONS ) {
      double newThresh = 1.25*(newPeptideThreshold + 1e-6);
      Array<BasicBigraph> completelyFragmented = iterativePartitionSubgraphs(preResult[k], newThresh);
      result.append( completelyFragmented );
    } else if (logNumConfig > LOG_MAX_ALLOWED_CONFIGURATIONS) {
      // the graph cannot become pruned to the desired efficiency;
      // prune as much as possible
      double largest = Vector(preResult[k].PSMsToProteins.weights).max();
      Array<BasicBigraph> completelyFragmented = iterativePartitionSubgraphs(preResult[k], largest);
      result.append( completelyFragmented );
    } else {
      // the graph is already pruned to the desired degree
      result.add( preResult[k] );
    }
  }
  
  return result;
}

void GroupPowerBigraph::read(Scores* fullset) {
  BasicBigraph basicBigraph(psmThreshold_, peptideThreshold_, proteinThreshold_, 
                  peptidePrior_);
  basicBigraph.read(fullset, addPeptideDecoyLabel_);  
  initialize(basicBigraph);
}

void GroupPowerBigraph::read(istream& is) {
  BasicBigraph basicBigraph(psmThreshold_, peptideThreshold_, proteinThreshold_, 
                  peptidePrior_);
  basicBigraph.read(is, addPeptideDecoyLabel_);
  initialize(basicBigraph);
}

void GroupPowerBigraph::initialize(BasicBigraph& basicBigraph) {
  severedProteins_ = Array<string>();
  if (noPartitioning_) {
    basicBigraph.prune();
    severedProteins_.append(basicBigraph.severedProteins);
    subgraphs_ = Array<BasicGroupBigraph>(1, BasicGroupBigraph(basicBigraph, noClustering_, trivialGrouping_));
  } else {
    Array<BasicBigraph> subBasic;
    subBasic = Array<BasicBigraph>();
    
    if (noPruning_)
      subBasic = iterativePartitionSubgraphs(basicBigraph, -1);
    else
      subBasic = iterativePartitionSubgraphs(basicBigraph, peptideThreshold_);

    subgraphs_ = Array<BasicGroupBigraph>(subBasic.size());

    for (int k = 0; k < subBasic.size(); k++) {
      subgraphs_[k] = BasicGroupBigraph(subBasic[k], noClustering_, trivialGrouping_);
    }
  }
  getGroupProtNames();
}

ostream & operator <<(ostream & os, pair<double,double> rhs) {
  os << "("<< rhs.first << ", " << rhs.second << ")";
  return os;
}

Array<string> GroupPowerBigraph::peptideNames() const {
  Array<string> pepNames;
  for (int k = 0; k < subgraphs_.size(); k++) {
    pepNames.append(subgraphs_[k].PSMsToProteins.names);
  }
  return pepNames;
}

pair<Array<Array<string> >, Array<double> > GroupPowerBigraph::getDescendingProteinsAndWeights() const {
  Array<double> sorted = probsPresentProteins_;
  for (int k = 0; k < sorted.size(); k++) {
    double tmp = (double) sorted[k];
    if (std::isnan(tmp)) {
      cerr << "error: found nan in GroupPowerBigraph::getDescendingProteinWeights" << endl;
      cerr << "\tprotein " << groupProtNames_[k] << endl;
    }
  }
  
  Array<int> indices = sorted.sort();

  Array<Array<string> > groupNames;
  Array<double> probabilities;
  for (int k = 0; k < sorted.size(); k++) {
    probabilities.add(sorted[k]);
    groupNames.add(groupProtNames_[ indices[k] ] );
  }
  return pair<Array<Array<string> >, Array<double> >(groupNames, probabilities);
}
