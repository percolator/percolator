// Written by Oliver Serang 2009
// see license for more information

#include "GroupPowerBigraph.h"

double GroupPowerBigraph::LOG_MAX_ALLOWED_CONFIGURATIONS = 18;

GroupPowerBigraph::~GroupPowerBigraph()
{
  /*delete gm;
  delete zeroChecker;
  FreeAll(severedProteins);
  FreeAll(probabilityR);
  for(unsigned i = 0 ; i < groupProtNames.size(); i++)
  {
    FreeAll(groupProtNames[i]);
  }
  FreeAll(groupProtNames);
  for(unsigned i = 0; i < subgraphs.size(); i++)
  {
    delete subgraphs[i];
  }
  FreeAll(subgraphs);*/
}

Array<double> GroupPowerBigraph::proteinProbs()
{
  Array<double> result;
  
  for (int k=0; k<subgraphs.size(); k++)
    {
      subgraphs[k].getProteinProbs(gm);
      result.append( subgraphs[k].proteinProbabilities() );
    }

  return result;
}

void GroupPowerBigraph::getProteinProbs()
{
  probabilityR = proteinProbs();
}

double GroupPowerBigraph::logLikelihoodAlphaBetaGivenD(const GridModel & myGM) const
{

  double sum = 0.0;

  for (int k=0; k<subgraphs.size(); k++)
    {
      sum += subgraphs[k].logLikelihoodAlphaBetaGivenD(myGM);
    }

  return sum - numberClones * log2( 1-myGM.spontaneousEmission() );
}

void GroupPowerBigraph::getGroupProtNames()
{
  int k,j;
  for (k=0; k<subgraphs.size(); k++)
    {
      const BasicGroupBigraph & bgb = subgraphs[k];

      for (j=0; j<bgb.proteinsToPSMs.size(); j++)
	{
	  groupProtNames.add( bgb.proteinGroupNames()[j] );
	}
    }
}

void GroupPowerBigraph::printProteinWeights() const
{
  //NOTE this is priting out PEPs not posterior probabilities
  Array<double> sorted = probabilityR;
  Array<int> indices = sorted.sort();
  for (int k=0; k<sorted.size(); k++)
  {
    cout << double(1 - sorted[k]) << " " << groupProtNames[ indices[k] ] << endl;
  }
  if(severedProteins.size()!=0)
    cout << "1.0 " << severedProteins << endl;
}

/*return a map of PEPs and their respectives proteins */
std::multimap<double, std::vector<std::string> > GroupPowerBigraph::getProteinProbsPercolator() const
{
  std::multimap<double, std::vector<std::string> > pepProteins;
  Array<double> sorted = probabilityR;
  Array<int> indices = sorted.sort();
  for (int k=0; k<sorted.size(); k++)
  {
    double pep = (1.0 - (double)sorted[k]);
    if(pep < 0.0)pep = 0.0;
    if(pep > 1.0)pep = 1.0;
 
    pepProteins.insert(std::make_pair<double,std::vector<std::string> >(pep,groupProtNames[ indices[k] ].getVector()));
  }
  if(severedProteins.size()!=0)
  {
    pepProteins.insert(std::make_pair<double,std::vector<std::string> >(1.0,severedProteins.getVector()));
  }
 
  return pepProteins;
}

void GroupPowerBigraph::getProteinProbsAndNames(std::vector<std::vector<std::string> > &names, std::vector<double> &probs) const
{
  names.clear();
  probs.clear();
  Array<double> sorted = probabilityR;
  Array<int> indices = sorted.sort();
  for (int k=0; k<sorted.size(); k++)
  {
    double pep = (1.0 - (double)sorted[k]);
    if(pep <= 0.0)pep = 0.0;
    if(pep >= 1.0)pep = 1.0;
    names.push_back(groupProtNames[ indices[k] ].getVector());
    probs.push_back(pep);
  }
  if(severedProteins.size()!=0)
  {
    names.push_back(severedProteins.getVector());
    probs.push_back(1.0);
  }
  return;
}


double GroupPowerBigraph::getLogNumberStates() const
{
  double total = 0;
  for (int k=0; k<subgraphs.size(); k++)
    {
      total = Numerical::logAdd(total, subgraphs[k].logNumberOfConfigurations());
    }
  // add one because each peptide needs to be estimated once
  return total + 1;
}

Array<BasicBigraph> GroupPowerBigraph::iterativePartitionSubgraphs(BasicBigraph & bb, double newPeptideThreshold )
{

  bb.PeptideThreshold = newPeptideThreshold;
  bb.prune();
  severedProteins.append( bb.severedProteins );
  numberClones += bb.numberClones;

  Array<BasicBigraph> preResult = bb.partitionSections();
  Array<BasicBigraph> result;

  for (int k=0; k<preResult.size(); k++)
    {
      double logNumConfig = BasicGroupBigraph( preResult[k] ).logNumberOfConfigurations();
      if ( logNumConfig > LOG_MAX_ALLOWED_CONFIGURATIONS )
	{
	  double newThresh = 1.25*(newPeptideThreshold + 1e-6);
	  Array<BasicBigraph> completelyFragmented = iterativePartitionSubgraphs(preResult[k], newThresh);
	  result.append( completelyFragmented );
	}
      else
	{
	  result.add( preResult[k] );
	}
    }

  return result;
}

void GroupPowerBigraph::read(Scores* fullset){
  
  BasicBigraph bb;
  bb.read(fullset);

  //TODO when noseparate active it does not work
  if(!noseparate)
  {
    numberClones = 0;
    severedProteins = Array<string>();

    Array<BasicBigraph> subBasic;
    subBasic = Array<BasicBigraph>();
    
    if(!noprune)
      subBasic = iterativePartitionSubgraphs(bb, PeptideThreshold);
    else
      subBasic = iterativePartitionSubgraphs(bb, -1);

    subgraphs = Array<BasicGroupBigraph>(subBasic.size());

    for (int k=0; k<subBasic.size(); k++)
      {
	subgraphs[k] = BasicGroupBigraph(subBasic[k],groupProteins);
      }
  }
  else
  {
    if(!noprune)
    {
      bb.prune();
      numberClones = bb.numberClones;
    }
    else
    {
      numberClones = 0;
    }
    severedProteins = Array<string>();
    subgraphs = Array<BasicGroupBigraph>(1, BasicGroupBigraph(bb,groupProteins));
  }
  
  initialize();

}



void GroupPowerBigraph::initialize()
{
  getGroupProtNames();
}


ostream & operator <<(ostream & os, pair<double,double> rhs)
{
  os << "("<< rhs.first << ", " << rhs.second << ")";
  return os;
}

void GroupPowerBigraph::outputPivdo(ostream & os) const
{
  for (int k=0; k<subgraphs.size(); k++)
    PivdoSplitter(subgraphs[k]).outputPivdo(os);
}

void GroupPowerBigraph::setAlphaBetaGamma(double alpha, double beta, double gamma)
{
  gm.setAlphaBetaGamma(alpha,beta,gamma);
}

Array<string> GroupPowerBigraph::peptideNames() const
{
  Array<string> pepNames;

  for (int k=0; k<subgraphs.size(); k++)
    {
      pepNames.append( subgraphs[k].PSMsToProteins.names );
    }

  return pepNames;
}

pair<Array<Array<string> >, Array<double> > GroupPowerBigraph::getDescendingProteinsAndWeights() const
{
  Array<double> sorted = probabilityR;
  int k;
  for (k=0; k<sorted.size(); k++)
    {
      double tmp = (double) sorted[k];
      if ( isnan( tmp ) )
	{
	  cerr << "error: found nan in GroupPowerBigraph::getDescendingProteinWeights" << endl;
	  cerr << "\tprotein " << groupProtNames[k] << endl;
	}
    }
  
  Array<int> indices = sorted.sort();

  Array<Array<string> > groupNames;
  Array<double> probabilities;
  for (k=0; k<sorted.size(); k++)
    {
      probabilities.add(sorted[k]);
      groupNames.add(groupProtNames[ indices[k] ] );
    }
  return pair<Array<Array<string> >, Array<double> >(groupNames, probabilities);
}

Array<std::string> GroupPowerBigraph::getSeveredProteins()
{
  return severedProteins;
}