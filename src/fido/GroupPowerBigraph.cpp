// Written by Oliver Serang 2009
// see license for more information

#include "GroupPowerBigraph.h"

double GroupPowerBigraph::LOG_MAX_ALLOWED_CONFIGURATIONS = 18;

Array<double> GroupPowerBigraph::proteinProbs(const GridModel & myGM)
{
  Array<double> result;
  
  for (int k=0; k<subgraphs.size(); k++)
    {
      subgraphs[k].getProteinProbs(myGM);
      result.append( subgraphs[k].proteinProbabilities() );
    }

  return result;
}

void GroupPowerBigraph::getProteinProbs()
{
  GridModel local = gm;
  local.start();
  probabilityR = proteinProbs(local);
}

double GroupPowerBigraph::likelihoodAlphaBetaGivenD(const GridModel & myGM) const
{
  return pow(2.0, logLikelihoodAlphaBetaGivenD(myGM) );

  double prod = 1.0;

  for (int k=0; k<subgraphs.size(); k++)
    {
      prod *= subgraphs[k].likelihoodAlphaBetaGivenD(myGM);
    }

  return prod / pow( 1-myGM.spontaneousEmission() , numberClones);
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

double GroupPowerBigraph::probabilityAlphaBetaGivenD(const GridModel & myGM) const
{
  // using cached functor
  return pow(2.0, logLikelihoodAlphaBetaGivenD(myGM) - sumLogLikelihoodOverAllAlphaBetaCachedFunctor(myGM, this) );
}

double GroupPowerBigraph::sumLogLikelihoodOverAllAlphaBeta(const GridModel & myGM) const
{
  GridModel local = myGM;

  double result = 0.0;
  bool starting = true;

  for (local.start(); local.inRange(); local.advance())
    {
      double logLike = logLikelihoodAlphaBetaGivenD(local);
     
      if ( starting ) 
	{
	  starting = false;
	  result = logLike;
	}
      else
	{
	  result = Numerical::logAdd(result, logLike);
	}
    }

  return result;
}

void GroupPowerBigraph::getProteinProbsOverAllAlphaBeta()
{
  probabilityR = proteinProbsOverAllAlphaBeta();
}

Array<double> GroupPowerBigraph::proteinProbsOverAllAlphaBeta()
{
  Vector cumulative;

  GridModel local = gm;
  for ( local.start(); local.inRange(); local.advance() )
    {
      double prob = probabilityAlphaBetaGivenD(local);

      // hack for efficiency
      if ( prob > 1e-5 )
	{
	  Array<double> protProbs = proteinProbs( local );
	  Vector posteriorsForCurrentAlphaBeta = prob * Vector( protProbs );

	  if ( cumulative.size() == 0 )
	    cumulative = posteriorsForCurrentAlphaBeta;
	  else
	    cumulative += posteriorsForCurrentAlphaBeta;
	}

    }

  return cumulative.unpack();
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
  Array<double> sorted = probabilityR;
  Array<int> indices = sorted.sort();
  cout << "\nProtein level probabilities:\n";
  for (int k=0; k<sorted.size(); k++)
  {
    cout << sorted[k] << " " << groupProtNames[ indices[k] ] << endl;
  }
  if(severedProteins.size()!=0)
    cout << "0.0 " << severedProteins << endl;
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
    std::vector<std::string> proteins;
    for ( int j = 0; j < groupProtNames[ indices[k] ].size(); j++)
    {
      std::string protein = groupProtNames[ indices[k] ][j];
      proteins.push_back(protein);
    }
    pepProteins.insert(std::make_pair<double,std::vector<std::string> >(pep,proteins));
  }
  if(severedProteins.size()!=0)
  {
    double pep = 1.0;
    std::vector<std::string> proteins;
    for(int i = 0; i < severedProteins.size(); i++)
    {
      std::string protein = severedProteins[i];
      proteins.push_back(protein);
    }
    pepProteins.insert(std::make_pair<double,std::vector<std::string> >(pep,proteins));
  }
    
  return pepProteins;
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
  //  cerr << "Iter partition... @ threshold = " << newPeptideThreshold << endl;

  bb.PeptideThreshold = newPeptideThreshold;
  bb.prune();
  severedProteins.append( bb.severedProteins );
  numberClones += bb.numberClones;

  Array<BasicBigraph> preResult = bb.partitionSections();
  Array<BasicBigraph> result;

  //  cerr << "Using threshold " << newPeptideThreshold << endl;
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

  //TODO when noseparate it does nor work
  if(!noseparate)
  {
    numberClones = 0;
    severedProteins = Array<string>();

    Array<BasicBigraph> subBasic;
    subBasic = Array<BasicBigraph>();
    
    if(!noprune)
      subBasic = iterativePartitionSubgraphs(bb, 0.0);
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
    bb.prune();

    numberClones = bb.numberClones;

    subgraphs = Array<BasicGroupBigraph>(1, BasicGroupBigraph(bb) );
    
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
  this->alpha = alpha;
  this->beta = beta;
  this->gamma = gamma;
  gm.setalphaRange(RealRange(alpha,1,alpha));
  gm.setbetaRange(RealRange(beta,1,beta));
  gm.setGamma(gamma);
 
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