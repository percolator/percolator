// Written by Oliver Serang 2009
// see license for more information

#include "GroupPowerBigraph.h"

double GroupPowerBigraph::LOG_MAX_ALLOWED_CONFIGURATIONS = 18;

Array<double> GroupPowerBigraph::proteinProbs(const GridModel & myGM)
{
  Array<double> result;
  
  for (int k=0; k<subgraphs.size(); k++)
    {
      //      cout << "Working on subgraph with size " << subgraphs[k].proteinsToPSMs.size() << endl;

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
  //  cout << "\tGetting logLike(alpha, beta | D) with " << (Model)myGM << endl;

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
  cerr << "\nProtein level probabilities:\n";
  for (int k=0; k<sorted.size(); k++)
  {
    cerr << sorted[k] << " " << groupProtNames[ indices[k] ] << endl;
  }
  if(severedProteins.size()!=0)
    cerr << "0.0 " << severedProteins << endl;
}

void GroupPowerBigraph::readFromMCMC(istream & graph, istream & pepProph)
{
  BasicBigraph bb;
  bb.readFromMCMC(graph, pepProph);

  //cout << "Separating cleverly ;)" << endl;

  numberClones = 0;
  severedProteins = Array<string>();

#ifndef NOSEPARATE
#ifndef NOPRUNE
  Array<BasicBigraph> subBasic = iterativePartitionSubgraphs(bb, 0.0);
#else
  Array<BasicBigraph> subBasic = iterativePartitionSubgraphs(bb, -1.0);  
#endif
#else
  Array<BasicBigraph> subBasic = Array<BasicBigraph>(1, bb);
#endif

  subgraphs = Array<BasicGroupBigraph>(subBasic.size());

  for (int k=0; k<subBasic.size(); k++)
    {
      subgraphs[k] = BasicGroupBigraph(subBasic[k]);
    }

  initialize();
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

  /***
  bool bigFlag = false;
  if ( newPeptideThreshold > 1.0 )
    {
      bb.displayDotty("Trying");
      bigFlag = true;
    }
  ***/

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
	  //	  cout << LOG_MAX_ALLOWED_CONFIGURATIONS << " < " << logNumConfig << endl;

	  //	  double smallest = Vector(preResult[k].PSMsToProteins.weights).min();
	  double newThresh = 1.25*(newPeptideThreshold + 1e-6);

	  //	  cerr << "\tworking on subdividing a group of size " << logNumConfig << " with threshold " << newThresh << endl; 

	  Array<BasicBigraph> completelyFragmented = iterativePartitionSubgraphs(preResult[k], newThresh);

	  result.append( completelyFragmented );
	}
      else
	{
	  //	  cout << LOG_MAX_ALLOWED_CONFIGURATIONS << " >= " << logNumConfig << endl;

	  //	  if ( logNumConfig > 3 )
	  //	    	    cerr << "\tadding section with size " << logNumConfig << endl;
	  result.add( preResult[k] );
	}

      //      if ( bigFlag ) 
      //	preResult[k].displayDotty("Sub");
    }

  return result;
}

// Mattia Tomasoni
void GroupPowerBigraph::read(Scores& fullset){
	//  cout << "Reading GroupPowerBigraph" << endl;
	//  cout << "\tReading BasicBigraph" << endl;
	BasicBigraph bb;
	bb.read(fullset);
	speedUp(bb);
}

void GroupPowerBigraph::read(istream & is)
{
	//  cout << "Reading GroupPowerBigraph" << endl;
	//  cout << "\tReading BasicBigraph" << endl;
	BasicBigraph bb;
	is >> bb;
	speedUp(bb);
}

#ifndef NOSEPARATE
void GroupPowerBigraph::speedUp(BasicBigraph& bb)
{
  numberClones = 0;
  severedProteins = Array<string>();

  //  cout << "Starting partition..." << endl;

  #ifndef NOPRUNE
  Array<BasicBigraph> subBasic = iterativePartitionSubgraphs(bb, 0.0);
  #else
  Array<BasicBigraph> subBasic = iterativePartitionSubgraphs(bb, -1);
  #endif
  
  //  cout << "Finished partitioning into " << subBasic.size() << " partitions" << endl;

  // useful for debugging the iterative partition for small graphs
  //  for (int cc=0; cc<subBasic.size(); cc++)
  //    {
  //subBasic[cc].displayDotty("Test");
  //    }

  //  cout << "Number of clones was " << numberClones << endl << endl;

  //  cout << "\tPartition complete" << endl;

  //  cout << "Total number of proteins is " << bb.proteinsToPSMs.size() << endl;
  //  cout << "Total number of PSMs is " << bb.PSMsToProteins.size() << endl;

  //  exit(0);

  subgraphs = Array<BasicGroupBigraph>(subBasic.size());

  for (int k=0; k<subBasic.size(); k++)
    {
      subgraphs[k] = BasicGroupBigraph(subBasic[k]);
    }

  initialize();
  //  cout << "Initialized" << endl;
}

#else

 // for brute force: put into one subgraph
 void GroupPowerBigraph::speedUp()
 {
 bb.prune();

 numberClones = bb.numberClones;

 subgraphs = Array<BasicGroupBigraph>(1, BasicGroupBigraph(bb) );
 initialize();
 }
#endif

void GroupPowerBigraph::initialize()
{
  getGroupProtNames();
  //  getSumLikelihoodOverAllAlphaBetaGivenD();

  //  cout << "Initialized" << endl;
}

void GroupPowerBigraph::gridScan() const
{
  GridModel local = gm;

  Array<Array<double> > mat( local.maxRowIndex(), Array<double>( local.maxColIndex(), -1) );
  double bestProb = -1;
  Model bestModel;
  for (local.start(); local.inRange(); local.advance())
    {
      double prob = probabilityAlphaBetaGivenD(local);
      //      double prob = logLikelihoodAlphaBetaGivenD(local);

      if ( bestProb < prob )
	{
	  bestProb = prob;
	  bestModel = Model(local);

	  //	  cout << "\tFound new best " << bestProb << endl << bestModel << endl;
	}

      mat[ local.getRowIndex() ][ local.getColIndex() ] = prob;
      //      mat[ local.getRowIndex() ][ local.getColIndex() ] = unprunedProbabilityAlphaBetaGivenD(local);
    }

  cerr << "Best model: " << bestModel << endl;
  cerr << "\thad prob: " << bestProb << endl;

  Matrix(mat).displayMatrix();
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
