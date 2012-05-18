// Written by Oliver Serang 2009
// see license for more information

#include "BasicGroupBigraph.h"

BasicGroupBigraph::~BasicGroupBigraph()
{
  /*for(unsigned i = 0; i < groupProtNames.size(); i++)
  {
    FreeAll(groupProtNames[i]);
  }
  FreeAll(groupProtNames);
  FreeAll(probabilityR);
  for(unsigned i = 0; i < originalN.size(); i++)
  {
    delete originalN[i];
  }
  FreeAll(originalN);
  delete logLikelihoodConstantCachedFunctor;*/
}


void BasicGroupBigraph::printProteinWeights() const
{
  Array<double> sorted = proteinsToPSMs.weights;
  Array<int> indices = sorted.sort();

  for (int k=0; k<proteinsToPSMs.size(); k++)
    {
      cout << sorted[k] << " " << groupProtNames[ indices[k] ] << endl;
    }
}

void BasicGroupBigraph::trivialGroupProteins()
{
  groupProtNames = Array<Array<string> >(proteinsToPSMs.size());
  originalN = Array<Counter>(proteinsToPSMs.size());

  for (int k=0; k<proteinsToPSMs.size(); k++)
    {
      groupProtNames[k] = Array<string>(1, proteinsToPSMs.names[k]);
      originalN[k] = Counter(1);
    }
}

void BasicGroupBigraph::groupProteinsBy(const Array<Set> & groups)
{
  // remake the list of names and then remake the graph with them collapsed
  groupProtNames = Array<Array<string> > (groups.size());
  originalN = Array<Counter> (groups.size());

  int k;
  for (k=0; k<groups.size(); k++)
    {
      groupProtNames[k] = proteinsToPSMs.names[ groups[k] ];
      originalN[k] = Counter( groups[k].size() );
    }

  // remove all but the first of each group from the graph
  for (k=0; k<groups.size(); k++)
    {
      Set reps = groups[k].without( Set::SingletonSet(groups[k][0]) );
      for (Set::Iterator iter = reps.begin(); iter != reps.end(); iter++)
	{
	  disconnectProtein(*iter);
	}
    }

  reindex();
}

void BasicGroupBigraph::groupProteins()
{
  Array<Set> groups = ReplicateIndexer<Set>::replicates(Set::sumSetElements, proteinsToPSMs.associations );

  groupProteinsBy(groups);
}

int BasicGroupBigraph::numberAssociatedProteins(int indexEpsilon) const
{
  int tot = 0;

  const Set & s = PSMsToProteins.associations[indexEpsilon];
  for (int k=0; k<s.size(); k++)
    {
      tot += originalN[ s[k] ].size;
    }
  
  return tot;
}

int BasicGroupBigraph::numberActiveAssociatedProteins(int indexEpsilon, const Array<Counter> & n) const
{
  int tot = 0;

  const Set & s = PSMsToProteins.associations[indexEpsilon];
  for (int k=0; k<s.size(); k++)
    {
      tot += n[ s[k] ].state;
    }
  
  return tot;
}

Array<double> BasicGroupBigraph::probabilityE(const Model & m) const
{
  Array<double> result( PSMsToProteins.size() );

  for (int k=0; k<result.size(); k++)
    {
      result[k] = probabilityEEpsilon(m, k);
    }

  return result;
}

double BasicGroupBigraph::probabilityEEpsilon(const Model & m, int indexEpsilon) const
{
  double tot = 0.0;
  int a = numberAssociatedProteins(indexEpsilon);
  
  for (int k=0; k <= a; k++)
    {
      tot += probabilityEEpsilonGivenActiveAssociatedProteins(m, k) * probabilityNumberAssociatedProteins(m, a, k);
    }

  return tot;
}

double BasicGroupBigraph::probabilityEEpsilonGivenActiveAssociatedProteins(const Model & m, int active) const
{
  return 1 - m.probabilityNoEmissionFrom(active);
}

double BasicGroupBigraph::probabilityNumberAssociatedProteins(const Model & m, int total, int active) const
{
  return m.probabilityProteins(total, active);
}

double BasicGroupBigraph::probabilityEEpsilonGivenN(const Model & m, int indexEpsilon, const Array<Counter> & n) const
{
  int active = numberActiveAssociatedProteins(indexEpsilon, n);

  return probabilityEEpsilonGivenActiveAssociatedProteins(m, active);
}

double BasicGroupBigraph::likelihoodNGivenD(const Model & m, const Array<Counter> & n) const
{
  double prod = 1.0;
  for (int k=0; k<PSMsToProteins.size(); k++)
    {
      double probEGivenD = PSMsToProteins.weights[k];
      double probEGivenN = probabilityEEpsilonGivenN(m, k, n);
      double probE = PeptidePrior;
      // NOTE peptidePrior is 0.07
      double termE = probEGivenD / probE * probEGivenN;
      double termNotE = (1-probEGivenD) / (1-probE) * (1-probEGivenN);
      double term = termE + termNotE;
      
      prod *= term;
    }

  return prod;
}

double BasicGroupBigraph::logLikelihoodNGivenD(const Model & m, const Array<Counter> & n) const
{
  double logProd = 0.0;

  for (int k=0; k<PSMsToProteins.size(); k++)
    {
      double probEGivenD = PSMsToProteins.weights[k];
      double probEGivenN = probabilityEEpsilonGivenN(m, k, n);
      double probE = PeptidePrior;
      // NOTE peptidePrior is 0.07
      double termE = probEGivenD / probE * probEGivenN;
      double termNotE = (1-probEGivenD) / (1-probE) * (1-probEGivenN);
      double term = termE + termNotE;
      
      logProd += log2(term);
    }

  return logProd;
}

double BasicGroupBigraph::probabilityN(const Model & m, const Array<Counter> & n) const
{
  double prod = 1.0;

  for (int k=0; k<n.size(); k++)
    {
      prod *= probabilityNNu(m, n[k]);
    }

  return prod;
}

double BasicGroupBigraph::logNumberOfConfigurations() const
{
  double result = 0.0;

  for (int k=0; k<originalN.size(); k++)
    {
      result += log2(originalN[k].size+1);
    }

  return result;
}

double BasicGroupBigraph::probabilityNNu(const Model & m, const Counter & nNu) const
{
  return m.probabilityProteins(nNu.size, nNu.state);
}

double BasicGroupBigraph::probabilityNGivenD(const Model & m, const Array<Counter> & n) const
{
  // working version, but numerically unstable
  //  double like = likelihoodNGivenD(m, n) * probabilityN(m, n) / likelihoodConstantCachedFunctor(m, this);
  // log version
  double logLike= logLikelihoodNGivenD(m,n) + log2(probabilityN(m,n)) - logLikelihoodConstantCachedFunctor(m,this);
  return pow(2.0, logLike);
}

double BasicGroupBigraph::logLikelihoodConstant(const Model & m) const
{
  double result = 0.0;
  bool starting = true;

  Array<Counter> n = originalN;

  for (Counter::start(n); Counter::inRange(n); Counter::advance(n))
    {
      double L = logLikelihoodNGivenD(m, n);
      double p = log2(probabilityN(m, n));
      double logLikeTerm = L+p;

      if ( starting )
	{
	  starting = false;
	  result = logLikeTerm;
	}
      else
	{
	  result = Numerical::logAdd(result, logLikeTerm);
	}
    }

  return result;
}

double BasicGroupBigraph::likelihoodConstant(const Model & m) const
{
  double result = 0.0;

  Array<Counter> n = originalN;
  for (Counter::start(n); Counter::inRange(n); Counter::advance(n))
    {
      double L = likelihoodNGivenD(m, n);
      double p = probabilityN(m, n);
      double term = L*p;
      result += term;
    }
    
  return result;
}

Array<double> BasicGroupBigraph::probabilityEGivenD(const Model & m)
{
  Array<Counter> n = originalN;
  Vector result;

  for (Counter::start(n); Counter::inRange(n); Counter::advance(n))
    {
      // NOTE: should this be logged?
      Vector term = pow(2.0, logLikelihoodNGivenD(m, n) - logLikelihoodConstantCachedFunctor(m, this) ) * Vector( eCorrection(m, n) );
      //Vector term = likelihoodNGivenD(m, n)/likelihoodConstantCachedFunctor(m, this) * Vector( eCorrection(m, n) );

      if ( result.size() == 0 )
	{
	  result = term;
	}
      else
	{
	  result += term;
	}
    }

  return result.unpack();
}

Array<double> BasicGroupBigraph::eCorrection(const Model & m, const Array<Counter> & n)
{
  Array<double> result(PSMsToProteins.size());

  for (int k=0; k<result.size(); k++)
    {
      // NOTE: check this

      double probEGivenD = PSMsToProteins.weights[k];
      double probEGivenN = probabilityEEpsilonGivenN(m, k, n);
      double probE = PeptidePrior;
      // NOTE peptidePrior is 0.07
      //      double probE = .5;
      double termE = probEGivenD / probE * probEGivenN;
      double termNotE = (1-probEGivenD) / (1-probE) * (1-probEGivenN);
      double term = termE + termNotE;
    }

  return result;
}

Array<double> BasicGroupBigraph::probabilityRGivenD(const Model & m)
{
  Array<Counter> n = originalN;
  Vector result;

  for (Counter::start(n); Counter::inRange(n); Counter::advance(n))
    {

      Vector term = probabilityNGivenD(m, n) * Vector( probabilityRGivenN(n) );

      if ( result.size() == 0 )
	{
	  result = term;
	}
      else
	{
	  result += term;
	}
    }

  return result.unpack();
}

Array<double> BasicGroupBigraph::probabilityRGivenN(const Array<Counter> & n)
{
  Array<double> result(n.size());

  for (int k=0; k<result.size(); k++)
    {
      result[k] = probabilityRRhoGivenN(k, n);
    }

  return result;
}

double BasicGroupBigraph::probabilityRRhoGivenN(int indexRho, const Array<Counter> & n)
{
  const Counter & c = n[indexRho];

  return double(c.state) / c.size;
}

void BasicGroupBigraph::getProteinProbs(const Model & m)
{
  probabilityR = probabilityRGivenD(m);
}

double BasicGroupBigraph::probabilityEEpsilonOverAllAlphaBeta(const GridModel & gm, int indexEpsilon) const
{
  GridModel localModel( gm );

  double val = 0.0;
  for ( localModel.start(); localModel.inRange(); localModel.advance() )
    {
      val += probabilityEEpsilon(localModel, indexEpsilon);
    }

  val /= localModel.getCount();

  return val;
}

Array<double> BasicGroupBigraph::probabilityEOverAllAlphaBeta(const GridModel & gm) const
{
  Array<double> result( PSMsToProteins.size() );

  for ( int k=0; k<result.size(); k++)
    {
      result[k] = probabilityEEpsilonOverAllAlphaBeta(gm, k);
    }

  return result;
}

double BasicGroupBigraph::logLikelihoodAlphaBetaGivenD(const GridModel & gm) const
{
  return logLikelihoodConstantCachedFunctor(gm, this);
}
