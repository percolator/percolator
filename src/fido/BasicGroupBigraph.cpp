// Written by Oliver Serang 2009
// see license for more information

#include "BasicGroupBigraph.h"

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

  //  cout << "\tL(N=n | D) for n = " << n << " is ";
  for (int k=0; k<PSMsToProteins.size(); k++)
    {
      double probEGivenD = PSMsToProteins.weights[k];
      double probEGivenN = probabilityEEpsilonGivenN(m, k, n);

      //      double probE = probabilityE(m)[k];

      // using cached functor
      //      double probE = probabilityECachedFunctor(m, this)[ k ];

      double probE = PeptideProphetPrior;
      //      double probE = .5;

      double termE = probEGivenD / probE * probEGivenN;
      double termNotE = (1-probEGivenD) / (1-probE) * (1-probEGivenN);
      double term = termE + termNotE;
      
      prod *= term;
    }

  //  cout << prod << endl;

  return prod;
}

double BasicGroupBigraph::logLikelihoodNGivenD(const Model & m, const Array<Counter> & n) const
{
  double logProd = 0.0;

  //  cout << "\tL(N=n | D) for n = " << n << " is ";
  for (int k=0; k<PSMsToProteins.size(); k++)
    {
      double probEGivenD = PSMsToProteins.weights[k];
      double probEGivenN = probabilityEEpsilonGivenN(m, k, n);

      //      double probE = probabilityE(m)[k];

      // using cached functor
      //      double probE = probabilityECachedFunctor(m, this)[ k ];

      double probE = PeptideProphetPrior;
      //      double probE = .5;

      double termE = probEGivenD / probE * probEGivenN;
      double termNotE = (1-probEGivenD) / (1-probE) * (1-probEGivenN);
      double term = termE + termNotE;
      
      logProd += log2(term);
    }

  //  cout << prod << endl;

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

  //  cout << "\tGetting logNumber configs for " << originalN << endl;
  for (int k=0; k<originalN.size(); k++)
    {
      result += log2(originalN[k].size+1);
    }

  //  cout << "\t\twas = " << result << endl;

  return result;
}

double BasicGroupBigraph::probabilityNNu(const Model & m, const Counter & nNu) const
{
  return m.probabilityProteins(nNu.size, nNu.state);
}

double BasicGroupBigraph::probabilityNGivenD(const Model & m, const Array<Counter> & n) const
{
  //  return likelihoodNGivenD(m, n) * probabilityN(m, n) / likelihoodConstant(m);

  // using cached functor
  //  cout << "Calling functor" << endl;
  //  cout << likelihoodConstantCachedFunctor.name << endl;
  //  cout << "correct object address is " << this << endl;

  // working version, but numerically unstable
  //  double like = likelihoodNGivenD(m, n) * probabilityN(m, n) / likelihoodConstantCachedFunctor(m, this);

  // log version
  double logLike= logLikelihoodNGivenD(m,n) + log2(probabilityN(m,n)) - logLikelihoodConstantCachedFunctor(m,this);

  /***
  cout << "L = " << like << endl;
  cout << "2^logL = " << pow(2.0, logLike) << endl << endl;

  cout << "L(N|D) = " << likelihoodNGivenD(m,n) << endl;
  cout << "2^logL(N|D) = " << pow( 2.0, logLikelihoodNGivenD(m,n) ) << endl << endl;
  ***/

  //return like;
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
      //      cout << "\tL and p are: " << L << ", " << p << endl;
      //      cout << "\t\tterm is " << term << endl;

      //      cout << "\tif you'd used logL: " << logLikelihoodNGivenD(m,n) << endl;
      result += term;
    }

  //  if ( isinf( result ) )
  //    print();
  //    displayDotty("Checker");

  //  cout << "\tWith result " << result << endl;
  return result;
}

Array<double> BasicGroupBigraph::probabilityEGivenD(const Model & m)
{
  Array<Counter> n = originalN;

  //  cout << "originalN = " << originalN << endl;

  Vector result;

  for (Counter::start(n); Counter::inRange(n); Counter::advance(n))
    {
      // NOTE: should this be logged?
      Vector term = pow(2.0, logLikelihoodNGivenD(m, n) - logLikelihoodConstantCachedFunctor(m, this) ) * Vector( eCorrection(m, n) );

      //      Vector term = likelihoodNGivenD(m, n)/likelihoodConstantCachedFunctor(m, this) * Vector( eCorrection(m, n) );

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

      //      double probE = probabilityE(m)[k];

      // using cached functor
      //      double probE = probabilityECachedFunctor(m, this)[ k ];

      double probE = PeptideProphetPrior;
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

  //  cout << "originalN = " << originalN << endl;

  Vector result;

  for (Counter::start(n); Counter::inRange(n); Counter::advance(n))
    {
      /***
      cout << "Current iter: " << endl;
      cout << '\t' << n << endl;
      cout << '\t' << probabilityNGivenD(m,n);
      cout << '\t' << probabilityRGivenN(n) << endl << endl;
      ***/

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
  //  cout << m << endl;
  probabilityR = probabilityRGivenD(m);
  //  cout << "BasicGroupBigraph::probR = " << probabilityR << endl;
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

double BasicGroupBigraph::likelihoodAlphaBetaGivenD(const GridModel & gm) const
{
  return likelihoodConstantCachedFunctor(gm, this);
}

double BasicGroupBigraph::logLikelihoodAlphaBetaGivenD(const GridModel & gm) const
{
  return logLikelihoodConstantCachedFunctor(gm, this);
}

