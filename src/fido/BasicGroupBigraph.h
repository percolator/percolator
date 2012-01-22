// Written by Oliver Serang 2009
// see license for more information

#ifndef _BasicGroupBigraph_H
#define _BasicGroupBigraph_H

#include "ProteinIdentifier.h"
#include "ReplicateIndexer.h"
#include "BasicBigraph.h"
#include "Model.h"
#include "Cache.h"

//TODO parameter??
//#define TRUE_BRUTE
#ifdef TRUE_BRUTE
#define NOCACHE
#endif

class Counter
{
 public:
  Counter()
    {
      size = 0;
      start();
    }
  Counter(int s)
    {
      size = s;
      start();
    }
  int size;
  int state;

  void start()
  {
    state = 0;
  }
  bool inRange() const
  {
    return state <= size;
  }
  void advance()
  {
    if ( size == 0 )
      cerr << "Warning: advancing group with size 0" << endl;

    state++;
  }

  friend ostream & operator <<(ostream & os, const Counter & c)
  {
    return os << " Counter(" << c.state << "/" << c.size << ") ";
  }
  
  static void start(Array<Counter> & cA)
  {
    for (int k=0; k<cA.size(); k++)
      cA[k].start();
  }

  static bool inRange(Array<Counter> & cA)
  {
    return cA.back().inRange();
  }

  static void advance(Array<Counter> & cA)
  {
    for (int k=0; k<cA.size(); k++)
      {
	cA[k].advance();
	if ( ! cA[k].inRange() && k != cA.size() - 1 )
	  {
	    // do not reset the last one so that you know when it is not inRange
	    cA[k].start();
	  }
	else
	  {
	    // no need to carry, so you can stop
	    break;
	  }
      }
  }
};


class BasicGroupBigraph : public BasicBigraph
{
  // hack should be protected
 public:
  // protected:

  Array<Counter> originalN;
  Array<Array<string> > groupProtNames;
  Array<double> probabilityR;

  // protected construction functions
  void groupProteins();
  void trivialGroupProteins();
  void groupProteinsBy(const Array<Set> & groups);
 
  // protected utility functions
  int numberAssociatedProteins(int indexEpsilon) const;
  int numberActiveAssociatedProteins(int indexEpsilon, const Array<Counter> & n) const;
 
  // protected mathematical functions
  double probabilityEEpsilon(const Model & m, int indexEpsilon) const;
  Array<double> probabilityE(const Model & m) const;
  double probabilityEEpsilonGivenActiveAssociatedProteins(const Model & m, int active) const;
  double probabilityNumberAssociatedProteins(const Model & m, int total, int active) const;
  double probabilityEEpsilonGivenN(const Model & m, int indexEpsilon, const Array<Counter> & n) const;

  double likelihoodNGivenD(const Model & m, const Array<Counter> & n) const;
  double logLikelihoodNGivenD(const Model & m, const Array<Counter> & n) const;

  double probabilityN(const Model & m, const Array<Counter> & n) const;
  double probabilityNNu(const Model & m, const Counter & nNu) const;
  double probabilityNGivenD(const Model & m, const Array<Counter> & n) const;

  double logLikelihoodConstant(const Model & m) const;
  double likelihoodConstant(const Model & m) const;

  Array<double> probabilityRGivenD(const Model & m);
  Array<double> probabilityRGivenN(const Array<Counter> & n);
  double probabilityRRhoGivenN(int indexRho, const Array<Counter> & n);

  Array<double> probabilityEGivenD(const Model & m);
  Array<double> eCorrection(const Model & m, const Array<Counter> & n);

  // for unknown alpha, beta
  double probabilityEEpsilonOverAllAlphaBeta(const GridModel & gm, int indexEpsilon) const;
  Array<double> probabilityEOverAllAlphaBeta(const GridModel & gm) const;

  // cache functors
  // note that these will need to be updated if the object is copied
  LastCachedMemberFunction<BasicGroupBigraph, double, Model> likelihoodConstantCachedFunctor;
  LastCachedMemberFunction<BasicGroupBigraph, Array<double>, Model> probabilityECachedFunctor;
  LastCachedMemberFunction<BasicGroupBigraph, Array<double>, GridModel> probabilityEOverAllAlphaBetaCachedFunctor;
  LastCachedMemberFunction<BasicGroupBigraph, double, Model> logLikelihoodConstantCachedFunctor;

 public:
  // for partitioning
  double logNumberOfConfigurations() const;

  void getProteinProbs(const Model & m);
  void printProteinWeights() const;

  double logLikelihoodAlphaBetaGivenD(const GridModel & gm) const;
  double likelihoodAlphaBetaGivenD(const GridModel & gm) const;

  const Array<double> & proteinProbabilities() const
  {
    return probabilityR;
  }
  const Array<Array<string> > & proteinGroupNames() const
  {
    return groupProtNames;
  }

 BasicGroupBigraph(bool __groupProtein = true) :
  likelihoodConstantCachedFunctor( & BasicGroupBigraph::likelihoodConstant, "likelihoodConstant"), 
  probabilityECachedFunctor( & BasicGroupBigraph::probabilityE, "probabilityE"), 
  probabilityEOverAllAlphaBetaCachedFunctor( & BasicGroupBigraph::probabilityEOverAllAlphaBeta, "probabilityEOverAllAlphaBeta"), 
  logLikelihoodConstantCachedFunctor( & BasicGroupBigraph::logLikelihoodConstant, "logLikelihoodConstant"),
  groupProtein(__groupProtein)
  {
  }

  void read(Scores* fullset)
   {
     BasicBigraph::read(fullset);
     if(groupProtein)
	groupProteins();
     else
        trivialGroupProteins();
   }

 BasicGroupBigraph(const BasicBigraph & rhs,bool __groupProtein = true) :
 BasicBigraph(rhs),   
 likelihoodConstantCachedFunctor( & BasicGroupBigraph::likelihoodConstant, "likelihoodConstant"), 
 probabilityECachedFunctor( & BasicGroupBigraph::probabilityE, "probabilityE"), 
 probabilityEOverAllAlphaBetaCachedFunctor( & BasicGroupBigraph::probabilityEOverAllAlphaBeta, "probabilityEOverAllAlphaBeta"), 
 logLikelihoodConstantCachedFunctor( & BasicGroupBigraph::logLikelihoodConstant, "logLikelihoodConstant"),
 groupProtein(__groupProtein)
  {
    if(groupProtein)
      groupProteins();
    else
      trivialGroupProteins();
  }
  
protected: 
  bool groupProtein;
};


#endif

