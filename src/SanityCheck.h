/*******************************************************************************
 * Percolator unofficial version
 * Copyright (c) 2006-7 University of Washington. All rights reserved.
 * Written by Lukas Käll (lukall@u.washington.edu) in the 
 * Department of Genome Science at the University of Washington. 
 *
 * $Id: SanityCheck.h,v 1.3 2008/04/19 21:40:33 lukall Exp $
 *******************************************************************************/
#ifndef SANITYCHECK_H_
#define SANITYCHECK_H_
class Scores;
class Normalizer;

class SanityCheck
{
public:
  SanityCheck();
  virtual ~SanityCheck();
    
  void readWeights(istream & weightStream, double * w);
  int getInitDirection(Scores * testset,Scores * trainset, Normalizer * pNorm,double *w,double test_fdr);
  virtual bool validateDirection(double* w);
  void resetDirection(double* w);

  static void setInitWeightFN(string fn) {initWeightFN=fn;}
  static void setOverrule(bool orl) {overRule=orl;}
protected:
  virtual void getDefaultDirection(double *w);
  int initPositives;
  double fdr;
  static bool overRule;
  static string initWeightFN;
  Scores *pTestset, *pTrainset;
};

#endif /*SANITYCHECK_H_*/
