/*******************************************************************************
 * Percolator unofficial version
 * Copyright (c) 2006-9 University of Washington. All rights reserved.
 * Written by Lukas Käll (lukall@u.washington.edu) in the
 * Department of Genome Sciences at the University of Washington.
 *
 * $Id: SanityCheck.h,v 1.7 2009/01/09 14:40:59 lukall Exp $
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

  void readWeights(istream & weightStream, vector<double>& w);
  int getInitDirection(vector<Scores>& testset,vector<Scores> &trainset, Normalizer * pNorm,vector<vector<double> >& w,double test_fdr);
  virtual bool validateDirection(vector<vector<double> >& w);
  void resetDirection(vector<vector<double> >& w);

  static void setInitWeightFN(string fn) {initWeightFN=fn;}
  static void setInitDefaultDir(int dir) {initDefaultDir=dir;}
  static void setOverrule(bool orl) {overRule=orl;}
protected:
  virtual void getDefaultDirection(vector<vector<double> >& w);
  int initPositives;
  double fdr;
  static bool overRule;
  static string initWeightFN;
  static int initDefaultDir; // Default Direction, 0=do not use,
                             // positive integer = feature number,
                             // negative integer = lower score better
  vector<Scores> *pTestset, *pTrainset;
};

#endif /*SANITYCHECK_H_*/
