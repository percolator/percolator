/*******************************************************************************
 Copyright (c) 2008 Lukas Käll

 Permission is hereby granted, free of charge, to any person
 obtaining a copy of this software and associated documentation
 files (the "Software"), to deal in the Software without
 restriction, including without limitation the rights to use,
 copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the
 Software is furnished to do so, subject to the following
 conditions:

 The above copyright notice and this permission notice shall be
 included in all copies or substantial portions of the Software. 
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 OTHER DEALINGS IN THE SOFTWARE.
 
 $Id: PosteriorEstimator.h,v 1.6 2008/05/20 00:24:43 lukall Exp $
 
 *******************************************************************************/

#ifndef POSTERIORESTIMATOR_H_
#define POSTERIORESTIMATOR_H_
#include<vector>
#include<utility>
using namespace std;
#include "LogisticRegression.h"

class PosteriorEstimator
{
public:
  PosteriorEstimator();
  virtual ~PosteriorEstimator();
  bool parseOptions(int argc, char **argv);   
  void run();
  static void estimatePEP( vector<pair<double,bool> >& combined, double pi0, vector<double>& peps);
  static void estimate( vector<pair<double,bool> >& combined, LogisticRegression& lr, double pi0);
  static void getPValues(const vector<pair<double,bool> >& combined, vector<double>& p);
  static void getQValues(double pi0,const vector<pair<double,bool> >& combined, vector<double>& q);
  static double estimatePi0(vector<pair<double,bool> >& combined, const unsigned int numBoot=100);
  static void setReversed(bool status) {reversed = status;}
protected:
  void finishStandalone(vector<pair<double,bool> >& combined, const vector<double>& peps, double pi0);
  static void binData(const vector<pair<double,bool> >& combined, vector<double>& medians, 
               vector<unsigned int>& negatives, vector<unsigned int>& sizes);
  string targetFile,decoyFile;
  static bool reversed;
};

#endif /*POSTERIORESTIMATOR_H_*/
