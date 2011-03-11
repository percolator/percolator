/*
 * ProteinProbEstimatorHelper.h
 *
 *  Created on: Feb 25, 2011
 *      Author: tomasoni
 */

/*******************************************************************************
 Copyright 2006-2009 Lukas Käll <lukas.kall@cbr.su.se>

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

 *******************************************************************************/

#ifndef PROTEINPROBESTIMATORHELPER_H_
#define PROTEINPROBESTIMATORHELPER_H_

#include <math.h>
#include <string>
#include "Vector.h"
#include <set>
using namespace std;
#include <ext/hash_set>

namespace __gnu_cxx {
template<> struct hash< std::string > {
    size_t operator()( const std::string & x ) const {
      return hash< const char* >()( x.c_str() );
    }
};
}

int matchCount(const __gnu_cxx::hash_set<string> & positiveNames,
    const Array<string> & atThreshold) {
  int count = 0;
  for (int k=0; k<atThreshold.size(); k++) {
    if ( positiveNames.count( atThreshold[k] ) > 0 )
      count++;
  }
  return count;
}

Array<string> matches(const __gnu_cxx::hash_set<string> & positiveNames,
    const Array<string> & atThreshold) {
  Array<string> result;
  for (int k=0; k<atThreshold.size(); k++) {
    if ( positiveNames.count( atThreshold[k] ) > 0 )
      result.add( atThreshold[k] );
  }
  return result;
}

/**
 * calculates empirical and estimated FDRs and stores the results in the
 * estimatedFdr and empiricalFdr Arrays for use in calculateMSE_FDR()
 */
void calculateFDRs(
    const fidoOutput output,
    const Array<string>& truePositives, const Array<string>& falsePositives,
    Array<double>& estimatedFdrs, Array<double>& empiricalFdrs) {
  __gnu_cxx::hash_set<string> truePosSet(truePositives.size()),
      falsePosSet(falsePositives.size());
  int k;
  for (k=0; k<truePositives.size(); k++)
    truePosSet.insert(truePositives[k]);
  for (k=0; k<falsePositives.size(); k++)
    falsePosSet.insert(falsePositives[k]);
  Array<string> protsAtThreshold;
  string line;
  double prob, lastProb=-1;
  int fpCount = 0, tpCount = 0;
  int numScored = 0;
  Array<string> observedProteins;
  double estFDR = 0.0;
  double empiricalFDR = 0.0;
  double totalFDR = 0.0;
  bool scheduledUpdate = false;

  for(k=0; k<output.peps.size(); k++){
    prob = output.peps[k];
    protsAtThreshold = output.protein_ids[k];
    numScored += protsAtThreshold.size();
    observedProteins.append(protsAtThreshold);
    int fpChange = matchCount(falsePosSet, protsAtThreshold);
    int tpChange = matchCount(truePosSet, protsAtThreshold);

    if ( prob != lastProb && lastProb != -1 ){
      scheduledUpdate = true;
    }
    if ( scheduledUpdate ) {
      if ( fpChange > 0 || tpChange > 0) {
        estimatedFdrs.add(estFDR);
        empiricalFdrs.add(empiricalFDR);
        scheduledUpdate = false;
      }
    }

    fpCount += fpChange;
    tpCount += tpChange;
    totalFDR += (1-prob) * (fpChange + tpChange);
    estFDR = totalFDR / (fpCount + tpCount);
    empiricalFDR = double(fpCount) / (fpCount + tpCount);
    lastProb = prob;
  }
  lastProb = prob;
  {
    estimatedFdrs.add(estFDR);
    empiricalFdrs.add(empiricalFDR);
  }
  // uncomment the following lines to output the results to cout and file
  //cout.precision(10);
  //cout << estimatedList << " " << empericalList << endl;
  //ofstream fout("/tmp/fido/rlistFDROut.txt");
  //fout << estimatedList << endl << empericalList << endl;
  //fout.close();
}

double squareAntiderivativeAt(double m, double b, double xVal) {
  double u = m*m;
  double v = 2*m*b;
  double t = b*b;
  return u*xVal*xVal*xVal/3.0 + v*xVal*xVal/2.0 + t*xVal;
}

double antiderivativeAt(double m, double b, double xVal) {
  return m*xVal*xVal/2.0 + b*xVal;
}

double area(double x1, double y1, double x2, double y2, double threshold) {
  double m = (y2-y1)/(x2-x1);
  double b = y1-m*x1;
  double area = squareAntiderivativeAt(m, b, min(threshold, x2) )
      - squareAntiderivativeAt(m, b, x1);
  return area;
}


/**
 * calculates the FDR Mean Square Error
 *
 * @return FDR Mean Square Error
 */
double calculateMSE_FDR(double threshold,
    const Array<double>& estimatedFdr, const Array<double>& empiricalFdr) {
  assert(estimatedFdr.size() == empiricalFdr.size());
  Vector diff = Vector(estimatedFdr) - Vector(estimatedFdr);
  double tot = 0.0;
  int k;
  for (k=0; k<diff.size()-1; k++) {
    // stop if no part of the estFDR is < threshold
    if (estimatedFdr[k] >= threshold) {
      if (k == 0)
        tot = 1.0 / 0.0;
      break;
    }
    tot += area(estimatedFdr[k],diff[k],estimatedFdr[k+1],diff[k+1],threshold);
  }
  double xRange = min(threshold, estimatedFdr[k]) - estimatedFdr[0];

  if (isinf(tot))
    return tot;
  else return (tot/xRange);
}

/**
 * calculates the roc curve and stores the results in the fps and tps Arrays
 * for use in calculateROC50
 *
 */
void calculateRoc(const fidoOutput output,
    const Array<string>& truePositives, const Array<string>& falsePositives,
    Array<int>& fps, Array<int>& tps) {

  __gnu_cxx::hash_set<string> truePosSet(truePositives.size()),
      falsePosSet(falsePositives.size());
  int k;
  for (k=0; k<truePositives.size(); k++) {
    truePosSet.insert( truePositives[k] );
  }
  for (k=0; k<falsePositives.size(); k++) {
    falsePosSet.insert( falsePositives[k] );
  }
  Array<string> protsAtThreshold;
  string line;
  double prob, lastProb=-1;
  int fpCount = 0, tpCount = 0;
  int numScored = 0;
  Array<string> observedProteins;
  fps.add(0);
  tps.add(0);
  bool scheduledUpdate = false;
  double totalFDR = 0.0, estFDR = 0.0;

  for(k=0; k<output.peps.size(); k++) {
    prob = output.peps[k];
    protsAtThreshold = output.protein_ids[k];
    numScored += protsAtThreshold.size();
    observedProteins.append( protsAtThreshold );
    int fpChange = matchCount(falsePosSet, protsAtThreshold);
    int tpChange = matchCount(truePosSet, protsAtThreshold);
    if ( prob != lastProb && lastProb != -1 ) {
      scheduledUpdate = true;
    }
    if ( scheduledUpdate ) {
      fps.add( fpCount );
      tps.add( tpCount );
      scheduledUpdate = false;
      totalFDR += (1-prob) * (fpChange + tpChange);
      estFDR = totalFDR / (fpCount + tpCount);
    }
    fpCount += fpChange;
    tpCount += tpChange;
    lastProb = prob;
  }
  lastProb = prob;
  fps.add( fpCount );
  tps.add( tpCount );
  fps.add( falsePosSet.size() );
  tps.add( truePosSet.size() );

  // uncomment the following lines to output the results to cout and file
  //cout.precision(10);
  //cout << fps << " " << tps << endl;
  //ofstream fout("/tmp/fido/rlistROCOut.txt");
  //fout << fps << endl << tps << endl;
  //fout.close();
}

/**
 * calculates the area under the roc curve up to 50 false positives
 *
 * @return roc50
 */
double calculateROC50(int N, const Array<int>& fps, const Array<int>& tps){
  double rocN = 0.0;
  if ( fps.back() < N ) {
    cerr << "There are not enough false positives; needed " << N
        << " and was only given " << fps.back() << endl << endl;
    exit(1);
  }
  for (int k=0; k<fps.size()-1; k++) {
    // find segments where the fp value changes
    if ( fps[k] >= N )
      break;
    if ( fps[k] != fps[k+1] ) {
      // this line segment is a function
      double currentArea = area(fps[k], tps[k], fps[k+1], tps[k+1], N);
      rocN += currentArea;
    }
  }
  double roc50 = rocN / (N * tps.back());
  return roc50;
}

#endif /* PROTEINPROBESTIMATORHELPER_H_ */
