/*******************************************************************************
 Copyright 2006-2012 Lukas Käll <lukas.kall@scilifelab.se>

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

 *****************************************************************************/
#include <algorithm>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <cmath>
#include <float.h>
#include "Globals.h"
#include "LTSRegression.h"

// initialize static variables
int LTSRegression::noSubsets = 500;
// maximum difference between q2 and q1 to achieve convergence
double LTSRegression::epsilon = 0.0001;
// cardinal of |H|(percentage of the total number of points used to build the regression line)
double LTSRegression::percentageH = 0.95;

LTSRegression::LTSRegression() :
  h(-1), regCoefficients(make_pair(0.0, 0.0)) {
}

LTSRegression::~LTSRegression() {
}

// function to compare 2 data points according to r
bool compareDataPoints(dataPoint x, dataPoint y) {
  return (x.absr < y.absr);
}

// set the data points
void LTSRegression::setData(vector<double> & x, vector<double> & y) {
  dataPoint tmp;
  for (int i = 0; i < x.size(); ++i) {
    tmp.x = x[i];
    tmp.y = y[i];
    tmp.absr = -1.0;
    data.push_back(tmp);
  }
  // since we definitely expect less than 25% contamination, we set h as 0.75*n
  h = (int)round(percentageH * x.size());
}

// for our case, constructing a random p-subset is equivalent to build the equation of a line through 2 randomly
// chosen points
vector<dataPoint> LTSRegression::getInitialHSubset() {
  double a, b;
  int i1, i2;
  int n = data.size();
  // generate two random indices
  i1 = PseudoRandom::lcg_rand() % n;
  i2 = PseudoRandom::lcg_rand() % n;
  while (i2 == i1) {
    i2 = PseudoRandom::lcg_rand() % n;
  }
  // calculate the a and b of the equation of the line going through the two points selected above (y = ax + b)
  a = (data[i2].y - data[i1].y) / (data[i2].x - data[i1].x);
  b = data[i1].y - (data[i1].x * a);
  // calculate the residuals and sort the data according to these values
  fillResiduals(make_pair(a, b));
  partial_sort(data.begin(),
               data.begin() + h,
               data.end(),
               compareDataPoints);
  return vector<dataPoint> (data.begin(), data.begin() + h);
}

// fill the absolute values of the residuals
void LTSRegression::fillResiduals(pair<double, double> par) {
  for (int i = 0; i < data.size(); ++i) {
    data[i].absr = abs(data[i].y - (par.first * data[i].x + par.second));
  }
}

// fit a line to the points in h using the least-squares method
pair<double, double> LTSRegression::fitLSLine(vector<dataPoint> hdata) {
  double sumxy = 0.0, sumx = 0.0, sumy = 0.0, sumxsq = 0.0;
  double a, b;
  for (int i = 0; i < h; ++i) {
    sumx += hdata[i].x;
    sumy += hdata[i].y;
    sumxy += hdata[i].x * hdata[i].y;
    sumxsq += pow(hdata[i].x, 2);
  }
  a = ((h * sumxy) - (sumx * sumy)) / ((h * sumxsq) - (pow(sumx, 2)));
  b = (sumy / (double)h) - (a * (sumx / (double)h));
  //cout << "a , b = " << a << ", " << b <<  endl;
  return make_pair(a, b);
}

// perform a C step (fit LS line, calculate residuals, return the h points with the lowest abs(residual))
vector<dataPoint> LTSRegression::performCstep(vector<dataPoint> hOld) {
  pair<double, double> par;
  // fit a least-squares regression line using data points from hold
  par = fitLSLine(hOld);
  // calculate and fill the residuals for all the data points
  fillResiduals(par);
  // sort the data points according to absolute values of residuals
  partial_sort(data.begin(),
               data.begin() + h,
               data.end(),
               compareDataPoints);
  return vector<dataPoint> (data.begin(), data.begin() + h);
}

// calculate the sum of squared residuals
double LTSRegression::calculateQ() {
  double res = 0.0;
  for (int i = 0; i < h; ++i) {
    res += pow(data[i].absr, 2);
  }
  return res;
}

void LTSRegression::runLTS() {
  double q1, q2, bestq = DBL_MAX, largestQ = 0.0;
  pair<double, double> par;
  vector<dataPoint> hold, hnew, besth(data.begin(), data.begin() + h);
  vector<vector<dataPoint> > best10Subsets;
  int noBestSubsets = 0, indexLargestQ = 0;
  vector<double> Q;
  if (VERB > 3) {
    cerr << "Regression parameters: " << endl;
    cerr << "   h = " << h << " = " << percentageH * 100
        << "%, no_initial_subsets = " << noSubsets << ", epsilon = "
        << epsilon << endl;
  }
  srand(time(NULL));
  for (int i = 0; i < noSubsets; ++i) {
    // get the first subset (data will be sorted according to abs(residuals)
    hold = getInitialHSubset();
    // apply 2 C-steps
    hnew = performCstep(hold);
    hold = performCstep(hnew);
    // compute Q
    par = fitLSLine(hold);
    fillResiduals(par);
    q1 = calculateQ();
    // include the current subset in the list of best 10 subsets if necessary
    if (noBestSubsets < 10) {
      best10Subsets.push_back(hold);
      Q.push_back(q1);
      if (q1 > largestQ) {
        largestQ = q1;
        indexLargestQ = noBestSubsets;
      }
      noBestSubsets++;
    } else {
      if (q1 < largestQ) {
        best10Subsets[indexLargestQ] = hold;
        Q[indexLargestQ] = q1;
        largestQ = Q[0];
        indexLargestQ = 0;
        for (int i = 1; i < Q.size(); ++i)
          if (Q[i] > largestQ) {
            largestQ = Q[i];
            indexLargestQ = i;
          }
      }
    }
  }
  // for the best 10 h-subset perform C-steps until convergence
  for (int j = 0; j < best10Subsets.size(); ++j) {
    hold = best10Subsets[j];
    par = fitLSLine(hold);
    fillResiduals(par);
    q2 = calculateQ();
    do {
      q1 = q2;
      hnew = performCstep(hold);
      par = fitLSLine(hnew);
      fillResiduals(par);
      q2 = calculateQ();
      hold = hnew;
    } while (abs(q2 - q1) > epsilon);
    if (q2 < bestq) {
      regCoefficients = par;
      bestq = q2;
      besth = hnew;
    }
  }
  /*if (VERB > 3)
   {
   cerr << "Data points used to generate the regression line: " << endl;
   printVector(besth);
   cerr << "Q = " << bestq << endl;
   }*/
  if (VERB > 2) {
    cerr << "Final LTS equation: y = " << regCoefficients.first
        << " * x + " << regCoefficients.second << endl;
  }
}

void LTSRegression::printVector(vector<dataPoint> v) {
  for (int i = 0; i < v.size(); ++i) {
    cerr << v[i].x << " " << v[i].y << endl;
  }
  cerr << endl;
}

void LTSRegression::printDataPoints() {
  for (int i = 0; i < h; ++i) {
    cout << data[i].x << " " << data[i].y << " " << data[i].absr << endl;
  }
  cout << "----" << endl;
  for (int i = h; i < data.size(); ++i) {
    cout << data[i].x << " " << data[i].y << " " << data[i].absr << endl;
  }
}
