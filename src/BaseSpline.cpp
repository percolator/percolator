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

#include<iterator>
#include<vector>
#include<algorithm>
#include<numeric>
#include<functional>
#include<cmath>

using namespace std;

#include "ArrayLibrary.h"
#include "BaseSpline.h"
#include "Globals.h"

class SplinePredictor {
    BaseSpline *bs;
  public:
    SplinePredictor(BaseSpline *b) {
      bs = b;
    }
    double operator()(double x) {
      return bs->predict(x);
    }
};

double BaseSpline::convergeEpsilon = 1e-4;
double BaseSpline::stepEpsilon = 1e-8;

double BaseSpline::splineEval(double xx) {
  xx = transf(xx);
  size_t n = x.size();
  vector<double>::iterator left, right = lower_bound(x.begin(),
                                                     x.end(),
                                                     xx);
  if (right == x.end()) {
    double derl = (g[n - 1] - g[n - 2]) / (x[n - 1] - x[n - 2])
        + (x[n - 1] - x[n - 2]) / 6 * gamma[n - 3];
    double gx = g[n - 1] + (xx - x[n - 1]) * derl;
    return gx;
  }
  size_t rix = right - x.begin();
  if (*right == xx) return g[rix];
  if (rix > 0) {
    left = right;
    left--;
    double dr = *right - xx;
    double dl = xx - *left;
    double gamr = (rix < (n - 1) ? gamma[rix - 1] : 0.0);
    double gaml = (rix > 1 ? gamma[rix - 1 - 1] : 0.0);
    double h = *right - *left;
    double gx = (dl * g[rix] + dr * g[rix - 1]) / h - dl * dr / 6 * ((1.0
        + dl / h) * gamr + (1.0 + dr / h) * gaml);
    return gx;
  }
  // if (rix==0)
  double derr = (g[1] - g[0]) / (x[1] - x[0]) - (x[1] - x[0]) / 6
      * gamma[0];
  double gx = g[0] - (x[0] - xx) * derr;
  return gx;
}

static double tao = 2 / (1 + sqrt(5.0)); // inverse of golden section

void BaseSpline::iterativeReweightedLeastSquares() {

  Numerical::epsilon = 1e-15;
  unsigned int n = x.size(), alphaIter = 0;
  initiateQR();
  double alpha = .05, step, cv = 1e100;
  initg();
  do {
    int iter = 0;
    do {
      g = gnew;
      calcPZW();
      PackedMatrix aWiQ = alpha * diagonalPacked(Vec(n, 1) / w) * Q;
      PackedMatrix M = R + Qt * aWiQ;
      gamma = Qt * z;
      //      limitgamma();
      solveEquation<double> (M, gamma);
      gnew = z - aWiQ * gamma;
      limitg();
      step = norm(g - gnew) / n;
      if (VERB > 2) cerr << "step size:" << step << endl;
    } while ((step > stepEpsilon || step <= 0.0) && (++iter < 20));
    double p1 = 1 - tao;
    double p2 = tao;
    pair<double, double> res =
        alphaLinearSearch(0.0,
                          1.0,
                          p1,
                          p2,
                          crossValidation(-log(p1)),
                          crossValidation(-log(p2)));

    if (VERB > 2) cerr << "Alpha=" << res.first << ", cv=" << res.second
        << endl;
    assert(isfinite(res.second));

    if ((cv - res.second) / cv < convergeEpsilon || alphaIter++ > 100) {
      // Reject our last attempt to set alpha,
      // Return with the alpha used when setting g
      break;
    }
    cv = res.second;
    alpha = res.first;
  } while (true);
  if (VERB > 2) cerr << "Alpha selected to be " << alpha << endl;
  g = gnew;
}

pair<double, double> BaseSpline::alphaLinearSearch(double min_p,
                                                   double max_p,
                                                   double p1, double p2,
                                                   double cv1, double cv2) {
  double oldCV;
  if (cv1 > cv2) {
    min_p = p1;
    p1 = p2;
    oldCV = cv1;
    cv1 = cv2;
    p2 = min_p + tao * (max_p - min_p);
    cv2 = crossValidation(-log(p2));
    if (VERB > 3) cerr << "New point with alpha=" << -log(p2)
        << ", giving cv=" << cv2 << " taken in consideration" << endl;
  } else {
    max_p = p2;
    p2 = p1;
    oldCV = cv2;
    cv2 = cv1;
    p1 = min_p + (1 - tao) * (max_p - min_p);
    cv1 = crossValidation(-log(p1));
    if (VERB > 3) cerr << "New point with alpha=" << -log(p1)
        << ", giving cv=" << cv1 << " taken in consideration" << endl;
  }
  if ((oldCV - min(cv1, cv2)) / oldCV < 1e-6 || (abs(p2 - p1) < 1e-10)) return (cv1
      > cv2 ? make_pair(-log(p2), cv2) : make_pair(-log(p1), cv1));
  return alphaLinearSearch(min_p, max_p, p1, p2, cv1, cv2);
}

void BaseSpline::initiateQR() {
  int n = x.size();
  dx.resize(n - 1, 0.0);
  for (int ix = 0; ix < n - 1; ix++) {
    dx[ix] = x[ix + 1] - x[ix];
    assert(dx[ix] > 0);
  }
  R.resize(n - 2);
  Q.resize(n);

  //Fill Q
  Q[0].push_back(0, 1 / dx[0]);
  Q[1].push_back(0, -1 / dx[0] - 1 / dx[1]);
  Q[1].push_back(1, 1 / dx[1]);
  for (int j = 2; j < n - 2; j++) {
    Q[j].push_back(j - 2, 1 / dx[j - 1]);
    Q[j].push_back(j - 1, -1 / dx[j - 1] - 1 / dx[j]);
    Q[j].push_back(j, 1 / dx[j]);
  }
  Q[n - 2].push_back(n - 4, 1 / dx[n - 3]);
  Q[n - 2].push_back(n - 3, -1 / dx[n - 3] - 1 / dx[n - 2]);
  Q[n - 1].push_back(n - 3, 1 / dx[n - 2]);
  //Fill R
  for (int i = 0; i < n - 3; i++) {
    R[i].push_back(i, (dx[i] + dx[i + 1]) / 3);
    R[i].push_back(i + 1, dx[i + 1] / 6);
    R[i + 1].push_back(i, dx[i + 1] / 6);
  }
  R[n - 3].push_back(n - 3, (dx[n - 3] + dx[n - 2]) / 3);
  Qt = transpose(Q);
}

double BaseSpline::crossValidation(double alpha) {

  int n = R.size();
  //  Vec k0(n),k1(n),k2(n);
  vector<double> k0(n), k1(n), k2(n);

  PackedMatrix B = R + alpha * Qt * diagonalPacked(Vec(n + 2, 1.0) / w)
      * Q;
  // Get the diagonals from K
  // ka[i]=B[i,i+a]=B[i+a,i]
  for (int row = 0; row < n; ++row) {
    for (int rowPos = B[row].packedSize(); rowPos--;) {
      int col = B[row].index(rowPos);
      if (col == row) {
        k0[row] = B[row][rowPos];
      } else if (col + 1 == row) {
        k1[row] = B[row][rowPos];
      } else if (col + 2 == row) {
        k2[row] = B[row][rowPos];
      }
    }
  }

  // LDL decompose Page 26 Green Silverman
  // d[i]=D[i,i]
  // la[i]=L[i+a,i]
  //  Vec d(n),l1(n),l2(n);
  vector<double> d(n), l1(n), l2(n);
  d[0] = k0[0];
  l1[0] = k1[0] / d[0];
  d[1] = k0[1] - l1[0] * l1[0] * d[0];

  for (int row = 2; row < n; ++row) {
    l2[row - 2] = k2[row - 2] / d[row - 2];
    l1[row - 1] = (k1[row - 1] - l1[row - 2] * l2[row - 2] * d[row - 2])
        / d[row - 1];
    d[row] = k0[row] - l1[row - 1] * l1[row - 1] * d[row - 1]
        - l2[row - 2] * l2[row - 2] * d[row - 2];
  }
  // Find diagonals of inverse Page 34 Green Silverman
  // ba[i]=B^{-1}[i+a,i]=B^{-1}[i,i+a]
  //  Vec b0(n),b1(n),b2(n);
  vector<double> b0(n), b1(n), b2(n);
  for (int row = n; --row;) {
    if (row == n - 1) {
      b0[n - 1] = 1 / d[n - 1];
    } else if (row == n - 2) {
      b0[n - 2] = 1 / d[n - 2] - l1[n - 2] * b1[n - 2];
    } else {
      b0[row] = 1 / d[row] - l1[row] * b1[row] - l2[row] * b2[row];
    }
    if (row == n - 1) {
      b1[n - 2] = -l1[n - 2] * b0[n - 1];
    } else if (row >= 1) b1[row - 1] = -l1[row - 1] * b0[row] - l1[row]
        * b1[row];
    if (row >= 2) b2[row - 2] = -l1[row - 2] * b0[row];
  }

  // Calculate diagonal elements a[i]=Aii p35 Green Silverman
  // (expanding q according to p12)
  //  Vec a(n+2),c(n+1);
  vector<double> a(n), c(n - 1);
  for (int ix = 0; ix < n - 1; ix++)
    c[ix] = 1 / dx[ix];

  for (int ix = 0; ix < n; ix++) {
    if (ix > 0) {
      a[ix] += b0[ix - 1] * c[ix - 1] * c[ix - 1];
      if (ix < n - 1) {
        a[ix] += b0[ix] * (-c[ix - 1] - c[ix]) * (-c[ix - 1] - c[ix]);
        a[ix] += 2 * b1[ix] * c[ix] * (-c[ix - 1] - c[ix]);
        a[ix] += 2 * b1[ix - 1] * c[ix - 1] * (-c[ix - 1] - c[ix]);
        a[ix] += 2 * b2[ix - 1] * c[ix - 1] * c[ix];
      }
    }
    if (ix < n - 1) a[ix] += b0[ix + 1] * c[ix] * c[ix];
  }

  // Calculating weighted cross validation as described in p
  double cv = 0.0;
  for (int ix = 0; ix < n; ix++) {
    double f = (z[ix] - gnew[ix]) * w[ix] / (alpha * a[ix]);
    //    double f =(z[ix]-gnew[ix])/(alpha*alpha*a[ix]*a[ix]);
    cv += f * f * w[ix];
  }

  return cv;
}

void BaseSpline::predict(const vector<double>& xx, vector<double>& predict) {
  predict.clear();
  transform(xx.begin(),
            xx.end(),
            back_inserter(predict),
            SplinePredictor(this));
}

void BaseSpline::setData(const vector<double>& xx) {
  x.clear();
  double minV = *min_element(xx.begin(), xx.end());
  double maxV = *max_element(xx.begin(), xx.end());
  if (minV >= 0.0 && maxV <= 1.0) {
    if (VERB > 1) cerr
        << "Logit transforming all scores prior to PEP calculation"
        << endl;
    transf = Transform(minV > 0.0 ? 0.0 : 1e-20,
                       maxV < 1.0 ? 0.0 : 1e-10,
                       true);
  } else if (minV >= 0.0) {
    if (VERB > 1) cerr
        << "Log transforming all scores prior to PEP calculation" << endl;
    transf = Transform(minV > 0.0 ? 0.0 : 1e-20, 0.0, false, true);
  }
  transform(xx.begin(), xx.end(), back_inserter(x), transf);
}
