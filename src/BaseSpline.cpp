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

 *******************************************************************************/

#include<iterator>
#include<vector>
#include<algorithm>
#include<numeric>
#include<functional>
#include<cmath>
#include "BaseSpline.h"
#include "Globals.h"

using namespace std;

class SplinePredictor {
    BaseSpline* bs;
  public:
    SplinePredictor(BaseSpline* b) {
      bs = b;
    }
    double operator()(double x) {
      return bs->predict(x);
    }
};

double BaseSpline::convergeEpsilon = 1e-4;
double BaseSpline::stepEpsilon = 1e-8;
double BaseSpline::weightSlope = 1e1;
double BaseSpline::scaleAlpha = 1;

double BaseSpline::splineEval(double xx) {
  xx = transf(xx);
  std::size_t n = x.size();
  vector<double>::iterator left, right = lower_bound(x.begin(),
                                                     x.end(),
                                                     xx);
  if (right == x.end()) {
    double derl = (g[n - 1] - g[n - 2]) / (x[n - 1] - x[n - 2])
        + (x[n - 1] - x[n - 2]) / 6 * gamma[n - 3];
    double gx = g[n - 1] + (xx - x[n - 1]) * derl;
    return gx;
  }
  int rix = static_cast<int>(right - x.begin());
  if (*right == xx) {
    return g[rix];
  }
  if (rix > 0) {
    left = right;
    left--;
    double dr = *right - xx;
    double dl = xx - *left;
    double gamr = (rix < static_cast<int>(n - 1) ? gamma[rix - 1] : 0.0);
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

void BaseSpline::roughnessPenaltyIRLS() {
  initiateQR();
  initg();
  double p1 = 1 - tao;
  double p2 = tao;
  double alpha = alphaLinearSearchBA(0.0,
                          1.0,
                          p1,
                          p2,
                          evaluateSlope(-scaleAlpha*log(p1)),
                          evaluateSlope(-scaleAlpha*log(p2)));
  if (VERB > 2) {
    cerr << "Alpha selected to be " << alpha << endl;
  }
  iterativeReweightedLeastSquares(alpha);
}

void BaseSpline::iterativeReweightedLeastSquares(double alpha) {
  double step = 0.0;
  int iter = 0;
  unsigned int n = static_cast<unsigned int>(x.size());
  do {
    g = gnew;
    calcPZW();
    PackedMatrix diag = PackedMatrix::packedDiagonalMatrix(PackedVector(static_cast<int>(n), 1) / w).packedMultiply(alpha);
    PackedMatrix aWiQ = (diag).packedMultiply(Q);
    PackedMatrix M = R.packedAdd(Qt.packedMultiply(aWiQ));
    gamma = Qt.packedMultiply(z);
    solveInPlace(M,gamma);
    gnew = z.packedSubtract(aWiQ.packedMultiply(gamma));
    limitg();
    PackedVector difference = g.packedSubtract(gnew);
    step = packedNorm(difference) / n;
    if (VERB > 3) {
      cerr << "step size:" << step << endl;
    }
  } while ((step > stepEpsilon || step < 0.0) && (++iter < 20));
  g = gnew;
}

double BaseSpline::alphaLinearSearchBA(double min_p,
                                       double max_p,
                                       double p1, double p2,
                                       double cv1, double cv2) {
  // Minimize Slope score
  // Use neg log of 0<p<1 so that we allow for searches 0<alpha<inf
  double oldCV = 0.0;
  if (cv2 < cv1) {
    // keep point 2
    min_p = p1;
    p1 = p2;
    p2 = min_p + tao * (max_p - min_p);
    oldCV = cv1;
    cv1 = cv2;
    cv2 = evaluateSlope(-scaleAlpha*log(p2));
    if (VERB > 3) {
      cerr << "New point with alpha=" << -scaleAlpha*log(p2) << ", giving slopeScore=" << cv2 << endl;
    }
  } else {
    // keep point 1
    max_p = p2;
    p2 = p1;
    p1 = min_p + (1 - tao) * (max_p - min_p);
    oldCV = cv2;
    cv2 = cv1;
    cv1 = evaluateSlope(-scaleAlpha*log(p1));
    if (VERB > 3) {
      cerr << "New point with alpha=" << -scaleAlpha*log(p1) << ", giving slopeScore=" << cv1 << endl;
    }
  }
  if ((oldCV - min(cv1, cv2)) / oldCV < 1e-5 || (abs(p2 - p1) < 1e-10)) {
    return (cv1 < cv2 ? -scaleAlpha*log(p1) : -scaleAlpha*log(p2));
  }
  return alphaLinearSearchBA(min_p, max_p, p1, p2, cv1, cv2);
}

void BaseSpline::initiateQR() {
  // NOTE: at least 4 data points are needed to construct a quadratic spline!
  assert(x.size() >= 4);
  int n = static_cast<int>(x.size());
  dx.resize(n-1);
  for (std::size_t ix = 0; static_cast<int>(ix) < n - 1; ix++) {
    dx.addElement(static_cast<int>(ix), x[ix + 1] - x[ix]);
    assert(dx[ix] > 0);
  }
  Q = PackedMatrix(n,n-2);
  R = PackedMatrix(n-2, n-2);
  //Fill Q
  Q[0].packedAddElement(0, 1 / dx[0]);
  Q[1].packedAddElement(0, -1 / dx[0] - 1 / dx[1]);
  Q[1].packedAddElement(1, 1 / dx[1]);
  for (int j = 2; j < n - 2; j++) {
    Q[j].packedAddElement(j-2, 1 / dx[j - 1]);
    Q[j].packedAddElement(j-1, -1 / dx[j - 1] - 1 / dx[j]);
    Q[j].packedAddElement(j, 1 / dx[j]);
  }
  Q[n - 2].packedAddElement(n-4, 1 / dx[n - 3]);
  Q[n - 2].packedAddElement(n - 3, -1 / dx[n - 3] - 1 / dx[n - 2]);
  Q[n - 1].packedAddElement(n - 3,  1 / dx[n - 2]);
  //Fill R
  for (int i = 0; i < n - 3; i++) {
    R[i].packedAddElement(i, (dx[i] + dx[i + 1]) / 3);
    R[i].packedAddElement(i + 1, dx[i + 1] / 6);
    R[i + 1].packedAddElement(i, dx[i + 1] / 6);
  }
  R[n - 3].packedAddElement(n - 3, (dx[n - 3] + dx[n - 2]) / 3);
  Qt = PackedMatrix(n-2,n);
  Qt = PackedMatrix::packedTranspose(Q);
}

double BaseSpline::evaluateSlope(double alpha) {
  // Calculate a spline for current alpha
  iterativeReweightedLeastSquares(alpha);
  // Find highest point (we only want to evaluate things to the right of that point)
  int n = g.numberEntries();
  int mixg=1; // Ignore 0 and n-1
  double maxg = g[mixg];
  for (int ix=mixg;ix<n-1;++ix) {
    assert(ix==g.index(ix)); //This should be a filled vector
    if (g[ix]>=maxg) {
      maxg = g[ix];
      mixg = ix;
    }
  }
  double maxSlope = -10e6;
  int slopeix=-1;
  for (int ix=mixg+1;ix<n-2; ++ix) {
    double slope=g[ix-1]-g[ix];
    if (slope>maxSlope) {
      maxSlope = slope;
      slopeix = ix;
    }
  }
  // Now score the fit based on a linear combination between
  // The bump area and alpha

  if (VERB > 3) {
    cerr << "mixg=" << mixg << ", maxg=" << maxg << ", maxBA=" << maxSlope << " at ix=" << slopeix << ", alpha=" << alpha << endl;
  }
  return maxSlope*weightSlope + alpha;
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
    if (VERB > 1) {
      cerr << "Logit transforming all scores prior to PEP calculation"
        << endl;
    }
    transf = Transform(minV > 0.0 ? 0.0 : 1e-20,
                       maxV < 1.0 ? 0.0 : 1e-10,
                       true);
  } else if (minV >= 0.0) {
    if (VERB > 1) {
      cerr << "Log transforming all scores prior to PEP calculation"
          << endl;
    }
    transf = Transform(minV > 0.0 ? 0.0 : 1e-20, 0.0, false, true);
  }
  transform(xx.begin(), xx.end(), back_inserter(x), transf);
}

void BaseSpline::solveInPlace(PackedMatrix& mat, PackedVector& res) {
  res = res.makeSparse();
  // Current implementation requires a quadratic mat
  int nCol = mat.numCols();
  PackedVector nonEmpty;
  int col, row, rowPos;
  for (col = 0; col < nCol; col++) {
    //int stop = 1;
    //if(col==stop-1){
//      cout << "*********************"<<endl;
//      cout << col << endl;
//      cout << "*********************"<<endl;
//      mat.displayMatrix();
//      cout << "RES: " <<res.values<<endl;
    //}

    // find the non-null elements in this column (below row "col")
    nonEmpty = PackedVector();
    int pivotPos(-1);
    for (row = col; row < mat.numRows(); row++) {
      int rowEntries = mat[row].numberEntries();
      for (rowPos = rowEntries; rowPos--;) {
        int entryIndex = mat[row].index(rowPos);
        if (entryIndex == col) {
          nonEmpty.packedAddElement(row, mat[row][rowPos]);
          if (row == col) {
            pivotPos = nonEmpty.numberEntries()-1;
          }
        }
      }
    }
//    if(col==stop-1){
//      cout<<"nonEmpty " << nonEmpty.values<< endl;
//      cout<<"pivot " << pivotPos<< endl<<endl;
//    }

    //find most significant row
    double maxVal(0.0);
    int maxRow(-1), maxRowPos(-1);
    for (rowPos = nonEmpty.numberEntries(); rowPos--;) {
      double val = nonEmpty[rowPos];
      assert(isfinite(val));
      if (fabs(val) > fabs(maxVal)) {
        maxVal = val;
        maxRow = nonEmpty.index(rowPos);
        maxRowPos = rowPos;
      }
    }
    //int imaxRow = maxRow;
    //int imaxRowPos = maxRowPos;
//    if(col==stop-1){
//      cout << "imaxRow " << imaxRow << endl;
//      cout << "imaxRowPos " << imaxRowPos << endl<<endl;
//    }

    // Put the most significant row at row "col"
    if (maxVal != 0.0) {
      if (maxRow != col) {
        swap(mat[col], mat[maxRow]);
        res.swapElements(col, maxRow);

        if (pivotPos >= 0) {
          nonEmpty.swapElements(maxRowPos, pivotPos);
        }
      }
//      if(col==stop-1){
//        cout << "SWAP\n";
//        mat.displayMatrix();
//      }

      // Divide the row with maxVal
      mat[col].packedDiv(maxVal);
      double value = res[col] / maxVal;
      res.packedReplace(col,value);
//      if(col==stop-1){
//        cout << "\nDIVIDE ROW\n";
//        mat.displayMatrix();
//        cout << "res " << res.values<<endl;
//        cout << "nonEmpty " << nonEmpty.values<<endl;
//        cout << "\n\n";
//      }
    }
    // subtract the row from other rows
    for (rowPos = nonEmpty.numberEntries(); rowPos--;) {
      row = nonEmpty.index(rowPos);
      if (row == col) {
        continue;
      }
      // If the pivotRow was empty (prior to swap) at col=row do not process this row
      if (pivotPos < 0 && row == maxRow) {
        continue;
      }
      double val = nonEmpty[rowPos];
      PackedVector prodVector = mat[col];
      mat[row] = mat[row].packedSubtract(prodVector.packedProd(val));
      double value = res[row] - (val * res[col]);
      res.packedReplace(row, value);

//      if(col==stop-1){
//        cout << "SUBTRACT ROW\n";
//        mat.displayMatrix();
//        cout << "res " << res.values<<endl;
//        cout << "nonEmpty " << nonEmpty.values<<endl;
//        cout << "\n\n";
//      }
    }
  }
  // Go bottom up and clear upper halfmatrix
//  mat.displayMatrix();
//  res.displayVector();
//  nonEmpty.displayVector();
//  cout << "\n\n\n\n";
  for (col = mat.numCols(); col--;) {
    nonEmpty = PackedVector();
    for (row = 0; row < col; row++) {
      for (rowPos = mat[row].numberEntries(); rowPos--;) {
        if (mat[row].index(rowPos) == col) {
          nonEmpty.packedAddElement(row, mat[row][rowPos]);
        }
      }
    }
//    cout << nonEmpty.values<<endl;
    // subtract the row from other rows
    for (rowPos = nonEmpty.numberEntries(); rowPos--;) {
      row = nonEmpty.index(rowPos);
      double val = nonEmpty[rowPos];
      mat[row] = mat[row].packedSubtract(mat[col].packedProd(val));
      double value = res[row] - (val * res[col]);
      res.packedReplace(row, value);
    }
//    cout << res.values<<endl;
  }
  return;
}


