/*******************************************************************************

 *******************************************************************************/

#ifndef KERNELLOGISTICREGRESSION_H_
#define KERNELLOGISTICREGRESSION_H_
#include<vector>
using namespace std;
#include <assert.h>
#include "Transform.h"
#include "Numerical.h"

class KernelLogisticRegression{
  public:
    KernelLogisticRegression(){};
    virtual ~KernelLogisticRegression(){};
    void setData(const vector<double>& xx);
    void setData(const vector<double>& xx, const vector<unsigned int>& yy,
                 const vector<unsigned int>& mm) {
      setData(xx);
      y = yy;
      m = mm;
    }
    double kernel(double x1, double x2, double h);
    static double bandwidth;
    static double stepEpsilon;
    static double convergeEpsilon;
    void predict(const vector<double>& x, vector<double>& predict);
    double predict(double xx);
    void IRLSKernelLogisticRegression();
    
  protected:
    Transform transf;
    vector<double> x,alphas;
    vector<unsigned int> y, m;
    static const double gRange;
    double beta;
//    ublas::matrix<double> w(n,n);
//    ublas::vector<double> alphas(n),z(n),g(n);
    
};

#endif /*KERNELLOGISTICREGRESSION_H_*/
