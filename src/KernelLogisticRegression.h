/*******************************************************************************

 *******************************************************************************/

#ifndef KERNELLOGISTICREGRESSION_H_
#define KERNELLOGISTICREGRESSION_H_
#include<vector>
using namespace std;
#include <assert.h>
#include "Transform.h"
#include "Numerical.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
namespace ublas = boost::numeric::ublas;

class KernelLogisticRegression{
  public:
    KernelLogisticRegression(){};
    virtual ~KernelLogisticRegression(){};
    double PEPeval(double xx);
    static double bandwidth;
    static double Epsilon;

    void predict(const vector<double>& x, vector<double>& predict);
    double predict(double xx) {
	return PEPeval(double xx);
    }
    void IRLSKernelLogisticRegression();
    void setDatafromBase(const vector<double>& xx)
    void setData(const vector<double>& xx, const vector<unsigned int>& yy,
                 const vector<unsigned int>& mm) {
      setDatafromBase(xx);
      y = yy;
      m = mm;
    }


  protected:
    virtual void calcPZW();
    Transform transf;
    vector<unsigned int> y, m;
    static const double gRange;
    double beta;
    ublas::matrix<double> w,;
    ublas::vector<double> alphas,z,g;
};

#endif /*KERNELLOGISTICREGRESSION_H_*/
