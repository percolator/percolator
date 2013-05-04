//May 2013, Xiao.
//
#include <iostream>
#include<vector>
#include<algorithm>
#include<numeric>

#include<assert.h>
#include "Numerical.h"
#define PI 3.14159265
#include "Globals.h"

#include <boost/numeric/ublas/blas.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include "KernelLogisticRegression.h"

extern "C"{
#include <atlas/clapack.h>
#include <atlas/cblas.h>
}

namespace ublas = boost::numeric::ublas;
namespace ublas1 = boost::numeric::ublas::blas_1;
namespace ublas2 = boost::numeric::ublas::blas_2;
namespace ublas3 = boost::numeric::ublas::blas_3;
using namespace std;

const double KernelLogisticRegression::gRange = 35.0;
double KernelLogisticRegression::bandwidth = 0.9;
double KernelLogisticRegression::Epsilon = 1e-8;



class Predictor {
	KernelLogisticRegression* bs;
	public:
	Predictor(KernelLogisticRegression* b){
		bs = b;}
	double operator() (double x ) {
		return bs -> predict(x);
	}
};

void KernelLogisticRegression::setDatafromBase(const vector<double>& xx) {
// copy from Basespline.cpp, reset x.
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

double KernelLogisticRegression::kernel(double x1, double x2,double h) {
        double ker;
        double u = (x1-x2)/h;
        ker = 1/sqrt(2* PI )*exp(-u*u/2);//Guassian kernel. could be fixed.
        return ker;
}

void KernelLogisticRegression::IRLSKernelLogisticRegression()
{
        double step = 0.0,beta_old=0.0,gamma=1.0,h=bandwidth;
        int iter = 0,nrhs=1;
	size_t N=x.size();
        ublas::vector<double> ones(N),xi(N),zeta(N);
        ublas::matrix<double> M(N,N),K(N,N),identity(N,N);

        for (int i = 0; i<N; i++){
                identity(i,i)=1;
                ones(i)=1;
		alphas(i)=0;
                for (int j=0;j<N;j++){
                        K(i , j )= kernel(x[i],x[j],h);
                }
        }
        do {
                g =  ublas2::gmv(g,0,gamma,K,alphas);
                calcPZW();
                M=K;
                ublas3::gmm(M,gamma,gamma,identity,w);//M = I*w+gamma*K
                ublas::matrix<double> M1 = M;
                xi = ones;
                clapack_dposv(CblasRowMajor,CblasUpper,N,nrhs,&M1(0,0), N,&xi(0),N);
                zeta = z;
                clapack_dposv(CblasRowMajor,CblasUpper,N,nrhs,&M(0,0), N,&zeta(0),N);
                beta = ublas1::dot(ones,zeta)/ublas1::dot(ones,xi);
                ublas2::gmv(xi,-beta,gamma,identity,zeta);//xi=I*zeta-beta*xi
                alphas = xi;
                step = (beta_old-beta)*(beta_old-beta);
                beta_old = beta;
        }while((step > Epsilon || step <0.0) && (++iter < 20) );
}


void KernelLogisticRegression::calcPZW(){//same as the function in LogisticRegression.cpp
        size_t N= y.size();
        for (int ix = 0; ix < N; ix++){
                g(ix) = g(ix)+beta;
                double e = exp(g(ix));
                assert(isfinite(e));
                Numerical num(1e-15);
                double p = min(max(e / (1 + e), num.epsilon), 1- num.epsilon);
                assert(isfinite(p));
                double mu = m[ix] * p;
                double sigma = mu * (1-p);
                w(ix,ix)= max(1/sigma, num.epsilon);
                assert(isfinite(w(ix,ix)));
                z(ix)= min(gRange, max(-gRange, g(ix) + (((double)y[ix]) -mu ) / sigma));
                assert(isfinite(z(ix)));
                }
}


//calculate all peps of vector xx and stored in vector predict.
void KernelLogisticRegression::predict(const vector<double> &xx, 
					vector<double> & predict){
	predict.clear();
	transform( xx.begin(),xx.end(),back_inserter(predict),Predictor(this));
} 

double KernelLogisticRegression::PEPeval(double xx){
        size_t N=x.size();
        double gx=beta,ker1=0;
        for(int i =0; i<N; i++){
                ker1= kernel(xx, x[i],bandwidth);
                gx = gx+alphas[i]*ker1;
        }
        double pep = exp(gx);
        return pep;
}
