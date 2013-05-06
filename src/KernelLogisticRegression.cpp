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

#include "KernelLogisticRegression.h"

extern "C"{
#include <atlas/clapack.h>
#include <atlas/cblas.h>
}

using namespace std;

const double KernelLogisticRegression::gRange = 35.0;
double KernelLogisticRegression::bandwidth = 0.9;
double KernelLogisticRegression::stepEpsilon = 1e-8;
double KernelLogisticRegression::convergeEpsilon=1e-4;


class Predictor { 
	KernelLogisticRegression* bs;
	public:
	Predictor(KernelLogisticRegression* b){
		bs = b;}
	double operator() (double x ) {
		return bs -> predict(x);
	}
};

void KernelLogisticRegression::setData(const vector<double>& xx) {
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
        ker = 1/sqrt(2* PI )*exp(-u*u/2);//Guassian kernel. could be changed.
        return ker;
}

void KernelLogisticRegression::IRLSKernelLogisticRegression()
{
        double step = 0.0,beta_old=0.0,gamma=1.0,h=bandwidth;
        int iter = 0,nrhs=1,incx=1,noncoef=0;
	size_t N=x.size();
        vector<double> ones(N),g(N),z(N);
        vector<double> K(N*N),identity(N*N),w(N*N);
	double p,mu,sigma;
	beta=0;
        for (int i = 0; i<N; i++){
                identity[i*N+i]=1;
                ones[i]=1;
		alphas.push_back(0);
                for (int j=0;j<N;j++){
                        K[i*N+j]= kernel(x[i],x[j],h);
                }
        }
        do {
		//blas subroutine:  http://www.netlib.org/blas/dgemv.f
                cblas_dgemv(CblasRowMajor,CblasNoTrans,N, N,1.0,&K[0],N,&alphas[0],incx,0.0,&g[0],incx);//g=K*alphas + 0*g
 		for (int ix = 0; ix < N; ix++){
                	g[ix] = g[ix]+beta;
			assert(isfinite(g[ix]));
                	p = 1 / (1 + exp(-g[ix]));
			assert(isfinite(p));
                	mu = m[ix] * p;
			if (mu == 0){sigma =1.0;}
			else{sigma = mu * (1-p);}
                	w[ix*N+ix]= 1/sigma;
                	z[ix]= g[ix] + (((double)y[ix]) -mu ) / sigma;
                	assert(isfinite(z[ix]));
                }
                vector<double> M=K;
		//blas subroutine: http://www.netlib.org/blas/dgemm.f
                cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, N,N,N,1,&identity[0],N,&w[0],N,gamma,&M[0],N);//M = I*w+gamma*K
                vector<double> M1 = M;
                vector<double> xi = ones;
		//lapack subroutine:  http://www.netlib.org/lapack/double/dposv.f 
                clapack_dposv(CblasRowMajor,CblasUpper,N,nrhs,&M1[0], N,&xi[0],N);
                vector<double> zeta = z;
                clapack_dposv(CblasRowMajor,CblasUpper,N,nrhs,&M[0], N,&zeta[0],N);
                beta = cblas_ddot(N, &ones[0],incx,&zeta[0],incx)/cblas_ddot(N, &ones[0],incx,&xi[0],incx);
		cblas_dgemv(CblasRowMajor,CblasNoTrans,N, N,1.0,&identity[0],N,&zeta[0],incx,-beta,&xi[0],incx);//xi=I*zeta-beta*xi
                alphas = xi;
                step = (beta_old-beta)*(beta_old-beta);
                beta_old = beta;
        }while((step > stepEpsilon || step <0.0) && (++iter < 100) );
}

//calculate all peps of vector xx and stored in vector predict. from BaseSpline.cpp
void KernelLogisticRegression::predict(const vector<double> &xx, 
					vector<double> & predict){
	predict.clear();
	transform( xx.begin(),xx.end(),back_inserter(predict),Predictor(this));
} 

double KernelLogisticRegression::predict(double xx){
        xx = transf(xx);
	size_t N=x.size();
        double gx=beta,ker1=0;
        for(int i =0; i<N; i++){
                ker1= kernel(xx, x[i],bandwidth);
                gx = gx+alphas[i]*ker1;
        }
        double pep = gx;
        return pep;
}


