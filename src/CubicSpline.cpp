#include <iostream>
#include <vector>
using namespace std;
#include "CubicSpline.h"

CubicSpline::CubicSpline():yp0(0.),ypn(0.),derivCalculated(false),natural0(true),naturaln(true),logged(false)
{
}

CubicSpline::~CubicSpline()
{
}


CubicSpline::CubicSpline(vector<double> &xx, vector<double> &yy, bool l) {
  setData(xx,yy,l);
}


void CubicSpline::setData(vector<double> &xx, vector<double> &yy, bool l) {
  logged = l;
  if(xx.size()!=yy.size()) {
    cerr << "Arrays x and y are of different length " << xx.size() 
         << " " << yy.size();
    exit(-1);
  }
  if(xx.size()<3) {
    cerr << "A minimum of three data points is needed";
    exit(-1);
  }
  x = xx;
  y = yy;
  d2y.assign(xx.size(),0.0);
  removeDuplicates();
  vector<double>::iterator val=y.begin();
  while(val!=y.end()) {
    if (*val<=0.0) *val = 1e-10;
    if (*val>=1.0) *val = 1-1e-10;
    *val=logit(*val);
  }
}

void CubicSpline::removeDuplicates() {
  vector<double>::iterator oldx,currentx,oldy,currenty;
  oldx = x.begin();
  oldy = y.begin();
  currentx = x.begin() + 1;
  currenty = y.begin() + 1;
  while (currentx<x.end()) {
    if (*oldx == *currentx && *oldy == *currenty) {
      x.erase(currentx);
      y.erase(currenty);
    } else {
      oldx++;oldy++;
      currentx++;currenty++;
    }
  }
}


void CubicSpline::calcDeriv() {
  if (derivCalculated) return;

  double	p=0.0,qn=0.0,sig=0.0,un=0.0;
  vector<double> u;
  u.assign(x.size(),0.0);

  if (natural0) {
    d2y[0]=u[0]=0.0;
  } else {
    d2y[0] = -0.5;
    u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp0);
  }

  for(unsigned int i=1;i<=x.size()-2;i++){
    sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
    p=sig*d2y[i-1]+2.0;
    d2y[i]=(sig-1.0)/p;
    u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
  }

  if (naturaln){
    qn=un=0.0;
  } else {
    qn=0.5;
    un=(3.0/(x[x.size()-1]-x[x.size()-2]))*(ypn-(y[y.size()-1]-y[y.size()-2])/(x[x.size()-1]-x[x.size()-2]));
  }

  d2y[x.size()-1]=(un-qn*u[x.size()-2])/(qn*d2y[x.size()-2]+1.0);
  for(int k=x.size()-2;k>=0;k--){
    d2y[k]=d2y[k]*d2y[k+1]+u[k];
  }
  derivCalculated = true;
}

double CubicSpline::interpolate(double xx){
  if (xx<x[0] || xx>x[x.size()-1]){
    cerr << "x (" << xx << ") is outside the range of data points (" << x[0] << " to " << x[x.size()-1]<< endl;
    exit(-1);
  }

  calcDeriv();

  double h=0.0,b=0.0,a=0.0, yy=0.0;
  int k=0;
  int klo=0;
  int khi=x.size()-1;
  while (khi-klo > 1){
    k=(khi+klo) >> 1;
    if(x[k] > xx){
      khi=k;
    }
    else{
      klo=k;
    }
  }
  h=x[khi]-x[klo];

  if (h == 0.0){
    cerr << "Two values of x are identical: point " << klo << " (" << x[klo] <<") and point " 
	 << khi << " (" << x[khi] << ")" << endl ;
    exit(-1);
  }
  else{
    a=(x[khi]-xx)/h;
    b=(xx-x[klo])/h;
    yy=a*y[klo]+b*y[khi]+((a*a*a-a)*d2y[klo]+(b*b*b-b)*d2y[khi])*(h*h)/6.0;
  }
  if (logged) 
    return invlogit(yy);
  return yy;
}

