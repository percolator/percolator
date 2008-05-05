#include<vector>
#include<utility>
#include<cstdlib>
#include<fstream>
#include<iterator>
#include<algorithm>
#include<numeric>
#include<functional>
using namespace std;
#include "Option.h"
#include "ArrayLibrary.h"
#include "LogisticRegression.h"
#include "PosteriorEstimator.h"

static unsigned int noIntevals = 500;
//static unsigned int noIntevals = 20;
static unsigned int numLambda = 20;
static double maxLambda = 0.9;

PosteriorEstimator::PosteriorEstimator()
{
}

PosteriorEstimator::~PosteriorEstimator()
{
}

pair<double,bool> make_my_pair(double d,bool b) {
  return make_pair(d,b);
}

class IsDecoy {public: bool operator() (const pair<double,bool>& aPair) 
    {return not aPair.second; }};

bool isMixed(const pair<double,bool>& aPair) {
  return aPair.second;
}

template<class T> void bootstrap(const vector<T>& in, vector<T>& out) {
  out.clear();
  double n = in.size();
  for (size_t ix=0;ix<n;++ix) {
    size_t draw = (size_t)((double)rand()/((double)RAND_MAX+(double)1)*n);
    out.push_back(in[draw]);
  }
  // sort in desending order
  sort(out.begin(),out.end());  
} 


void PosteriorEstimator::estimate( vector<pair<double,bool> >& combined, LogisticRegression& lr, double pi0) {
  // sorting in accending order
  sort(combined.begin(),combined.end());

  vector<double> medians;
  vector<unsigned int> negatives,sizes;
  
  binData(combined,medians,negatives,sizes);

  lr.setData(medians,negatives,sizes);
  lr.iterativeReweightedLeastSquares();

  // sorting combined in score decending order
  reverse(combined.begin(),combined.end());
  lr.setCutOff(pi0);
}

// Estimates q-values and prints
void PosteriorEstimator::finishStandalone(const vector<pair<double,bool> >& combined, LogisticRegression& lr, double pi0) { 
  vector<double> q;
  getQValues(pi0,combined,q);
  
  sort(combined.begin(),combined.end(),greater<pair<double,bool> >());

  for(size_t ix=0;ix<combined.size();++ix) {
    if (not (combined[ix].second))
      continue;
    double pred = lr.predict(combined[ix].first);
    cout << combined[ix].first << "\t" << pred << "\t" << q[ix] << endl;
  }
}


void PosteriorEstimator::binData(const vector<pair<double,bool> >& combined,  
  vector<double>& medians, vector<unsigned int>& negatives, vector<unsigned int>& sizes) {
  // Create bins and count number of negatives in each bin
  double binSize = max(combined.size()/(double)noIntevals,1.0); 
  vector<pair<double,bool> >::const_iterator combinedIter=combined.begin();  

  unsigned int binNo =0;
  size_t firstIx, pastIx = 0; 

  while (pastIx<combined.size()) {
    firstIx = pastIx; pastIx = min(combined.size(),(size_t)((++binNo)*binSize)); 
    int inBin = pastIx-firstIx;
    int negInBin = count_if(combinedIter,combinedIter+inBin,IsDecoy());
    combinedIter += inBin;
    double median = combined[firstIx+inBin/2].first;
    if (medians.size()>0 and *(medians.rbegin())==median) {
      *(sizes.rbegin()) += inBin;
      *(negatives.rbegin()) += negInBin;
    } else {
       medians.push_back(median);
       sizes.push_back(inBin);
       negatives.push_back(negInBin);
    }    
  }
  
}


double mymin(double a,double b) {return a>b?b:a;}

void PosteriorEstimator::getQValues(double pi0,
     const vector<pair<double,bool> >& combined, vector<double>& q) {
  // assuming combined sorted in decending order
  vector<pair<double,bool> >::const_iterator myPair = combined.begin();
  unsigned int nTargets = 0, nDecoys = 0;
  while(myPair != combined.end()) {
    if (!(myPair->second)) {
      ++nDecoys; 
    } else {
      ++nTargets;
      q.push_back(((double)nDecoys)/(double)nTargets);
    }
    ++myPair;
  }
  transform(q.begin(), q.end(), q.begin(), bind2nd(multiplies<double>(), pi0*((double)nTargets/(double)nDecoys)));
  partial_sum(q.rbegin(), q.rend(), q.rbegin(), mymin);
  return;  
}


void PosteriorEstimator::getPValues(
     const vector<pair<double,bool> >& combined, vector<double>& p) {
  // assuming combined sorted in decending order
  vector<pair<double,bool> >::const_iterator myPair = combined.begin();
  unsigned int nDecoys = 0;
  while(myPair != combined.end()) {
    if (!(myPair->second)) {
      ++nDecoys; 
    } else {
      p.push_back((double)nDecoys);
    }
    ++myPair;
  }
  transform(p.begin(), p.end(), p.begin(), bind2nd(divides<double>(), (double)nDecoys));
    
  return;  
}

double PosteriorEstimator::estimatePi0(vector<pair<double,bool> >& combined, 
                                       const unsigned int numBoot) {
  
  vector<double> p,pBoot,lambdas,pi0s,mse;
  vector<double>::iterator start;
  
  getPValues(combined,p);
  
  size_t n = p.size();
  // Calculate pi0 for different values for lambda    
  for(unsigned int ix=0; ix <= numLambda; ++ix) {
    double lambda = ((ix+1)/(double)numLambda)*maxLambda;
    lambdas.push_back(lambda);
    start = lower_bound(p.begin(),p.end(),lambda);
    double Wl = (double) distance(start,p.end());
    double pi0 = Wl/n/(1-lambda);
    pi0s.push_back(pi0);
  }
     
  double minPi0 = *min_element(pi0s.begin(),pi0s.end()); 
  fill_n(back_inserter(mse),pi0s.size(),0.0); 

  // Examine witch lambda level that is most stable under bootstrap   
  for (unsigned int boot = 0; boot< numBoot; ++boot) {
    bootstrap<double>(p,pBoot); 
    for(unsigned int ix=0; ix < lambdas.size(); ++ix) {
      start = lower_bound(pBoot.begin(),pBoot.end(),lambdas[ix]);
      double Wl = (double) distance(start,pBoot.end());
      double pi0Boot = Wl/n/(1-lambdas[ix]);
      mse[ix] += (pi0Boot-minPi0)*(pi0Boot-minPi0);
    }
  }   
  unsigned int minIx = distance(mse.begin(),min_element(mse.begin(),mse.end()));
  double pi0 = min(pi0s[minIx],1.0); 
  
  cerr << "pi_0=" << pi0 << endl;
   
  return pi0;
}

void PosteriorEstimator::run() {
  ifstream target(targetFile.c_str(),ios::in),decoy(decoyFile.c_str(),ios::in);
  istream_iterator<double> tarIt(target),decIt(decoy);

  // Merge a labeled version of the two lists into a combined list
  vector<pair<double,bool> > combined;  
  transform(tarIt,istream_iterator<double>(),back_inserter(combined),
            bind2nd(ptr_fun(make_my_pair),true)); 
  transform(decIt,istream_iterator<double>(),back_inserter(combined),
            bind2nd(ptr_fun(make_my_pair), false)); 
  
  // Estimate pi0
  double pi0 = estimatePi0(combined);
  
  // Logistic regression on the data
  LogisticRegression lr;
  estimate(combined,lr,pi0);

  finishStandalone(combined,lr,pi0);  
}

bool PosteriorEstimator::parseOptions(int argc, char **argv){
  // init
  CommandLineParser cmd("Posterior Estimation");
  // finally parse and handle return codes (display help etc...)
  cmd.parseArgs(argc, argv);

  if (cmd.arguments.size()>2) {
      cerr << "Too many arguments given" << endl;
      cmd.help();
  }
  if (cmd.arguments.size()==0) {
      cerr << "No arguments given" << endl;
      cmd.help();
  }
  targetFile = cmd.arguments[0];
  decoyFile = cmd.arguments[1];
  return true;
}

