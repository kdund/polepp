#define RANDOM_CXX
#include "Random.h"
/////////////////////////////////////////////////////////////////////

namespace RND {
  Random gRandom;

  Random::Random(unsigned int seed) {
    setSeed(seed);
  }
  Random::~Random() {
  }

  void Random::setSeed(unsigned int seed) {
    //  Set the random generator seed
    //  if seed is zero, the seed is set to the current  machine clock
    //  Note that the machine clock is returned with a precision of 1 second.
    //  If one calls setSeed(0) within a loop and the loop time is less than 1s, 
    //  all generated numbers will be identical!
    
    if( seed==0 ) {
      time_t curtime;      // Set 'random' seed number  if seed=0
      time(&curtime);      // Get current time in m_seed.
      m_seed = (unsigned int)curtime;
    } else {
      m_seed = seed;
    }
  }
  
  double Random::rndm() {
    //  Machine independent random number generator.
    //  Produces uniformly-distributed floating points between 0 and 1.
    //  Identical sequence on all machines of >= 32 bits.
    //  Periodicity = 10**8
    //  Universal version (Fred James 1985).
    //  generates a number in ]0,1]
    
    const double kCONS = 4.6566128730774E-10;
    const int kMASK24  = 2147483392;
    
    m_seed *= 69069;  
    unsigned int jy = (m_seed&kMASK24); // Set lower 8 bits to zero to assure exact float
    if (jy) return kCONS*jy;
    return rndm();
  }
  
  int Random::poisson(double mean)
  {
    // Generates a random integer N according to a Poisson law.
    // Coded from Los Alamos report LA-5061-MS
    // Prob(N) = exp(-mean)*mean^N/Factorial(N)
    //
    if (mean <= 0) return 0;
    // use a gaussian approximation for large values of mean
    if (mean > 88) return static_cast<int>(gauss(0,1)*sqrt(mean) + mean +0.5);
    //
    int N;
    double expmean = exp(-mean);
    double pir = 1;
    N = -1;
    while(1) {
      N++;
      pir *= rndm();
      if (pir <= expmean) break;
    }
    return N;
  }
  
  double Random::gamma(double mean, double sigma) {
    //See algo at http://en.wikipedia.org/wiki/Gamma_distribution#Generating_Gamma_random_variables
  }

  double Random::gauss(double mean, double sigma) {
    //      Return a number distributed following a gaussian with mean and sigma
    
    double x, y, z, result;
    
    y = rndm();
    z = rndm();
    x = z * 6.28318530717958623;
    result = mean + sigma*sin(x)*sqrt(-2*log(y));
    
    return result;
  }
  
  double Random::logNormal(double mean, double sigma) {
    // Return a number distributed following a lognormal with mean and sigma
    double nmean = log(mean*mean/sqrt(sigma*sigma + mean*mean));
    double nsigma = sqrt(log(sigma*sigma/mean/mean+1));
    return exp(gauss() * nsigma + nmean);
  }

  double Random::logNormalLN(double logMean, double logSigma) {
    // Return a number distributed following a lognormal with mean and sigma
    return exp(gauss() * logSigma + logMean);
  }
  
  double Random::flat(double mean, double sigma) {
    double dx=sigma*1.73205081; // == sqrt(12.0)*0.5 => sigma(flat) = (xmax-xmin)/sqrt(12)
    double xmin = mean-dx;
    return rndm()*dx*2.0 + xmin;
  }
  
  double Random::general(int npts, double *xvec, double xmin, double xmax, double *fvec, double fmin, double fmax) {
    double dx = xmax-xmin;
    double du = fmax-fmin;
    double x;
    double xn;
    double u;
    double fx;
    double rval=0;
    bool accepted=false;
    
    int nlow = 0;
    int nup  = npts-1;
    int n;
    int indx;
    //
    while(!accepted) {
      x = rndm()*dx + xmin;
      u = rndm()*du + fmin;
      // search for index of given x
      nlow = 0;
      n = nlow+(int)(((double)nup)*(x-xmin)/dx); // First approximation -> is correct if equidistant points
      while ((nup-nlow>1)) {
	xn = xvec[n];
	if (xn>x) {
	  nup = n;
	} else {
	  nlow = n;
	}
	n = (nup+nlow)/2;
      }
      indx = nlow; // the closest index from below
      fx = fvec[indx];
      //
      if (u<fx) {
	accepted=true;
	rval=x;
      }
    }
    return rval;
  }
};
