#include "Random.h"
/////////////////////////////////////////////////////////////////////

Random::Random(unsigned int seed) {
  SetSeed(seed);
}
Random::~Random() {
}

void Random::SetSeed(unsigned int seed) {
  //  Set the random generator seed
  //  if seed is zero, the seed is set to the current  machine clock
  //  Note that the machine clock is returned with a precision of 1 second.
  //  If one calls SetSeed(0) within a loop and the loop time is less than 1s, 
  //  all generated numbers will be identical!
   
  if( seed==0 ) {
    time_t curtime;      // Set 'random' seed number  if seed=0
    time(&curtime);      // Get current time in m_seed.
    m_seed = (unsigned int)curtime;
  } else {
    m_seed = seed;
  }
}

double Random::Rndm() {
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
  return Rndm();
}

int Random::Poisson(double mean)
{
  // Generates a random integer N according to a Poisson law.
  // Coded from Los Alamos report LA-5061-MS
  // Prob(N) = exp(-mean)*mean^N/Factorial(N)
  //
  if (mean <= 0) return 0;
  // use a gaussian approximation for large values of mean
  if (mean > 88) return static_cast<int>(Gaus(0,1)*sqrt(mean) + mean +0.5);
  //
  int N;
  double expmean = exp(-mean);
  double pir = 1;
  N = -1;
  while(1) {
    N++;
    pir *= Rndm();
    if (pir <= expmean) break;
  }
  return N;
}

double Random::Gaus(double mean, double sigma) {
  //      Return a number distributed following a gaussian with mean and sigma
  
  double x, y, z, result;
  
  y = Rndm();
  z = Rndm();
  x = z * 6.28318530717958623;
  result = mean + sigma*sin(x)*sqrt(-2*log(y));

  return result;
}

double Random::LogNormal(double mean, double sigma) {
  //      Return a number distributed following a lognormal with mean and sigma
  double nmean = log(mean*mean/sqrt(sigma*sigma + mean*mean));
  double nsigma = sqrt(log(sigma*sigma/mean/mean+1));
  return exp(Gaus() * nsigma + nmean);
}
