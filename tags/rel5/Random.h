#ifndef RANDOM_H
#define RANDOM_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <ctime>
//
// Local copy of what we need from TRandom, extended
//
namespace RND {
  class Random {
  protected:
    unsigned int m_seed;  //Random number generator seed
  public:
    Random(unsigned int seed=65539);
    virtual ~Random();
    void copy(const Random &rnd) { m_seed = rnd.getSeed(); }
    Random & operator=(const Random &rnd) { copy(rnd); return *this; }
    //
    double gauss(double mean=0.0, double sigma=1.0);
    double logNormal(double mean=1.0, double sigma=1.0);
    double logNormalLN(double logMean, double logSigma);
    int    poisson(double mean);
    double flat(double mean, double sigma=1.0);
    double flatRange(double xmin, double xmax);
    double gamma(double mean, double sigma);
    double general(int npts, double *x, double xmin, double xmax, double *f, double fmin, double fmax);
    const unsigned int getSeed() const {return m_seed;}
    void   setSeed(unsigned int seed=65539);
    //
    virtual double rndm();
  };
  
  
#ifndef RANDOM_CXX
  extern Random gRandom;
#endif
};

#endif
