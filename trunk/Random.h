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
class Random {
protected:
  unsigned int m_seed;  //Random number generator seed
  double m_mean;
  double m_sigma;
public:
  Random(unsigned int seed=65539);
  virtual ~Random();
  double gaus(double mean=0.0, double sigma=1.0);
  double logNormal(double mean=1.0, double sigma=1.0);
  int    poisson(double mean);
  double flat(double mean, double sigma=1.0);
  unsigned int getSeed() {return m_seed;}
  void   setSeed(unsigned int seed=65539);
  //
  virtual double rndm();
};

#endif
