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
  unsigned int   m_seed;  //Random number generator seed
  double m_mean;
  double m_sigma;
public:
  Random(unsigned int seed=65539);
  virtual ~Random();
  double Gaus(double mean=0.0, double sigma=1.0);
  double LogNormal(double mean=1.0, double sigma=1.0);
  int    Poisson(double mean);
  double Flat(double mean, double sigma=1.0);
  unsigned int   GetSeed() {return m_seed;}
  void   SetSeed(unsigned int seed=65539);
  //
  virtual double Rndm();

};

#endif
