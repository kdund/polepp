#ifndef PDF_H
#define PDF_H

#include <iostream>
#include <cmath>

namespace PDF {

class Poisson {
public:
  Poisson();
  ~Poisson();
  //
  void init(int nlambda=100000, int nn=60, double lmbmax=200.0);
  inline double getVal(int no, double s);
private:
  unsigned long m_nLambda;
  unsigned long m_nN;
  unsigned long m_nTot;
  double m_lambdaMax;
  double m_dL;
  double *m_data;
  //
  double rawPoisson(int n, double s);
};

class Gauss {
public:
  Gauss();
  ~Gauss();
  //
  void init(int ndata=10000, double nmumax=100.0);
  inline double getVal(double x, double mu, double s);
  // 2D gauss
  inline double getVal2D(double x1, double mu1, double s1, double x2, double mu2, double s2, double corr);
  inline double getVal2D(double x1, double mu1, double x2, double mu2, double sdetC, double seff1, double seff2, double veffc);
  inline double getDetC(double s1,double s2,double c) { return (s1*s1*s2*s2*(1.0-c*c)); }
  inline double getVeff(double detC, double s) { return (detC/(s*s)); }
  inline double getVeffCorr(double detC, double s1, double s2, double corr) { return (detC/(corr*s1*s2)); }
  // Log-Normal
  inline double getValLogN(double x, double nmean, double nsigma);
  inline double getLNMean(double mean,double sigma);
  inline double getLNSigma(double mean,double sigma);
private:
  unsigned long m_nData;
  double m_muMax;
  double m_dMu;
  double *m_data;
  //
  double rawGauss(double mu);
};

// class General {
// public:
//   General();
//   General(const char *name);
//   General(int n, double *x, double *f);
//   ~General();
//   //
//   void load(const char *name);
//   void init();
//   //
//   double getVal(double x);
// private:
//   std::double vector m_x;
//   std::double vector m_f;
//   std::double vector m_acc;
// };

inline double Gauss::getVal(double x, double mu0, double s) {
  double rval;
  double mu = fabs((x-mu0)/s); // symmetric around mu0
  if (mu>m_muMax) {
    rval = rawGauss(mu)/s;
  } else {
    unsigned long sind  = static_cast<int>(mu/m_dMu);
    if (sind>=m_nData) {
      rval = rawGauss(mu)/s;
    } else {
      rval = m_data[sind]/s;
    }
  }
  return rval;
}

inline double Gauss::getVal2D(double x1, double mu1, double s1, double x2, double mu2, double s2, double corr) {
  double sdetC = sqrt(getDetC(s1,s2,corr));
  double seff1 = sdetC/s2;
  double seff2 = sdetC/s1;
  double veffc = sdetC*sdetC/(corr*s1*s2);
  //
  double rval;
  rval  = getVal(x1,mu1,seff1)*seff1;
  rval *= getVal(x2,mu2,seff2)*seff2;
  rval *= exp((x1-mu1)*(x2-mu2)/veffc);
  rval *= 1.0/sdetC;
  return rval;
}

inline double Gauss::getVal2D(double x1, double mu1, double x2, double mu2, double sdetC, double seff1, double seff2, double veffc) {
  double rval;
  rval  = getVal(x1,mu1,seff1)*seff1;
  rval *= getVal(x2,mu2,seff2)*seff2;
  rval *= exp((x1-mu1)*(x2-mu2)/veffc);
  rval *= 1.0/sdetC;
  return rval;
}

inline double Gauss::getValLogN(double x, double nmean, double nsigma) {
  return (x>0.0 ? getVal(log(x),nmean,nsigma)/x:0.0);
}
inline double Gauss::getLNMean(double mean,double sigma) {
  return log(mean*mean/sqrt(sigma*sigma + mean*mean));
}
inline double Gauss::getLNSigma(double mean,double sigma) {
  return sqrt(log(sigma*sigma/mean/mean+1));
}

inline double Poisson::getVal(int no, double s) {
  double rval;
  double ests;
  unsigned long sind  = static_cast<int>(s/m_dL);
  unsigned long index = no+sind*m_nN;
  double df;
  //
  if (index<m_nTot) {
    rval = m_data[index];
    df = -1.0;
    if (sind>0) df = rval*((static_cast<double>(no)/s)-1.0);
    ests = static_cast<double>(sind)*m_dL;
    rval = df*(s-ests) + rval;
  } else {
    rval = rawPoisson(no,s);
  }
  return rval;
}

};

#endif
