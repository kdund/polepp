#ifndef TABULATED_H
#define TABULATED_H

#include <iostream>
#include <cmath>

class TabFun {
public:
  TabFun();
  virtual ~TabFun();
  //
  void init(double xmin, double xmax, int n);
  virtual double F(double x) = 0;
  virtual double getVal(double x);
  //
protected:
  int    m_nData;
  double m_xMin;
  double m_xMax;
  double m_dx;
  double *m_data;
};

class TabTrig:public TabFun {
public:
  TabTrig();
  virtual ~TabTrig();
  //
  virtual double F(double x);
  virtual double Cos(double x);
  virtual double Sin(double x);
};

class TabLog:public TabFun {
public:
  TabLog();
  virtual ~TabLog();
  //
  virtual double F(double x);
  virtual double Log(double x);
};

class Poisson {
public:
  Poisson();
  ~Poisson();
  //
  void init(int nlambda=10000, int nn=51, double lmbmax=100.0);
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
  unsigned long sind  = static_cast<int>(s/m_dL);
  unsigned long index = no+sind*m_nN;
  if (index<m_nTot) {
    rval = m_data[index];
  } else {
    rval = rawPoisson(no,s);
  }
  return rval;
}

// inline double LogNormal::getVal(double x) {
//   return (x>0.0 ? Gaus(log(x))*exp(-x):0.0);
// }

#endif

