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
  inline double getVal(double s, double mu, double sigma);
private:
  unsigned long m_nData;
  double m_muMax;
  double m_dMu;
  double *m_data;
  //
  double rawGauss(double mu);
};

inline double Gauss::getVal(double s, double mu0, double sigma) {
  double rval;
  double mu = fabs((s-mu0)/sigma); // symmetric around mu0
  if (mu>m_muMax) {
    rval = rawGauss(mu)/sigma;
  } else {
    unsigned long sind  = static_cast<int>(mu/m_dMu);
    if (sind>=m_nData) {
      rval = rawGauss(mu)/sigma;
    } else {
      rval = m_data[sind]/sigma;
    }
  }
  return rval;
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

