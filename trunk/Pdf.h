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
  //
  inline const double getVal(int no, double s) const;
  const double *getData() const { return m_data; }
  const int getNdata() const    { return m_nTot; }
private:
  unsigned long m_nLambda;
  unsigned long m_nN;
  unsigned long m_nTot;
  double m_lambdaMax;
  double m_dL;
  double *m_data;
  //
  const double rawPoisson(int n, double s) const;
};

class Gauss {
public:
  Gauss();
  ~Gauss();
  //
  void init(int ndata=10000, double nmumax=100.0);
  inline const double getVal(double x, double mu, double s) const;
  // 2D gauss
  inline const double getVal2D(double x1, double mu1, double s1, double x2, double mu2, double s2, double corr) const;
  inline const double getVal2D(double x1, double mu1, double x2, double mu2, double sdetC, double seff1, double seff2, double veffc) const;
  inline const double getDetC(double s1,double s2,double c) const { return (s1*s1*s2*s2*(1.0-c*c)); }
  inline const double getVeff(double detC, double s) const { return (detC/(s*s)); }
  inline const double getVeffCorr(double detC, double s1, double s2, double corr) const { return (detC/(corr*s1*s2));}
  // Log-Normal
  inline const double getValLogN(double x, double nmean, double nsigma) const;
  inline const double getLNMean(double mean,double sigma) const;
  inline const double getLNSigma(double mean,double sigma) const;
private:
  unsigned long m_nData;
  double m_muMax;
  double m_dMu;
  double *m_data;
  //
  const double rawGauss(double mu) const;
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

inline const double Gauss::getVal(double x, double mu0, double s) const {
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

inline const double Gauss::getVal2D(double x1, double mu1, double s1, double x2, double mu2, double s2, double corr) const {
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

inline const double Gauss::getVal2D(double x1, double mu1, double x2, double mu2, double sdetC, double seff1, double seff2, double veffc) const {
  double rval;
  rval  = getVal(x1,mu1,seff1)*seff1;
  rval *= getVal(x2,mu2,seff2)*seff2;
  rval *= exp((x1-mu1)*(x2-mu2)/veffc);
  rval *= 1.0/sdetC;
  return rval;
}

inline const double Gauss::getValLogN(double x, double nmean, double nsigma) const {
  return (x>0.0 ? getVal(log(x),nmean,nsigma)/x:0.0);
}
inline const double Gauss::getLNMean(double mean,double sigma) const {
  return log(mean*mean/sqrt(sigma*sigma + mean*mean));
}
inline const double Gauss::getLNSigma(double mean,double sigma) const {
  return sqrt(log(sigma*sigma/mean/mean+1));
}

inline const double Poisson::getVal(int no, double s) const {
  double rval;
  double ests;
  unsigned long sind  = static_cast<int>(s/m_dL);
  unsigned long index = no+sind*m_nN;
  double df,df2,ndl;
  //
  if (index<m_nTot) {
    rval = m_data[index];
    ests = static_cast<double>(sind)*m_dL;
    df = -(s-ests);
    df2=  0.0;
    if (sind>0) {
      ndl = static_cast<double>(no)/ests;
      df = rval*(ndl-1.0)*(s-ests);
      if (rval<0.01)
	df2 = 0.5*rval*((ndl-1.0)*(ndl-1.0) - (1.0+(ndl/ests)))*(s-ests)*(s-ests);
    }
    rval += df + df2;
  } else {
    rval = rawPoisson(no,s);
  }
  return rval;
}

extern Poisson gPoisson;
extern Gauss   gGauss;

};

#endif
