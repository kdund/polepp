#ifndef PDF_H
#define PDF_H
//
#include <iostream>
#include <vector>
#include <cmath>
#include "Random.h"

/*!
  
 */
namespace PDF {
  enum DISTYPE {
    DIST_NONE=0,   /*!< No distrubution */
    DIST_POIS,     /*!< Poisson */
    DIST_GAUS,     /*!< Gaussian */
    DIST_FLAT,     /*!< Flat */
    DIST_LOGN,     /*!< Log-Normal */
    DIST_GAUS2D,   /*!< Correlated gauss (eff,bkg) */
    DIST_UNDEF=999 /*!< No distrubution defined*/
  };
  /*!
    Returns a string corresponding to the given DISTYPE.
  */
  inline const std::string distTypeStr(DISTYPE dt) {
    std::string rval;
    switch (dt) {
    case DIST_UNDEF:
      rval = "Undefined";
      break;
    case DIST_NONE:
      rval = "None";
      break;
    case DIST_POIS:
      rval = "Poisson";
      break;
    case DIST_GAUS:
      rval = "Gauss";
      break;
    case DIST_FLAT:
      rval = "Flat";
      break;
    case DIST_LOGN:
      rval = "LogN";
      break;
    case DIST_GAUS2D:
      rval = "Gauss 2D";
      break;
    default:
      rval = "Unknown";
      break;
    }
    return rval;
  }
  //
  class Base {
  public:
    Base() { m_dist=DIST_UNDEF; }
    Base(const char *name, DISTYPE d=DIST_UNDEF, double m=0.0, double s=0.0):m_dist(d), m_mean(m), m_sigma(s) { if (name) m_name=name; }
    Base(const Base & other) { copy(other); }
    virtual ~Base() {};
    //
    virtual void setMean(double m)  { m_mean  = m; }
    virtual void setSigma(double s) { m_sigma = s; }
    //
    const std::string & getName() const { return m_name;}
    const DISTYPE getDist()      const { return m_dist;}
    const double  getMean()      const { return m_mean;}
    const double  getSigma()     const { return m_sigma;}
    
  protected:
    void copy(const Base & other) {
      if (this != &other) {
	m_name  = "Copy of " + other.getName();
	m_dist  = other.getDist();
	m_mean  = other.getMean();
	m_sigma = other.getSigma();
      }
    }
    void setDist(DISTYPE d) { m_dist = d; }
    void setName(std::string name) { m_name = name; }
    std::string m_name;
    DISTYPE     m_dist;
    double      m_mean;
    double      m_sigma;
  };

  template <typename T>
  class BaseType: public Base {
  public:
    BaseType():Base() {}
    BaseType(const char *name, DISTYPE d=DIST_NONE, double m=0.0, double s=0.0):Base(name,d,m,s) {}
    BaseType(const BaseType<T> & other):Base() {
      if (this!=&other) {
	Base::copy(other);
      }
    }
    virtual ~BaseType() {}
    //
    virtual const double F(T x) const=0;
    virtual const double getVal(T x, double mean, double sigma) const=0;
    inline const double operator()(T x) const { return F(x); }
    
  };

  class Gauss : public BaseType<double> {
  public:
    Gauss():BaseType<double>("Gaussian",DIST_GAUS,0.0,1.0) {}
    Gauss(double mean, double sigma):BaseType<double>("Gaussian",DIST_GAUS,mean,sigma) {}
    Gauss(const Gauss & other):BaseType<double>(other) {}
    virtual ~Gauss() {};
    //
    inline const double F(double val) const;
    inline const double phi(double mu) const;
    inline const double getVal(double x, double mean, double sigma) const;
  };

  class Gauss2D : protected Gauss {
  public:
    Gauss2D():Gauss() { m_name="Gauss2D"; m_dist=DIST_GAUS2D; }
    Gauss2D(const Gauss2D & other):Gauss(other) {}
    virtual ~Gauss2D() {};
    //
    inline const double getVal2D(double x1, double mu1, double s1, double x2, double mu2, double s2, double corr) const;
    inline const double getVal2D(double x1, double mu1, double x2, double mu2, double sdetC, double seff1, double seff2, double veffc) const;
    inline const double getDetC(double s1,double s2,double c) const { return (s1*s1*s2*s2*(1.0-c*c)); }
    inline const double getVeff(double detC, double s) const { return (detC/(s*s)); }
    inline const double getVeffCorr(double detC, double s1, double s2, double corr) const { return ((corr*s1*s2)/detC);}
  };

  class LogNormal : protected Gauss {
  public:
    LogNormal():Gauss(1.0,1.0)                             { m_name="LogNormal"; m_dist=DIST_LOGN; }
    LogNormal(double mean, double sigma):Gauss(mean,sigma) { m_name="LogNormal"; m_dist=DIST_LOGN; }
    LogNormal(const LogNormal & other):Gauss(other) {}
    virtual ~LogNormal() {};
    //
    void setMean( double m)  { m_mean = m;  m_logMean = calcLogMean(m,m_sigma); m_logSigma = calcLogSigma(m,m_sigma); }
    void setSigma( double m) { m_sigma = m; m_logMean = calcLogMean(m_mean,m);  m_logSigma = calcLogSigma(m_mean,m); }
    //
    inline const double calcLogMean(double mean,double sigma)  const { return log(mean*mean/sqrt(sigma*sigma + mean*mean)); }
    inline const double calcLogSigma(double mean,double sigma) const { return sqrt(log((sigma*sigma/(mean*mean))+1)); }
    inline const double getLogMean()  const { return m_logMean; }
    inline const double getLogSigma() const { return m_logSigma; }


    inline const double F(double x) const {return (x>0.0 ? Gauss::getVal(log(x),m_logMean,m_logSigma)/x:0.0);}
    inline const double getVal(double x, double m, double s) const {
      if (x<=0) return 0.0;
      return Gauss::getVal(log(x),calcLogMean(m,s), calcLogSigma(m,s))/x;
    }
    inline const double getValLogN(double x, double m, double s) const {
      if (x<=0) return 0.0;
      return Gauss::getVal(log(x),m, s)/x;
    }
  protected:
    double m_logMean;
    double m_logSigma;
  };

  class Poisson : public BaseType<int> {
  public:
    Poisson():BaseType<int>("Poisson",DIST_POIS,1.0,1.0) {}
    Poisson(double lambda):BaseType<int>("Poisson",DIST_POIS,lambda,sqrt(lambda)) {}
    Poisson(const Poisson & other):BaseType<int>(other) {}

    virtual ~Poisson() {}
    //
    void setMean(double mean)   { m_mean = mean; m_sigma = sqrt(mean); }
    void setSigma(double sigma) { m_mean = sigma*sigma; m_sigma = sigma; }
    //
    virtual inline const double F(int val) const;
    virtual inline const double getVal(int x, double mean, double sigma) const;
    inline const double getVal(int x, double mean) const;
    inline const double raw(int n, double s) const;
  protected:

  };
   
  template <typename T>
  class Tabulated : public BaseType<T> {
  public:
    Tabulated():BaseType<T>() { m_pdf=0; m_table=0; }
    Tabulated(const Tabulated<T> & other):BaseType<T>() {
      if (this != &other) {
	BaseType<T>::copy(other);
	m_pdf  = other.getPdf();
	
	m_xmin = other.getXmin();
	m_xmax = other.getXmax();
	m_dx   = other.getdX();
	m_nX   = other.getNX();
	
	m_mmin = other.getMmin();
	m_mmax = other.getMmax();
	m_dm   = other.getdM();
	m_nMean= other.getNM();
	
	m_smin   = other.getSmin();
	m_smax   = other.getSmax();
	m_ds     = other.getdS();
	m_nSigma = other.getNS();
	m_nTotal = other.getNTot();
	
	m_table = new double[m_nTotal];
	double *otab = other.getTable();
	if (otab) {
	  for (int i=0; i<m_nTotal; i++) {
	    m_table[i] = otab[i];
	  }
	}
      }
    }
    virtual ~Tabulated() {if (m_table) delete [] m_table;};
    //
    virtual void setMean(double m)  { m_pdf->setMean(m);  Base::setMean(m_pdf->getMean());}
    virtual void setSigma(double s) { m_pdf->setSigma(s); Base::setMean(m_pdf->getSigma());}
    //
    void setRangeX(     int npts, T      min, T      max) { m_nX     = npts; m_xmin = min; m_xmax = max; }
    void setRangeMean(  int npts, double min, double max) { m_nMean  = npts; m_mmin = min; m_mmax = max; }
    void setRangeSigma( int npts, double min, double max) { m_nSigma = npts; m_smin = min; m_smax = max; }
    void setBasePdf( BaseType<T> * pdf ) {
      m_pdf = pdf;
      Base::setMean(m_pdf->getMean());
      Base::setSigma(m_pdf->getSigma());
      Base::setDist(m_pdf->getDist());
      Base::setName("Tabulated " + m_pdf->getName());
    }
    
    void initTab() {
      if ((m_nX<1) || (m_nMean<1) || (m_nSigma<1)) return;
      m_nTotal = m_nX*m_nMean*m_nSigma;
      m_table = new double[m_nTotal];
      m_dx=0;
      m_ds=0;
      m_dm=0;
      if (m_nX>1)     m_dx = (m_xmax-m_xmin)/(m_nX-1);
      if (m_nMean>1)  m_dm = (m_mmax-m_mmin)/(m_nMean-1);
      if (m_nSigma>1) m_ds = (m_smax-m_smin)/(m_nSigma-1);
    }
    void clearTable() { if (m_table) delete [] m_table; m_table=0; m_nTotal=0;}
    //
    virtual void tabulate() {
      initTab();
      if (m_table==0) return;
      int ind=0;
      for (int ns=0; ns<m_nSigma; ns++) {
	double s = double(ns)*m_ds+m_smin;
	if (m_pdf->getDist()!=DIST_POIS)
	  m_pdf->setSigma(s); // skip this for poisson
	for (int nm=0; nm<m_nMean; nm++) {
	  double m = double(nm)*m_dm+m_mmin;
	  m_pdf->setMean(m);
	  for (int nx=0; nx<m_nX; nx++) {
	    T x = nx*m_dx+m_xmin;
	    m_table[ind] = m_pdf->F(x);
	    ind++;
	  }
	}
      }
    }

    virtual const double F(T x) const {
      return m_pdf->F(x);
    }

    virtual const double getVal(T x, double m, double s) const {
      if ((m_table!=0) &&
	  (x>=m_xmin) && (x<=m_xmax) &&
	  (m>=m_mmin) && (m<=m_mmax) &&
	  (s>=m_smin) && (s<=m_smax) ) {
	int sind, mind, xind, ind;
	sind = int(m_ds>0 ? (s-m_smin)/m_ds : 0);
	mind = int(m_dm>0 ? (m-m_mmin)/m_dm : 0);
	xind = int(m_dx>0 ? (x-m_xmin)/m_dx : 0);
	ind = xind + mind*m_nX + sind*m_nX*m_nMean;
	if (ind<m_nTotal)
	  return m_table[ind];
      }
      //
      if (m_pdf==0) {
	std::cerr << "ERROR in PDF::Tabulated - no pdf defined!" << std::endl;
	return 0;
      }
      return m_pdf->getVal(x,m,s);
    }

    const BaseType<T> *getPdf()  const { return m_pdf;  }
    const T        getXmin() const { return m_xmin; }
    const T        getXmax() const { return m_xmax; }
    const T        getdX()   const { return m_dx; }
    const double   getMmin() const { return m_mmin; }
    const double   getMmax() const { return m_mmax; }
    const double   getdM()   const { return m_dm; }
    const double   getSmin() const { return m_smin; }
    const double   getSmax() const { return m_smax; }
    const double   getdS()   const { return m_ds; }
    const int      getNX()   const { return m_nX; }
    const int      getNM()   const { return m_nMean; }
    const int      getNS()   const { return m_nSigma; }
    const int      getNTot() const { return m_nTotal; }

    const double *getTable() const { return m_table; }

  protected:
    BaseType<T> *m_pdf;
    int   m_nTotal;
    int   m_nX;
    int   m_nMean;
    int   m_nSigma;
    //
    T      m_xmin, m_xmax, m_dx;
    double m_mmin, m_mmax, m_dm;
    double m_smin, m_smax, m_ds;

    double *m_table;
  };

  class PoisTab : public Tabulated<int> {
  public:
    PoisTab():Tabulated<int>() {};
    PoisTab(Poisson *pdf):Tabulated<int>() {
      m_pdf = pdf;
      m_dist = DIST_POIS;
      m_mean  = pdf->getMean();
      m_sigma = pdf->getSigma();
    }
    virtual ~PoisTab() {};
    //
    void tabulate() {
      if (m_pdf==0) return;
      initTab();
      if (m_table==0) return;
      double mean;
      int ind0;
      for (int m=0; m<m_nMean; m++) {
	mean = m*m_dm+m_mmin;
	ind0 = m*m_nX + 0;
	m_table[ind0] = static_cast<Poisson *>(m_pdf)->raw(0,mean);
	for (int n=1; n<m_nX; n++) {
	  m_table[ind0+n] = m_table[ind0+n-1]*mean/static_cast<double>(n);
	}
      }
    }
    //
    void setMean(double mean)   { m_pdf->setMean(mean);   m_mean = mean;   m_sigma = m_pdf->getSigma(); }
    void setSigma(double sigma) { m_pdf->setSigma(sigma); m_sigma = sigma; m_mean  = m_pdf->getMean(); }

    virtual const double getVal(int x, double m) const {
      if ((m_table!=0) &&
	  (x>=m_xmin) && (x<=m_xmax) &&
	  (m>=m_mmin) && (m<=m_mmax)) {
	int  mind, xind, ind;
	mind = int(m_dm>0 ? (m-m_mmin)/m_dm : 0);
	xind = int(m_dx>0 ? (x-m_xmin)/m_dx : 0);
	//	ind = xind + sind*m_nX + mind*m_nX*m_nSigma;
	ind = xind + mind*m_nX;
	if (ind<m_nTotal)
	  return m_table[ind];
      }
      //
      if (m_pdf==0) {
	std::cerr << "ERROR in PDF::PoisTab - no pdf defined!" << std::endl;
	return 0;
      }
      return m_pdf->getVal(x,m,0); // Poisson ignores sigma
    }
    //
    virtual const double getVal(int x, double m, double s) const {
      return getVal(x,m);
    }
  };
    
  class GaussTab : public Tabulated<double> {
  public:
    GaussTab():Tabulated<double>() {}
    GaussTab(Gauss *pdf):Tabulated<double>() {
      m_pdf = pdf;
      m_dist = DIST_GAUS;
      m_mean  = pdf->getMean();
      m_sigma = pdf->getSigma();
    }
    virtual ~GaussTab() {}
    //
    void tabulate() {
      if (m_pdf==0) return;
      m_nSigma = 1; // force them to be unity - tabulate only for N(0,1)
      m_nMean  = 1;
      initTab();
      if (m_table==0) return;
      double x;
      for (int n=0; n<m_nX; n++) {
	x = double(n)*m_dx+m_xmin;
	m_table[n] = static_cast<Gauss *>(m_pdf)->phi(x);
      }
    }

    virtual const double getVal(double x, double m, double s) const {
      if (m_table!=0) {
	double mu = fabs((x-m)/s);
	if (mu>m_xmax)
	  return static_cast<Gauss *>(m_pdf)->phi(mu)/s;
	int muind = int(m_dx>0 ? (mu-m_xmin)/m_dx : 0);
	if (muind<m_nTotal)
	  return m_table[muind];
      }
      //
      if (m_pdf==0) return 0;
      return m_pdf->getVal(x,m,s);
    }
  };
  
  // class General : public Base<double> {
  // public:
  //   General():Base<double>("General") {}
  //   General(std::vector<double> & x, std::vector<double> & f):Base<double>("General") { setData(x,f); init(); }
  //   //  Gauss(const Gauss & other):<double>(other) {m_mean=other.getMean(); m_sigma=other.getSigma(); m_name=other.getName();}
  //   virtual ~Gauss() {};
  //   //
  //   void setData(std::vector<double> & x, std::vector<double> & f) { m_x = x; m_f = f; }
  //   void init();
  //   double getMean()  { return m_mean; }
  //   double getSigma() { return m_sigma; }
  //   //
  //   inline double F(double val);
  //   inline double phi(double mu);
  // private:
  //   std::vector<double> m_x;
  //   std::vector<double> m_f;
  //   double m_mean;
  //   double m_sigma;
  // };

  inline const double Gauss::phi(double mu) const {
    return (1.0L/sqrt(2.0*M_PIl))*exp(-0.5L*mu*mu);
  }
  inline const double Gauss::F(double x) const {
    double mu = fabs((x-m_mean)/m_sigma); // symmetric around mu0
    return phi(mu)/m_sigma;
  }
  inline const double Gauss::getVal(double x, double mean, double sigma) const {
    double mu = fabs((x-mean)/sigma); // symmetric around mu0
    return phi(mu)/sigma;
  }

  inline const double Gauss2D::getVal2D(double x1, double mu1, double s1, double x2, double mu2, double s2, double corr) const {
    double sdetC = sqrt(getDetC(s1,s2,corr));
    double seff1 = sdetC/s2;
    double seff2 = sdetC/s1;
    double veffc = (corr*s1*s2)/sdetC*sdetC;
    //
    double rval;
    rval  = getVal(x1,mu1,seff1)*seff1;
    rval *= getVal(x2,mu2,seff2)*seff2;
    rval *= exp((x1-mu1)*(x2-mu2)*veffc);
    rval *= 1.0/sdetC;
    return rval;
  }
  
  inline const double Gauss2D::getVal2D(double x1, double mu1, double x2, double mu2, double sdetC, double seff1, double seff2, double veffc) const {
    double rval;
    rval  = getVal(x1,mu1,seff1)*seff1;
    rval *= getVal(x2,mu2,seff2)*seff2;
    rval *= exp((x1-mu1)*(x2-mu2)*veffc);
    rval *= 1.0/sdetC;
    return rval;
  }
  inline const double Poisson::F(int x) const {
    return raw(x,m_mean);
  }
  inline const double Poisson::getVal(int x, double mean) const {
    return raw(x,mean);
  }
  inline const double Poisson::getVal(int x, double mean, double sigma) const {
    return raw(x,mean);
  }

  inline const double Poisson::raw(int n, double s) const {
    double prob;
    double nlnl,lnn,lnf;
    prob = 0.0;
    nlnl = double(n)*log(s);  // n*ln(s)
    lnn  = lgamma(n+1);       // ln(fac(n))
    lnf  = nlnl - lnn - s;
    if (isinf(lnf) || isnan(lnf)) {
      prob=(n==0 ? 1.0:0.0);
    } else {
      prob=exp(lnf);
    }
    if (isnan(prob)) {
      std::cout << "NaN in rawPoisson: " << n << ", " << s << ", " << prob << std::endl;
    }
    return prob;
  }
//   template <typename T>
//   const double Tabulated::getVal(T x, double m, double s) {
//     double rval;
//     if (!((m_table==0) ||
// 	  (x<xmin) || (x>xmax) ||
// 	  (m<mmin) || (m>mmax) ||
// 	  (s<smin) || (s>smax))) {
//       int sind, mind, xind, ind;
//       sind = (m_ds>0 ? (s-m_smin)/m_ds : 0);
//       mind = (m_dm>0 ? (m-m_mmin)/m_dm : 0);
//       xind = (m_dx>0 ? (x-m_xmin)/m_dx : 0);
//       ind = m_nX + sind*m_nX + mind*m_nSigma;
//       if (ind<m_nTotal)
// 	rval = m_table[ind];
//     } else {
//       m_pdf->setMean(m);
//       m_pdf->setSigma(s);
//       rval = m_pdf->F(x);
//     }
//     return rval;
//   }


//
#ifndef PDF_CXX
  extern Poisson  gPoisson;
  extern PoisTab  gPoisTab;
  extern Gauss    gGauss;
  extern GaussTab gGaussTab;

  extern Gauss2D   gGauss2D;
  extern LogNormal gLogNormal;
#endif
};

#endif
