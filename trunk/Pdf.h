#ifndef PDF_H
#define PDF_H
//
#include <iostream>
#include <vector>
#include <cmath>
#include "Random.h"

/*!
  
 */
#ifndef M_PIl
#define M_PIl 3.1415926535897932384626433832795029L
#endif



namespace PDF {
#ifdef PDF_CXX
  bool gPrintStat = false;
#else
  extern bool gPrintStat;
#endif
  enum DISTYPE {
    DIST_NONE=0,   /*!< No distrubution */
    DIST_POIS,     /*!< Poisson */
    DIST_GAUS,     /*!< Gaussian */
    DIST_FLAT,     /*!< Flat */
    DIST_LOGN,     /*!< Log-Normal */
    DIST_GAUS2D,   /*!< Correlated gauss (eff,bkg) */
    DIST_LAST,     /*!< Last element before UNDEF */
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
  // Help functions for Log Normal dist
  //
  inline const double calcLogMean(const double mean,const double sigma)  { return log(mean*mean/sqrt(sigma*sigma + mean*mean)); }
  inline const double calcLogSigma(const double mean,const double sigma) { return sqrt(log((sigma*sigma/(mean*mean))+1)); }
  //
  class Base {
  public:
    Base() { m_dist=DIST_UNDEF; m_statNraw = 0; m_statNtot=0; m_statNtab=0; }
    Base(const char *name, DISTYPE d=DIST_UNDEF, double m=0.0, double s=0.0)
      :m_dist(d), m_mean(m), m_sigma(s), m_statNraw(0), m_statNtot(0), m_statNtab(0)
    { if (name) m_name=name; setDist(d); }
    Base(const Base & other) { copy(other); m_statNraw=0; m_statNtot=0; m_statNtab=0; }
    virtual ~Base() { this->printStat(); }
    //
    virtual void setMean(double m)  { m_mean  = m; }
    virtual void setSigma(double s) { m_sigma = s; }
    //
    virtual const double getVal(const double x, const double mean, const double sigma) const {
      std::cerr << "ERROR: PDF::Base - Accessing getVal(x,m,s) - VERBOTEN!!!" << std::endl;
      exit(-1);
      return 0;
    }
    virtual const double getVal(const double x, const double mean) const {
      std::cerr << "ERROR: PDF::Base - Accessing getVal(double x,m) - VERBOTEN!!!" << std::endl;
      exit(-1);
      return 0;
    }
    virtual const double getVal(const int x, const double mean, const double sigma) const {
      std::cerr << "ERROR: PDF::Base - Accessing getVal(int x,m,s) - VERBOTEN!!!" << std::endl;
      exit(-1);
      return 0;
    }
    virtual const double getVal(const int x, const double mean) const {
      std::cerr << "ERROR: PDF::Base - Accessing getVal(int x,m) - VERBOTEN!!!" << std::endl;
      exit(-1);
      return 0;
    }
    const std::string & getName() const { return m_name;}
    const DISTYPE getDist()      const { return m_dist;}
    const double  getMean()      const { return m_mean;}
    const double  getSigma()     const { return m_sigma;}

    virtual void printStat() const {
      if (gPrintStat) {
        std::cout << "----- " << this->getName() << " -----" << std::endl;
        std::cout << " PDF called                    : " << this->m_statNtot << std::endl;
        std::cout << " PDF called using raw function : " << this->m_statNraw << std::endl;
        std::cout << " PDF called using table        : " << this->m_statNtab << std::endl;
        std::cout << "--------------------------------------------------" << std::endl;
      }
    }
  protected:
    void copy(const Base & other) {
      if (this != &other) {
	m_name  = "Copy of " + other.getName();
	m_dist  = other.getDist();
	m_mean  = other.getMean();
	m_sigma = other.getSigma();
      }
    }
    void setDist(DISTYPE d) { if (d>=DIST_LAST) d=DIST_UNDEF; m_dist = d; }
    void setName(std::string name) { m_name = name; }
    std::string m_name;
    DISTYPE     m_dist;
    double      m_mean;
    double      m_sigma;

    mutable int m_statNraw;
    mutable int m_statNtot;
    mutable int m_statNtab;
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
    virtual ~BaseType() { }
    //
    virtual const double F(T x) const=0;
    virtual const double getVal(const T x, const double mean, const double sigma) const {std::cerr << "ERROR: Accessing getVal(T x,m,s) - NOT IMPLEMENTED for this class" << std::endl; exit(-1); return 0;}
    virtual const double getVal(const T x, const double mean) const {std::cerr << "ERROR: Accessing getVal(T x,m) - NOT IMPLEMENTED for this class" << std::endl; exit(-1); return 0;}
    inline const double operator()(T x) const { return F(x); }

  };

  class Flat : public BaseType<double> {
  public:
    Flat():BaseType<double>("Flat",DIST_FLAT,1.0,0.1) { setMinMax(getMean(),getSigma()); }
    Flat(double mean, double sigma):BaseType<double>("Flat",DIST_FLAT,mean,sigma) { setMinMax(mean,sigma); }
    Flat(const Flat & other):BaseType<double>(other) { m_min = other.getMin(); m_max = other.getMax(); m_F = other.getF();}
    virtual ~Flat() {};
    //
    inline const double F(double val) const;
    inline const double getVal(const double x, const double mean, const double sigma) const;
    //
    void setMean(  const double m)  { Base::setMean(m);  setMinMax(m,getSigma()); }
    void setSigma( const double m)  { Base::setSigma(m); setMinMax(getMean(),m); }

    inline const double getMin() const { return m_min; }
    inline const double getMax() const { return m_max; }
    inline const double getF()   const { return m_F; }
    
  protected:
    double m_min;
    double m_max;
    double m_F;

    inline const double raw(const double x, const double f) const;
    inline const double raw(const double x, const double f, const double xmin, const double xmax) const;

    inline void setMinMax(double mean, double sigma) {
      calcMinMax(mean,sigma,m_min,m_max);
      m_F = calcF(m_min,m_max);
      //      std::cout << "Flat: " << m_F << " , " << m_min << " , " << m_max << std::endl;
    }

    inline void calcMinMax(double mean, double sigma, double & xmin, double & xmax) const {
      const double d = sqrt(12.0);
      double wh = 0.5*d*sigma;
      xmin = mean-wh;
      xmax = mean+wh;
      if (xmax<xmin) xmax=xmin;
    }

    inline const double calcF(double xmin, double xmax) const {
      double d=(xmax-xmin);
      return (d>0 ? 1.0/d:-1);
    }
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
    inline const double getVal(const double x, const double mean, const double sigma) const;
  };

  class Gauss2D : public Gauss {
  public:
    Gauss2D():Gauss() { this->m_name="Gauss2D"; this->m_dist=DIST_GAUS2D; }
    Gauss2D(const Gauss2D & other):Gauss(other) {}
    virtual ~Gauss2D() {};
    //
    inline const double getVal2D(const double x1, const double mu1, const double s1, const double x2, const double mu2, const double s2, const double corr) const;
    inline const double getVal2D(const double x1, const double mu1, const double x2, const double mu2, const double sdetC, const double seff1, const double seff2, const double veffc) const;
    inline const double getDetC(const double s1,const double s2,const double c) const { return (s1*s1*s2*s2*(1.0-c*c)); }
    inline const double getVeff(const double detC, const double s) const { return (detC/(s*s)); }
    inline const double getVeffCorr(const double detC, const double s1, const double s2, const double corr) const { return ((corr*s1*s2)/detC);}
  };

  class LogNormal : public Gauss {
  public:
    LogNormal():Gauss()                                    { m_name="LogNormal"; m_dist=DIST_LOGN; setMean(1.0);  setSigma(1.0); }
    LogNormal(const double mean, const double sigma):Gauss(mean,sigma) { m_name="LogNormal"; m_dist=DIST_LOGN; setMean(mean); setSigma(sigma);}
    LogNormal(const LogNormal & other):Gauss(other) {
      m_mean = other.getMean(); m_sigma = other.getSigma();
      m_logMean = other.getLogMean(); m_logSigma = other.getLogSigma();
    }
    virtual ~LogNormal() {};
    //
    void setMean( const double m)  { this->m_mean = m;  m_logMean = calcLogMean(m,this->m_sigma); m_logSigma = calcLogSigma(m,this->m_sigma); }
    void setSigma( const double m) { this->m_sigma = m; m_logMean = calcLogMean(this->m_mean,m);  m_logSigma = calcLogSigma(this->m_mean,m); }
    //

    inline const double getLogMean()  const { return m_logMean; }
    inline const double getLogSigma() const { return m_logSigma; }


    inline const double F(const double x) const {return (x>0.0 ? Gauss::getVal(log(x),m_logMean,m_logSigma)/x:0.0);}
    inline const double getVal(const double x, const double m, const double s) const {
      m_statNtot++;
      if (x<=0) return 0.0;
      m_statNraw++;
      return Gauss::getVal(log(x),calcLogMean(m,s), calcLogSigma(m,s))/x;
    }
    inline const double getValLogN(const double x, const double m, const double s) const {
      m_statNraw++;
      return Gauss::getVal(x, m, s)/exp(x);
    }
  protected:
    double m_logMean;
    double m_logSigma;
  };

  class Poisson : public BaseType<int> {
  public:
    Poisson():BaseType<int>("Poisson",DIST_POIS,1.0,1.0) {}
    Poisson(const double lambda):BaseType<int>("Poisson",DIST_POIS,lambda,sqrt(lambda)) {}
    Poisson(const Poisson & other):BaseType<int>(other) {}

    virtual ~Poisson() {}
    //
    void setMean(const double mean)   { this->m_mean = mean; this->m_sigma = sqrt(mean); }
    void setSigma(const double sigma) { this->m_mean = sigma*sigma; this->m_sigma = sigma; }
    //
    virtual inline const double F(int val) const;
    virtual inline const double getVal(const int x, const double mean, const double sigma) const;
    inline const double getVal(const int x, const double mean) const;
    virtual inline const double getVal(const double x, const double mean, const double sigma) const;
    inline const double getVal(const double x, const double mean) const;
    inline const double raw(const int n, const double s) const;
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
    virtual ~Tabulated() { if (m_table) delete [] m_table;};
    //
    virtual void setMean(double m)  { m_pdf->setMean(m);  Base::setMean(m_pdf->getMean());}
    virtual void setSigma(double s) { m_pdf->setSigma(s); Base::setMean(m_pdf->getSigma());}
    //
    virtual void setRangeX(     int npts, T      min, T      max) { m_nX     = npts; m_xmin = min; m_xmax = max; }
    virtual void setRangeMean(  int npts, double min, double max) { m_nMean  = npts; m_mmin = min; m_mmax = max; }
    virtual void setRangeSigma( int npts, double min, double max) { m_nSigma = npts; m_smin = min; m_smax = max; }
    virtual void setBasePdf( BaseType<T> * pdf ) {
      m_pdf = pdf;
      Base::setMean(m_pdf->getMean());
      Base::setSigma(m_pdf->getSigma());
      Base::setDist(m_pdf->getDist());
      Base::setName("Tabulated " + m_pdf->getName());
    }
    
    void initTab() {
      if ((m_nX<1) || (m_nMean<1) || (m_nSigma<1)) return;
      std::cout << "---------------------------------------------------------------------" << std::endl;
      std::cout << "--- Tabulating pdf <" << m_pdf->getName() << ">" << std::endl;
      std::cout << "--- N(X)      = " << m_nX << std::endl;
      std::cout << "---   min      = " << m_xmin << std::endl;
      std::cout << "---   max      = " << m_xmax << std::endl;
      std::cout << "--- N(mean)    = " << m_nMean << std::endl;
      std::cout << "---   min      = " << m_mmin << std::endl;
      std::cout << "---   max      = " << m_mmax << std::endl;
      std::cout << "--- N(sigma)   = " << m_nSigma << std::endl;
      std::cout << "---   min      = " << m_smin << std::endl;
      std::cout << "---   max      = " << m_smax << std::endl;
      std::cout << "---------------------------------------------------------------------" << std::endl;
      
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
      std::cout << "--- Tabulating ... be patient" << std::endl;
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
      std::cout << "--- Tabulating DONE!" << std::endl;
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
	if (ind<m_nTotal) {
          m_statNtab++;
	  return m_table[ind];
        }
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
      setBasePdf(pdf);
//       this->m_pdf = pdf;
//       m_dist = DIST_POIS;
//       this->m_mean  = pdf->getMean();
//       this->m_sigma = pdf->getSigma();
      m_nSigma = 1;
    }
    virtual ~PoisTab() { }
    //
    // X range for poisson is integer -> force one point per integer -> xmax = xmin + npts - 1
    //
    virtual void setRangeX( int npts, int min, int max=0) {
      if (npts<1) npts=0;
      m_nX = npts;
      m_xmin = min;
      m_xmax = min+npts-1;
    }
    virtual void setRangeSigma( int npts, double min, double max) { m_nSigma = 1; }
    //
    void tabulate() {
      if (this->m_pdf==0) return;
      initTab();
      if (m_table==0) return;
      double mean;
      int n0;
      int indn0;
      int ind0;
      std::cout << "--- Tabulating ... be patient" << std::endl;
      for (int m=0; m<m_nMean; m++) {
	mean = m*m_dm+m_mmin;
        // start at n0 = int(mean) -> important if mean is very large
        n0 = static_cast<int>(mean);
        if (n0<m_xmin) n0=m_xmin;
        if (n0>m_xmax) n0=m_xmax;
        indn0 = n0-m_xmin;
	ind0 = m*m_nX;
        // set reference point at maximum
	m_table[ind0+indn0] = static_cast<Poisson *>(this->m_pdf)->raw(n0,mean);

        // set for all n>n0
	for (int n=indn0+1; n<m_nX; n++) {
	  m_table[ind0+n] = m_table[ind0+n-1]*mean/static_cast<double>(n+m_xmin);
	}
        // set for all n<n0
	for (int n=indn0-1; n>=0; n--) {
	  m_table[ind0+n] = m_table[ind0+n+1]*static_cast<double>(n+m_xmin+1)/mean;
        }
      }
      std::cout << "--- Tabulating DONE!" << std::endl;
    }
    //
    void setMean(double mean)   { this->m_pdf->setMean(mean);   this->m_mean = mean;   this->m_sigma = this->m_pdf->getSigma(); }
    void setSigma(double sigma) { this->m_pdf->setSigma(sigma); this->m_sigma = sigma; this->m_mean  = this->m_pdf->getMean(); }

    virtual const double getVal(int x, double m) const {
      this->m_statNtot++;
      //
      // check if table is created and that the requested values are within the table
      // if not, the value will be calculated using raw()
      //
      if ((m_table!=0) &&
	  (x>=this->m_xmin) && (x<=this->m_xmax) &&
	  (m>=this->m_mmin) && (m<=this->m_mmax)) {
        //
        // calculate index in table
        //
	int  mind, xind, ind;
	mind = int(this->m_dm>0 ? (m-this->m_mmin)/this->m_dm : 0);
	xind = x-m_xmin;
	ind = xind + mind*this->m_nX;
        //
        // Make sure index is OK (should always be or else BUG!
        //
	if (ind<this->m_nTotal) {
          double lmb0 = double(mind)*this->m_dm + m_mmin;
          double dlmb = m - lmb0;
          double f0   = this->m_table[ind];
          //
          // Do Taylor expansion around lmb0 to second order
          //
          double alpha=0.0;
          double beta=0.0;
          if (lmb0>0.0) {
            alpha = (double(x)/lmb0)-1.0;
            beta  = double(x)/(lmb0*lmb0);
          }
          double corr1 = f0*alpha*dlmb;
          double corr2 = 0.5*f0*(alpha*alpha - beta)*dlmb*dlmb;
          this->m_statNtab++;
	  return f0 + corr1 + corr2;
        }
        // This line should never be called - a BUG trap
        std::cout << "PoisTab::Failed for whatever reason! BUG! mean = " << m << " and x = " << x << std::endl;
      }
      //
      // OK. Call the pdf, but first check if it's defined
      //
      if (this->m_pdf==0) {
	std::cerr << "ERROR in PDF::PoisTab - no pdf defined!" << std::endl;
	return 0;
      }
      //
      // Call the raw() function
      //
      this->m_statNraw++;
      return this->m_pdf->getVal(x,m,0); // Poisson ignores sigma
    }
    //
    virtual const double getVal(int x, double m, double s) const {
      return this->getVal(x,m);
    }
    virtual const double getVal(double x, double m, double s) const {
      return this->getVal(int(x+0.5),m);
    }
  };
    
  class GaussTab : public Tabulated<double> {
  public:
    GaussTab():Tabulated<double>() {}
    GaussTab(Gauss *pdf):Tabulated<double>() {
      this->m_pdf = pdf;
      this->m_dist = DIST_GAUS;
      this->m_mean  = pdf->getMean();
      this->m_sigma = pdf->getSigma();
    }
    virtual ~GaussTab() {}
    //
    void tabulate() {
      if (this->m_pdf==0) return;
      m_nSigma = 1; // force them to be unity - tabulate only for N(0,1)
      m_nMean  = 1;
      initTab();
      if (m_table==0) return;
      double x;
      for (int n=0; n<m_nX; n++) {
	x = double(n)*m_dx+m_xmin;
	m_table[n] = static_cast<Gauss *>(this->m_pdf)->phi(x);
      }
    }

    virtual const double getVal(double x, double m, double s) const {
      m_statNtot++;
      if (m_table!=0) {
	double mu = fabs((x-m)/s);
	if (mu>m_xmax)
	  return static_cast<Gauss *>(this->m_pdf)->phi(mu)/s;
	int muind = int(m_dx>0 ? (mu-m_xmin)/m_dx : 0);
	if (muind<m_nTotal) {
          m_statNtab++;
	  return m_table[muind];
        }
      }
      //
      if (this->m_pdf==0) {
        std::cout << "FATAL: undefined PDF!" << std::endl;
        exit(-1);
        return 0;
      }
      m_statNraw++;
      return this->m_pdf->getVal(x,m,s);
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
    double mu = fabs((x-this->m_mean)/this->m_sigma); // symmetric around mu0
    return phi(mu)/this->m_sigma;
  }
  inline const double Gauss::getVal(const double x, const double mean, const double sigma) const {
    m_statNraw++;
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
    return raw(x,this->m_mean);
  }
  inline const double Poisson::getVal(const int x, const double mean) const {
    return raw(x,mean);
  }
  inline const double Poisson::getVal(const int x, const double mean, const double sigma) const {
    return raw(x,mean);
  }
  inline const double Poisson::getVal(const double x, const double mean) const {
    return raw(int(x+0.5),mean);
  }
  inline const double Poisson::getVal(const double x, const double mean, const double sigma) const {
    return raw(int(x+0.5),mean);
  }

  inline const double Poisson::raw(const int n, const double s) const {
    m_statNtot++;
    m_statNraw++;
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


  inline const double Flat::F(double x) const {
    return raw(x,m_F);
  }
  inline const double Flat::getVal(const double x, const double mean, const double sigma) const {
    //    std::cout << "Flat: getVal() " << std::endl;
    double xmin,xmax,f;
    calcMinMax(mean,sigma,xmin,xmax);
    f = calcF(xmin,xmax);
    double rval = raw(x,f,xmin,xmax);
    //    std::cout << " : x = " << x << " , f = " << f << " , range = " << xmin << " - " << xmax << " ==> " << rval << std::endl;
    return rval;
  }

  inline const double Flat::raw(const double x, const double f) const {
    m_statNraw++;
    return (((x>=m_min) && (x<=m_max)) ? f:0);
  }

  inline const double Flat::raw(const double x, const double f, const double xmin, const double xmax) const {
    m_statNraw++;
    return (((x>=xmin) && (x<=xmax)) ? f:0);
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
//       this->m_pdf->setMean(m);
//       this->m_pdf->setSigma(s);
//       rval = this->m_pdf->F(x);
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
  extern Flat      gFlat;
#endif
};

#endif
