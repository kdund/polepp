#ifndef MEASUREMENT_H
#define  MEASUREMENT_H
//
// One measurement consists of N(observed) plus a number of nuisance parameters.
// The signal is a function of these parameters and N. It is implemented in getSignal().
//
#include <string>

//! Distribution type of nuisance parameters
enum DISTYPE {
  DIST_NONE=0,   /*!< No distrubution */
  DIST_GAUS,     /*!< Gaussian */
  DIST_FLAT,     /*!< Flat */
  DIST_LOGN,     /*!< Log-Normal */
  DIST_GAUSCORR  /*!< Correlated gauss (eff,bkg) */
};
/*!
  Returns a string corresponding to the given DISTYPE.
 */
inline const std::string distTypeStr(DISTYPE dt) {
  std::string rval;
  switch (dt) {
  case DIST_NONE:
    rval = "None";
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
  case DIST_GAUSCORR:
    rval = "GaussCorr";
    break;
  default:
    rval = "Unknown";
    break;
  }
  return rval;
}

class Measurement {
 public:
  Measurement():m_nObserved(0),
		m_effMeas(0),m_effSigma(0),m_effDist(DIST_NONE),
		m_bkgMeas(0),m_bkgSigma(0),m_bkgDist(DIST_NONE), m_beCorr(0), m_beCorrInv(0) {}
  Measurement(int nobs,
	      double effmeas, double effsigma, DISTYPE effdist,
	      double bkgmeas, double bkgsigma, DISTYPE bkgdist, double corr=0)
    :m_nObserved(nobs),m_effMeas(effmeas),m_effSigma(effsigma),m_effDist(effdist),
     m_bkgMeas(bkgmeas),m_bkgSigma(bkgsigma),m_bkgDist(bkgdist), m_beCorr(corr) { m_beCorrInv = sqrt(1.0-m_beCorr*m_beCorr);}

  Measurement(const Measurement & m) {}
  virtual ~Measurement() {}
  //
  void copy( const Measurement & m) {
    m_nObserved = m.getNObserved();
    m_effMeas   = m.getEffMeas();
    m_effSigma  = m.getEffSigma();
    m_effDist   = m.getEffDist();
    m_bkgMeas   = m.getBkgMeas();
    m_bkgSigma  = m.getBkgSigma();
    m_bkgDist   = m.getBkgDist();
    m_beCorr    = m.getBEcorr();
    m_beCorrInv = m.getBEcorrInv();
  }

  inline Measurement & operator=( const Measurement & m ) {
    if (this!=&m) {
      copy(m);
    }
    return *this;
  }

  inline bool operator==( const Measurement & m ) {
    return ( ( m_nObserved == m.getNObserved()) ||
	     ( m_effMeas   == m.getEffMeas()) ||
	     ( m_effSigma  == m.getEffSigma()) ||
	     ( m_effDist   == m.getEffDist()) ||
	     ( m_bkgMeas   == m.getBkgMeas()) ||
	     ( m_bkgSigma  == m.getBkgSigma()) ||
	     ( m_bkgDist   == m.getBkgDist()) ||
	     ( m_beCorr    == m.getBEcorr()) );
  }

  void setName(const char *name)   { m_name = name; }
  void setDescr(const char *descr) { m_description = descr; }
  void setEff(double effmeas, double effsigma, DISTYPE effdist) { m_effMeas = effmeas; m_effSigma = effsigma; m_effDist = effdist; }
  void setBkg(double bkgmeas, double bkgsigma, DISTYPE bkgdist) { m_bkgMeas = bkgmeas; m_bkgSigma = bkgsigma; m_bkgDist = bkgdist; }
  void setEffMeas( const double v )  { m_effMeas  = v; }
  void setEffSigma( const double v ) { m_effSigma = v; }
  void setEffDist( const DISTYPE v ) { m_effDist  = v; }
  void setBkgMeas( const double v )  { m_bkgMeas  = v; }
  void setBkgSigma( const double v ) { m_bkgSigma = v; }
  void setBkgDist( const DISTYPE v ) { m_bkgDist  = v; }
  void setBEcorr(double corr) { m_beCorr = corr; }
  void setNObserved(int nobs) { m_nObserved = nobs; }
  void dump() const;
  //
  const double  getEffMeas()   const { return m_effMeas; }
  const double  getEffSigma()  const { return m_effSigma; }
  const DISTYPE getEffDist()   const { return m_effDist; }
  const double  getBkgMeas()   const { return m_bkgMeas; }
  const double  getBkgSigma()  const { return m_bkgSigma; }
  const DISTYPE getBkgDist()   const { return m_bkgDist; }
  const double  getBEcorr()    const { return m_beCorr; }
  const double  getBEcorrInv() const { return m_beCorrInv; }
  const int     getNObserved() const { return m_nObserved; }
  //
  const double getSignal()     const { return (double(m_nObserved) - m_bkgMeas)/m_effMeas; }
  const double getSignalUnc()  const { 
    double s = getSignal();
    return sqrt(double(m_nObserved)+m_bkgSigma*m_bkgSigma+s*s*m_effSigma*m_effSigma)/m_effMeas;
  }
  //
  bool isFullyCorrelated(double eps=1e-16) const { return (((fabs(fabs(m_beCorr)-1.0)) < eps)); }
  bool isNotCorrelated(double eps=1e-16)   const { return (fabs(m_beCorr) < eps); }
  //
 protected:
  std::string m_name;
  std::string m_description;
  //
  int     m_nObserved;
  double  m_effMeas;
  double  m_effSigma;
  DISTYPE m_effDist;
  double  m_bkgMeas;
  double  m_bkgSigma;
  DISTYPE m_bkgDist;
  double  m_beCorr;
  double  m_beCorrInv;
};


#endif
