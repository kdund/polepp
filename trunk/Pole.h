#ifndef POLE_H
#define POLE_H
//
// The idea is to seperate out Pole-specific code from the Coverage class.
// class Pole
// ----------
// User input:
//   * Measured N(observed), efficiency and background
//   * Distribution of efficiency and background - include the possibility
//     to have arbitrary distributions (a 
//   * Optional: an array containing the values of the folded PDF (the double integral)
//   * Optional: true signal - use it to check if it's in- or outside the CL - <Coverage>
// Output:
//   * Lower and upper limit
//   * Optional: if true signal is given, then stop limit search when it's clear
//     that the true signal is in- or outside the range.
// Methods:
//   * Generate weight table (from dists of nuasance-parameters) - POSSIBILITY TO REUSE/SAVE
//   * Generate s_best and R(s,n) - reuse if eff/bkg not changed
//
// class Coverage
// --------------
//
// Input:
//   * Ranges of signal, efficiency, background
//   * Distributions of eff,bkg and correlation
//   * Later: An arbitrary number of observables with their distributions (class Observable)
//
//
//
//

#include <iostream>
#include <vector>
#include <cmath>
#include <ctime>
#include <utility>
#include <map>
#include <string>
#include <list>

#if __GNUC__ < 3
#include <cstdio>
#endif

#include "Range.h"
#include "Pdf.h"
#include "BeltEstimator.h"

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

/*! @class Pole
 *
 *  @brief The main class containing the methods for calculating limits.
 *  
 *  This class calculates confidence intervals for
 *  a Poisson process with background using a frequentist confidence belt
 *  construction.
 *  It assumes that the measurement can be described as follows:
 *  \f[N_{obs} = \epsilon s + b\f]
 *  where
 *  - \f$N_{obs}\f$ = number of observed events (Poisson)
 *  - \f$\epsilon\f$ = measured efficiency with a known (or assumed) distribution
 *  - \f$s\f$ = unknown signal for which the limits are to be found
 *  - \f$b\f$ = measured background with a known (or assumed) distribution
 *
 *  The PDF used for describing the number of observed events is:  
 *  \f[
 *  q(n)_{s+b} =  
 *  \frac{1}{2\pi\sigma_{b}\sigma_{\epsilon}} \times
 *  \intop_0^{\infty}\intop_0^{\infty}p(n)_{b+ 
 *  \epsilon's}\;\;
 *  f_{b,\sigma_b}(b')\;\;
 *  g_{\epsilon,\sigma_{\epsilon}}(\epsilon')\;\;
 *  db'd\epsilon'
 *  \f]
 *  where \f$f()\f$ and \f$g()\f$ are the PDF for the background and efficiency respectively.
 *
 *  In order to find the lower and upper limits, a so called confidence belt is constructed.
 *  For each value of \f$s\f$, \f$n_1\f$ and \f$n_2\f$ are found such that
 *  \f$\sum_{n=n_1}^{n_2} q(n)_{s+b} = 1 - \alpha\f$ where \f$1-\alpha\f$ is the confidence limit.
 *  The confidence belt is then given by the set \f$[n_1(s+b,\alpha),n_2(s+b,\alpha)]\f$.
 *  An upper and lower limit is found by finding the intersection of the vertical line \f$n = N_{obs}\f$
 *  and the boundaries of the belt (not exactly true due to the discreteness of a Poisson).
 *
 *  However, the confidence belt is not unambigously defined. A specific ordering scheme is required.
 *  The method used for selecting \f$n_1\f$ and
 *  \f$n_2\f$ is based on likelihood ratios (aka Feldman & Cousins). For each n, a \f$s_{best}\f$ is found that
 *  will maximise the likelihood \f$\mathcal{L}(n)_{s+b}\f$ (mathematically it's just \f$q(n)_{s+b}\f$,
 *  but here it's not used as a PDF but rather as a hypothesis). For a fixed s, a likelihood ratio is calculated
 *  \f[R(n,s)_{\mathcal{L}} = \frac{\mathcal{L}_{s+b}(n)}{\mathcal{L}_{s_{best}+b}(n)}\f]
 *  and the n are included in the sum starting with the highest rank (ratio) and continuing with decreasing rank until
 *  the sum (ref) equals the requested confidence.
 *
 *  \b SETUP
 *
 *  Basic
 *  - setCL() : Confidence limit, default 0.9
 *    the requested confidence [0.0,1.0]
 *  - setNobserved() : Number of observed events
 *  - setEffMeas() : Measured efficiency\n
 *    efficiency distribution (mean,sigma and distribution type ( DISTYPE ))
 *  - setBkgMeas() : Measured background\n
 *    background distribution (mean,sigma and distribution type ( DISTYPE ))
 *  - setEffBkgCorr() : Correlation between eff and bkg (if applicable)
 *    correlation coefficient [-1.0,1.0]
 *
 *  Integration
 *  - setEffInt() : Integration range of efficiency\n
 *    The integration range must cover the PDF such that the tails are
 *    negligable.
 *  - setBkgInt() : Dito, background
 *
 *  Belt construction
 *  - setBelt() : Set the maximum n in the belt construction\n
 *    For large signals and/or events, this value might be increased. 
 *    If the nBelt < 1, then this is automatically selected (see suggestBelt()).
 *  - findBelt() : Calculates the confidence belt [n1(s,b),n2(s,b)] for all (s,b).
 *  - calcBelt() : Dito but for a specific (s,b)
 *
 *  Finding \f$s_{best}\f$
 *  - setDmus() : Sets the precision in findBestMu().\n
 *    Default = 0.01 and it should normally be fine.
 *  - setNLR() : If true, use an alternate method (Gary Hill)
 *
 *  Hypothesis testing
 *  - setTestHyp() : test range for s+b\n
 *    This sets the range and step size (precision) for finding the limits.\n
 *    The default range is [0.0,35.0] with a step size of 0.01.
 *
 *  Coverage related
 *  - setTrueSignal() : Set the true signal.
 *  - setCoverage() : Set the coverage flag.\n
 *    Setting this flag causes the limit calculation to terminate the scan as soon as it is
 *    decided whether the true signal is inside or outside the confidence limit.
 *
 *  Tabulated PDFs
 *  - initPoisson() : Initialises a tabulated poisson.\n
 *    It is not required but it greatly speeds up the coverage calculations.
 *  - initGauss() : Dito but for a gauss pdf.\n
 *    The gain in speed is much less significant (REMOVE???)
 *
 *  \b RUNNING
 *
 *  Print
 *  - printSetup() : Prints the setup to stdout.
 *  - printLimit() : Prints the calculated limit.
 *
 *  General
 *  - analyseExperiment() : Calculates the limit using the current setup.
 *
 *  Debug
 *  - setVerbose() : Sets verbose level.
 *
 *  \b EXAMPLES
 *  see polelim.cxx (limit calculation) and polecov.cxx (coverage study)
 *
 * @author Fredrik Tegenfeldt (fredrik.tegenfeldt@cern.ch)
 */
//  e^{\frac{-(b-b')^2}{2\sigma_b^2}}\;\;
//  e^{\frac{- (1 -\epsilon')^2}{2\sigma_{\epsilon}^2}}
class Pole {
public:
  Pole();
  ~Pole();
  // Main parameters
  void setCL(double cl)    { m_cl = cl; }
  void setNobserved(int nobs) {m_nObserved = nobs;}

  //! distribution info on eff and bkg
  void setEffMeas(double mean,double sigma, DISTYPE dist=DIST_GAUS);
  void setBkgMeas(double mean,double sigma, DISTYPE dist=DIST_GAUS);
  void setEffBkgCorr(double corr) {m_beCorr = corr;}

  //
  bool checkEffBkgDists();
  bool isFullyCorrelated() { return (((fabs(fabs(m_beCorr)-1.0)) < 1e-16)); }
  bool isNotCorrelated()   { return (fabs(m_beCorr) < 1e-16); }

  // POLE construction
  void setEffInt(double scale=0,double step=0); // efficiency range for integral (7) [AFTER the previous]
  void setBkgInt(double scale=0,double step=0); // dito background

  // Checks the integration range for eff and bkg
  bool checkParams();

  // Set belt max value
  void setBelt(int v)    { m_nBelt = v; m_suggestBelt = (v<1); }
  void setBeltMax(int v) { m_nBeltMax = v; } //allocated maximum
  int  suggestBelt();                // will suggest a m_nBelt based on no. observed
  void setDmus(double dmus) { m_dmus = (dmus > 0.0 ? dmus:m_stepMin); }

  // POLE test hypothesis range
  void setTestHyp(double low, double high, double step); // test mu=s+b for likelihood ratio
  void setNuppLim(int n=-1) { m_nUppLim = n; }

  // POLE true signal, used only if coverage run
  void setTrueSignal(double s) { m_sTrue = s; } // true signal
  void setCoverage(bool flag) {m_coverage = flag;} // true if coverage run

  // Ordering scheme - if true it will use an optional likelihood ratio (Gary Hill)
  void setNLR(bool flag) { m_useNLR=flag;}
  bool usesNLR() {return m_useNLR;}

  // Debug
  void setVerbose(int v=0) { m_verbose=v; }

  // Tabulated PDF's
  void initPoisson(int nlambda=10000, int nn=51, double lmbmax=100); // will init poisson table
  void initGauss(int ndata  =10000, double mumax=1000.0); // will init gauss table
  ///////////////////////////////
  //
  void initIntArrays();   // will initialise integral arrays (if needed)
  void initBeltArrays();  // will initialise belt arrays (if needed)
  void initIntegral();    // calculates double integral kernal (eff*bkg*db*de) according to setup (7)

  // POLE
  inline double calcProb(int n, double s); // calculates probability (7)
  void findBestMu(int n); // finds the best fit (mu=s+b) for a given n. Fills m_bestMu[n] and m_bestMuProb[n].
  void findAllBestMu();   // dito for all n (loop n=0; n<m_nMuUsed)
  void calcConstruct(double s);
  double calcBelt(double s, int & n1, int & n2); // calculate (4) and find confidence belt
  double calcLimit(double s); // calculate (4) and find limits, returns probability for given signal hypothesis
  bool limitsOK(); // check if calculated limit is OK using the sum of probs.
  inline bool normOK(double p);
  void setNormMaxDiff(double dpmax=0.001) { m_normMaxDiff=dpmax; }
  void findConstruct();
  void findBelt();
  bool findLimits();        // finds CL limits
  bool findCoverageLimits();//  same as previous but stop if it's obvious that initial true mean is inside or outside
  void printLimit(bool doTitle=false);
  // 
  void printSetup();
  // POLE
  bool analyseExperiment();  // makes the double integral and finds the limits/coverage
  //
  // Access functions
  //
  double getStepMin()    { return m_stepMin; }
  double getCL()         { return m_cl; }
  double getSTrue()      { return m_sTrue; }
  bool   getCoverage()   { return m_coverage; }
  int    getNObserved()  { return m_nObserved; }
  // Efficiency
  double getEffMeas()    { return m_effMeas; }
  double getEffSigma()   { return m_effSigma; }
  DISTYPE getEffDist()   { return m_effDist; }
  // Background
  double  getBkgMeas()   { return m_bkgMeas; }
  double  getBkgSigma()  { return m_bkgSigma; }
  DISTYPE getBkgDist()   { return m_bkgDist; }
  double  getEffBkgCorr() { return m_beCorr; }
  // range and steps in double integral (7), in principle infinite
  double  getEffIntScale() { return m_effIntScale; }
  Range  *getEffRangeInt() { return &m_effRangeInt; }
  // range and steps in double integral (7)
  double  getBkgIntScale() { return m_bkgIntScale; }
  Range  *getBkgRangeInt() { return &m_bkgRangeInt; }
  // Test range for the likelihood ratio calculation (4)
  Range  *getHypTest()     { return &m_hypTest; }
  //
  bool    isValidInt()   { return m_validInt; }
  int     getNInt()      { return m_nInt; }
  int     getNIntMax()   { return m_nIntMax; }
  double *getWeightInt() { return m_weightInt; }
  double *getEffInt()    { return m_effInt; }
  double *getBkgInt()    { return m_bkgInt; }
  //
  double  getDmus() { return m_dmus; }
  int     getNBelt() { return m_nBelt; }
  int     getNBeltMax() { return m_nBeltMax; }
  bool    isValidBestMu()   { return m_validBestMu; }
  double *getBestMuProb() { return m_bestMuProb; }
  double *getBestMu() { return m_bestMu; }
  double *getMuProb() { return m_muProb; }
  double *getLhRatio() { return m_lhRatio; }
  //
  double getSumProb()    { return m_sumProb; }
  double getLowerLimit() { return m_lowerLimit; }
  double getUpperLimit() { return m_upperLimit; }
  double getLowerLimitNorm() { return m_lowerLimitNorm; }
  double getUpperLimitNorm() { return m_upperLimitNorm; }
  int    getNuppLim()    { return m_nUppLim; }

private:
  void setInt(double & low, double & high, double & step, double scale, double mean, double sigma, DISTYPE dt);

  PDF::Poisson m_poisson;
  PDF::Gauss   m_gauss;
  //
  int    m_verbose;
  //
  double m_stepMin;
  // CL, confidence limit
  double m_cl;
  // True signal - used in coverage studies
  double m_sTrue;
  bool   m_coverage;
  // Number of observed events
  int    m_nObserved;
  // Efficiency, gaussian
  double  m_effMeas;
  double  m_effSigma;
  bool    m_effNoDist;
  DISTYPE m_effDist;
  // Background, gaussian
  double  m_bkgMeas;
  double  m_bkgSigma;
  bool    m_bkgNoDist;
  DISTYPE m_bkgDist;
  // correlation between eff and bkg [-1.0..1.0]
  double  m_beCorr;
  //
  int    m_intNdef;
  // range and steps in double integral (7), in principle infinite
  Range  m_effRangeInt;
  double m_effIntScale;
  // range and steps in double integral (7)
  Range  m_bkgRangeInt;
  double m_bkgIntScale;
  // Test range for the likelihood ratio calculation (4)
  Range  m_hypTest;
  double m_hypTestMin;
  double m_hypTestMax;
  double m_hypTestStep;
  int    m_hypTestN; // obtained from range and stepsize
  //
  // Kernel of double integral (7)
  //
  bool    m_validInt;  // true if integral is valid
  int     m_nInt;      // == m_effIntN*m_bkgIntN
  int     m_nIntMax;   // allocated
  double *m_weightInt; // array containing weights (Gauss(e)*Gauss(b)*de*db), size = n_points
  double *m_effInt;    // array containg e used in integral (e)
  double *m_bkgInt;    // array containg b used in integral (b)
  //
  // Arrays of best fit and limits
  //
  // POLE
  double  m_dmus;       // step size in search for s_best (LHR)
  int     m_nBelt;      // how many Nobs are tested to find R (likelihood ratio)
  int     m_nBeltMax;   // dito allocated
  bool    m_suggestBelt;// if true, always call suggestBelt(); set to true if setBelt(v) is called with v<1
  //  std::vector<int> m_nBeltList; // list of suggested nBelt - filled in constructor
  bool    m_validBestMu;//
  double *m_bestMuProb; // prob. of best mu=e*s+b, index == (N observed)
  double *m_bestMu;     // best mu=e*s+b
  double *m_muProb;     // prob for mu
  double *m_lhRatio;    // likelihood ratio
  double m_sumProb;    // sum of probs for conf.belt construction
  bool   m_foundLower; // true if lower limit is found
  bool   m_foundUpper; // true if an upper limit is found
  double m_lowerLimit; // lowerlimit obtained from ordering (4)
  double m_upperLimit; // upperlimit
  double m_lowerLimitNorm; // sum of the probabilities in the belt given by the lower limit
  double m_upperLimitNorm; // dito upper limit (useful to check whether nbelt was enough or not)
  double m_maxNorm;  // max probability sum (should be very near 1)
  double m_normMaxDiff; // max(norm-1.0) allowed before giving a warning
  int    m_nUppLim;  // number of points to scan the hypothesis after the first upper limit is found
  //
  bool   m_useNLR; // Use Gary Hills likelihood ratio
};

inline bool Pole::normOK(double p) {
  return (fabs(p-1.0)<m_normMaxDiff);
}

inline double Pole::calcProb(int n, double s) {  
  double p = 0.0; 
  //
  for(int i=0;i<m_nInt;i++) {
    p += m_weightInt[i]*m_poisson.getVal(n,m_effInt[i]*s + m_bkgInt[i]);
  }
  return p;
}
//
// allow for gcc 2.96 fix
//
# if __GNUC__ > 2
inline void coutFixed(int precision, double val) {
  std::cout << std::fixed << std::setprecision(precision) << val;
}
inline void coutFixed(int precision, int val) {
  std::cout << std::fixed << std::setprecision(precision) << val;
}
# else
void coutFixed(int precision, double val) {
  static char fmt[32];
  sprintf(&fmt[0],"%%.%df",precision);
  printf(fmt,val);
}

void coutFixed(int precision, int val) {
  static char fmt[32];
  sprintf(&fmt[0],"%%%dd",precision);
  printf(fmt,val);
}
#endif

#endif

