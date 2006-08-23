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
#include <iomanip>
#include <vector>
#include <cmath>
#include <ctime>
#include <utility>
#include <map>
#include <string>
#include <list>
#include <fstream>
#include <sstream>

#if __GNUC__ < 3
#include <cstdio>
#endif

#include "Tools.h"
#include "Range.h"
#include "BeltEstimator.h"
#include "Measurement.h"



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
 *  - setNObserved() : Number of observed events
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
 *  - calcBelt() : Calculates the confidence belt [n1(s,b),n2(s,b)] for all (s,b).
 *  - calcBelt(s,n1,n2,v) : Dito but for a specific (s,b)
 *
 *  Finding \f$s_{best}\f$
 *  - setBestMuStep() : Sets the precision in findBestMu().\n
 *    Default = 0.01 and it should normally be fine.
 *  - setMethod() : sets the Likelihood ratio method (MBT or FHC2)
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

enum RLMETHOD {
  RL_NONE=0,
  RL_FHC2,
  RL_MBT
};

class Pole {
public:
  Pole();
 ~Pole();

  //! running with the current setup
  void execute();
  void exeEvent(bool first);
  void exeFromFile();

  //! Main parameters
  void setCL(double cl)    { m_cl = cl; if ((cl>1.0)||(cl<0.0)) m_cl=0.9;}

  void setMethod( RLMETHOD m ) { m_method = m; }
  void setMethod( int m ) { m_method = RLMETHOD(m); }

  //! Set input file
  void setInputFile( const char *s ) { m_inputFile = s; }
  void setInputFileLines( int nmax ) { m_inputFileLines = nmax; }

  //! Set measurement
  void setMeasurement( const MeasPoisEB & m ) { m_measurement.copy(m); }
  //
  void setNObserved(int nobs) { m_nBeltUsed = nobs; m_measurement.setObsVal(nobs); }
  //  {m_measurement.setNObserved(nobs); }
  //! distribution info on eff and bkg
  void setEffPdf(double mean,double sigma, PDF::DISTYPE dist=PDF::DIST_GAUS) {
    m_measurement.setEffPdf(mean,sigma,dist);
    m_validBestMu = false;
  }
  void setBkgPdf(double mean,double sigma, PDF::DISTYPE dist=PDF::DIST_GAUS) {
    m_measurement.setBkgPdf(mean,sigma,dist);
    m_validBestMu = false;
  }
  // set observed eff and bkg, not changing the PDF def
  void setEffObs(double mean) {
    m_measurement.setEffObs(mean);
    m_validBestMu = false;
  }
  void setBkgObs(double mean) {
    m_measurement.setBkgObs(mean);
    m_validBestMu = false;
  }
  // set ditto using the pdf mean
  void setEffObs() {
    m_measurement.setEffObs();
    m_validBestMu = false;
  }
  void setBkgObs() {
    m_measurement.setBkgObs();
    m_validBestMu = false;
  }
  // generate a random observation (observable + nuisance parameters)
  void generatePseudoExperiment() {
    m_measurement.generatePseudoExperiment();
    m_validBestMu = false;
  }
  // eff,bkg correlation...
  void setEffPdfBkgCorr(double corr)    { m_measurement.setBEcorr(corr); }
  ////////////////////////////////
  //
  bool checkEffBkgDists();
  bool isFullyCorrelated() { return false; } //m_measurement.isFullyCorrelated(); } // { return (((fabs(fabs(m_beCorr)-1.0)) < 1e-16)); }
  bool isNotCorrelated()   { return true;  } // m_measurement.isNotCorrelated(); }   // { return (fabs(m_beCorr) < 1e-16); }
 
  //////////  bool checkParams();

  // Set belt max value
  void setBelt(int v)    { m_nBeltMaxUsed = 0; m_nBeltMinUsed = v; m_nBelt = v; }
  int  suggestBelt();                // will suggest a m_nBelt based on no. observed

  // scan parameters for the search for s_best
  void setBestMuStep(double dmus,double stepmin=0.001) { m_bestMuStep = (dmus > 0.0 ? dmus:stepmin); }
  void setBestMuNmax(int n)                            { m_bestMuNmax = n; }

  // set minimum probability considered - TODO: need to check this against precision in hypothesis testing
  void setMinMuProb(double m=-1) {
    if (m<0.0) {
      double e=floor(log10(1-m_cl))-2.0;
      m_minMuProb = pow(10.0,e);
    } else {
      m_minMuProb = m;
    }
  }
  // POLE test hypothesis range
  void setLimitHypStep(double step=0.001) { m_limitHypStep=(step>0 ? step:0.001); } // set the hypothesis step size when scanning for limits
  void setTestHyp(double step=-1.0); // set test range based on the input measurement
  void setTestHyp(double low, double high, double step); // test mu=s+b for belt construction etc (NOT for limit scan)
  void setNuppLim(int n=-1) { m_nUppLim = n; }

  // POLE true signal, used only if coverage run
  void setTrueSignal(double s) { m_sTrue = s; } // true signal
  void setCoverage(bool flag) {m_coverage = flag;} // true if coverage run

  // Debug
  void setVerbose(int v=0) { m_verbose=v; }

  // Tabulated PDF's
  void setPoisson(  const PDF::PoisTab *pdf) {m_poisson=pdf;}
  void setGauss(    const PDF::Gauss   *pdf) {m_gauss=pdf;}
  void setGauss2D(  const PDF::Gauss2D *pdf) {m_gauss2d=pdf;}
  void setLogNormal(const PDF::LogNormal   *pdf) {m_logNorm=pdf;}
  ///////////////////////////////
  //
  void initIntArrays();   // will initialise integral arrays (if needed)
  void initBeltArrays();  // will initialise belt arrays (if needed)
  void initIntegral();    // calculates double integral kernal (eff*bkg*db*de) according to setup (7)
  void initIntegral(std::vector<double> & eff, std::vector<double> & bkg, std::vector<double> & weight);

  //
  double getObservedSignal() { return m_measurement.getSignal(); }
  // POLE
  void findBestMu(int n); // finds the best fit (mu=s+b) for a given n. Fills m_bestMu[n] and m_bestMuProb[n].
  void findAllBestMu();   // dito for all n (loop n=0; n<m_nMuUsed)
  void calcConstruct(double s, bool verb);
  double calcBelt(double s, int & n1, int & n2,bool verb);//,double muMinProb=1e-5); // calculate (4) and find confidence belt
  int  calcLimit(double s); // calculate (4) and find limits
  bool calcLimit(double s, bool scanDown); // calculate (4) and find limits
  void   calcLh(double s); // fills the likelihood array
  double calcLhRatio(double s, int & nb1, int & nb2);//, double minMuProb=1e-6); // fills the likelihood ratio array
  bool limitsOK(); // check if calculated limit is OK using the sum of probs.
  inline const bool normOK(double p) const;
  void setNormMaxDiff(double dpmax=0.001) { m_normMaxDiff=dpmax; }
  void calcPower();
  void calcConstruct();
  void calcNMin();
  void calcBelt();
  bool calcLimits();        // finds CL limits
  bool calcCoverageLimits();//  same as previous but stop if it's obvious that initial true mean is inside or outside
  void setPrintLimitStyle( int style ) { m_printLimStyle = style; }
  void printLimit(bool doTitle=false);
  // 
  void printSetup();
  //
  void printFailureMsg();
  //
  // POLE
  void initAnalysis(); // inits the vectors and calculates the double integral weights
  bool analyseExperiment();  // finds the limits/coverage
  //
  // Access functions
  //
  const bool getMethod() const { return m_method; }
  const bool usesMBT()   const { return (m_method==RL_MBT); }
  const bool usesFHC2()  const { return (m_method==RL_FHC2); }

  const int    getVerbose() const    { return m_verbose; }
  const double getCL() const         { return m_cl; }
  const double getSTrue() const      { return m_sTrue; }
  const bool   getCoverage() const   { return m_coverage; }
  //
  const MeasPoisEB & getMeasurement() const  { return m_measurement; }
  MeasPoisEB & getMeasurement()              { return m_measurement; }
  const int    getNObserved() const          { return m_measurement.getObsVal(); }
  // Efficiency
  const double  getEffObs()   const          { return m_measurement.getEffObs(); }
  const double  getEffPdfMean()  const       { return m_measurement.getEffPdfMean(); }
  const double  getEffPdfSigma() const       { return m_measurement.getEffPdfSigma(); }
  const PDF::DISTYPE getEffPdfDist()  const  { return m_measurement.getEffPdfDist(); }
  // Background
  const double  getBkgObs()  const  { return m_measurement.getBkgObs(); }
  const double  getBkgPdfMean() const  { return m_measurement.getBkgPdfMean(); }
  const double  getBkgPdfSigma() const  { return m_measurement.getBkgPdfSigma(); }
  const PDF::DISTYPE getBkgPdfDist()  const  { return m_measurement.getBkgPdfDist(); }
  const double  getEffPdfBkgCorr() const { return m_measurement.getBEcorr(); }
  // Indep. variable
//   const double  getSVar()     const  { return BeltEstimator::getT(m_measurement.getNObserved(),
// 								  m_measurement.getEffObs(),
// 								  m_measurement.getEffPdfSigma(),
// 								  m_measurement.getBkgObs(),
// 								  m_measurement.getBkgPdfSigma(),
// 								  m_measurement.getNuisanceIntNorm());}

  const double  getSVar()     const  { std::cout << "BELTESTIMATOR not functional!" << std::endl; return 0.0; }

  // Test range for the likelihood ratio calculation (4)
  const Range<double>  *getHypTest() const     { return &m_hypTest; }
  //
  double getEffIntMin() const { return m_measurement.getEff()->getIntXmin(); }
  double getEffIntMax() const { return m_measurement.getEff()->getIntXmax(); }
  double getEffIntN()   const { return m_measurement.getEff()->getIntN(); }
  double getBkgIntMin() const { return m_measurement.getBkg()->getIntXmin(); }
  double getBkgIntMax() const { return m_measurement.getBkg()->getIntXmax(); }
  double getBkgIntN()   const { return m_measurement.getBkg()->getIntN(); }
  double getEffIntNorm() const { return m_measurement.getEff()->getIntegral(); }
  double getBkgIntNorm() const { return m_measurement.getBkg()->getIntegral(); }
  double getIntNorm() const { return m_measurement.getNuisanceIntNorm(); }
  //
  const double  getLsbest(int n) const;
  const int     getNBelt() const { return m_nBelt; }
  const int     getNBeltUsed() const { return m_nBeltUsed; }
  const int     getNBeltMinUsed() const { return m_nBeltMinUsed; }
  const int     getNBeltMaxUsed() const { return m_nBeltMaxUsed; }
  const bool    isValidBestMu() const  { return m_validBestMu; }
  //
  const double  getBestMuStep() const { return m_bestMuStep; }
  const int     getBestMuNmax() const { return m_bestMuNmax; }
  const std::vector<double> & getBestMuProb() const { return m_bestMuProb; }
  const std::vector<double> & getBestMu() const { return m_bestMu; }
  const std::vector<double> & getMuProb() const { return m_muProb; }
  const std::vector<double> & getLhRatio() const { return m_lhRatio; }
  const double getMinMuProb() const { return m_minMuProb; }
  const double getMuProb(int n) const { if ((n>m_nBeltMaxUsed)||(n<m_nBeltMinUsed)) return 0.0; return m_muProb[n];}
  //
  const double getLimitHypStep() const { return m_limitHypStep; }
  const double getSumProb() const    { return m_sumProb; }
  const double getLowerLimit() const { return m_lowerLimit; }
  const double getUpperLimit() const { return m_upperLimit; }
  const double getLowerLimitNorm() const { return m_lowerLimitNorm; }
  const double getUpperLimitNorm() const { return m_upperLimitNorm; }
  const int    getNuppLim() const    { return m_nUppLim; }
  const double getRejS0P()  const { return m_rejs0P; }
  const int    getRejS0N1() const { return m_rejs0N1; }
  const int    getRejS0N2() const { return m_rejs0N2; }

  const std::string & getInputFile() const { return m_inputFile; }
  const int getInputFileLines() const { return m_inputFileLines; }

private:

  const PDF::PoisTab   *m_poisson;
  const PDF::Gauss     *m_gauss;
  const PDF::Gauss2D   *m_gauss2d;
  const PDF::LogNormal *m_logNorm;
  //
  int    m_verbose;
  int    m_printLimStyle;
  //
  // CL, confidence limit
  double m_cl;

  // RL method
  RLMETHOD m_method;

  // True signal - used in coverage studies
  double m_sTrue;
  bool   m_coverage;
  // Measurement
  MeasPoisEB m_measurement;
  // correlation between eff and bkg [-1.0..1.0]
  double  m_beCorr;
  ////////////////////////////////////////////////////
  //
  // Test range for belt construction (calcBelt()), power calculation (calcPower()) and construction (calcConstruct())
  Range<double>  m_hypTest;
  //
  // Arrays of best fit and limits
  //
  // POLE
  int     m_nBelt;        // beltsize used for explicit construction of conf. belt (calcBelt() etc)
  int     m_nBeltUsed;    // dynamic beltsize determined in calcLhRatio() - default = max(2,N(obs))
  int     m_nBeltMinUsed; // the minimum nBelt used for the calculation
  int     m_nBeltMaxUsed; // the maximum nBelt used for the calculation
  int     m_nBeltMinLast; // the minimum nBelt used in previous call to calcLhRatio()
  //  std::vector<int> m_nBeltList; // list of suggested nBelt - filled in constructor
  //
  bool    m_validBestMu;//
  double  m_bestMuStep;       // step size in search for s_best (LHR)
  int     m_bestMuNmax;   // maximum N in search for s_best (will locally nodify dmus)
  std::vector<double> m_bestMuProb; // prob. of best mu=e*s+b, index == (N observed)
  std::vector<double> m_bestMu;     // best mu=e*s+b
  std::vector<double> m_muProb;     // prob for mu
  std::vector<double> m_lhRatio;    // likelihood ratio
  double m_minMuProb;  // minimum probability accepted
  //
  double m_limitHypStep; // step size when searching for the limits (see calcLimits())
  double m_sumProb;    // sum of probs for conf.belt construction - set by calcLimit()
  double m_prevSumProb;// ditto from previous scan - idem
  double m_scanBeltNorm; // sum(p) for all n used in belt at the current s - idem
  bool   m_lowerLimitFound; // true if lower limit is found
  bool   m_upperLimitFound; // true if an upper limit is found
  double m_lowerLimit; // lowerlimit obtained from ordering (4)
  double m_upperLimit; // upperlimit
  double m_lowerLimitNorm; // sum of the probabilities in the belt given by the lower limit
  double m_upperLimitNorm; // dito upper limit (useful to check whether nbelt was enough or not)
  //
  double m_maxNorm;  // max probability sum (should be very near 1)
  double m_normMaxDiff; // max(norm-1.0) allowed before giving a warning
  int    m_nUppLim;  // number of points to scan the hypothesis after the first upper limit is found REMOVE!!!!
  double m_rejs0P;  // probability for belt at s=0
  int    m_rejs0N1; // min N at s=0
  int    m_rejs0N2; // ditto max N

  std::string m_inputFile; // input file with data
  int         m_inputFileLines; // number of lines to read; if < 1 => read ALL lines
  //


  // TO BE REMOVED
  double calcLimitOLD(double s); // calculate (4) and find limits, returns probability for given signal hypothesis
  bool calcLimitsOLD();        // finds CL limits
  bool calcCoverageLimitsOLD();
};

inline const double Pole::getLsbest(int n) const {
  double rval = 0.0;
  if (usesMBT()) {
    double g;
    if (n>m_measurement.getBkgObs()) {
      g = static_cast<double>(n);
    } else {
      g = m_measurement.getBkgObs();
    }
    rval = m_poisson->getVal(n,g);
  } else {
    rval = m_bestMuProb[n];
  }
  return rval;
}

inline const bool Pole::normOK(double p) const {
  return (fabs(p-1.0)<m_normMaxDiff);
}

#endif

