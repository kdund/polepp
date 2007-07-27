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
#include "Measurement.h"
#include "Integrator.h"


/*! @class Pole
 *
 *  @brief The main class containing the methods for calculating limits.
 *  
 *  \b SETUP
 *
 *  Basic
 *  - setMethod() : Sets the Likelihood ratio method (MBT=2 or FHC2=1)
 *  - setCL() : Confidence limit, default 0.9
 *    the requested confidence [0.0,1.0]
 *  - setNObserved() : Number of observed events
 *  - setEffPdf() : Measured efficiency\n
 *    efficiency distribution (mean,sigma and distribution type ( PDF::DISTYPE ))
 *  - setBkgPdf() : Measured background\n
 *    background distribution (mean,sigma and distribution type ( PDF::DISTYPE ))
 *  - setEffBkgCorr() : Correlation between eff and bkg (if applicable)
 *    correlation coefficient [-1.0,1.0]
 *  - setEffObs( double m ) : Set the observed efficiency\n
 *    if no argument is given, the mean of the pdf is used
 *  - setBkgObs( double m ) : Set the observed background\n
 *    see previous
 *
 *  Integration
 *  - setIntEffNSigma() : Integration range of efficiency\n
 *    The integration range must cover the PDF such that the tails are
 *    negligable.
 *  - setIntBkgNSigma() : Dito, background
 *
 *  Belt construction
 *  - calcBelt() : Calculates the confidence belt [n1(s,b),n2(s,b)]
 *
 *  Finding \f$s_{best}\f$
 *  - setBestMuStep() : Sets the precision in findBestMu().\n
 *    Default = 0.01 and it should normally be fine.
 *  - setMethod() : sets 
 *
 *  Hypothesis testing
 *  - setHypTestRange() : test hypothesis range, $s_{best}$\n
 *    This sets the range and step size (precision) when calculating the construction or confidence belt.\n
 *    For limit calculations, this setting has no effect. It is now dynamically set.
 *
 *  Coverage related
 *  - setTrueSignal() : Set the true signal.
 *  - setCoverage() : Set the coverage flag.\n
 *    Setting this flag causes the limit calculation to terminate the scan as soon as it is
 *    decided whether the true signal is inside or outside the confidence limit.
 *
 *  PDFs
 *  - setPoisson() : Initialises the poisson generator - use PDF::gPoisTab.\n
 *    Use the tabulated version, PDF::gPoisTab.
 *  - setGauss() : Ditto but for a gauss pdf.\n
 *    Use PDF::gGauss
 *  - setGauss2D() : Ditto but for a 2d gauss pdf.\n
 *    Use PDF::gGauss2D
 *  - setLogNormal() : Ditto but for a log normal pdf.\n
 *    Use PDF::gLogNormal
 *
 *  See the files argsPole.cxx and argsCoverage.cxx for examples of usage.
 *
 *  \b RUNNING
 *
 *  Print
 *  - printSetup() : Prints the setup to stdout.
 *  - printLimit() : Prints the calculated limit.
 *
 *  General
 *  - execute() : Main routine to call. Initialises and runs with the current setup.
 *  - analyseExperiment() : Calculates the limit using the current setup.
 *
 *  Debug
 *  - setVerbose() : Sets verbose level.
 *
 *  \b TOOLS
 *  - polelim.cxx : limit calculation
 *  - polecov.cxx : coverage calculator
 *  - polebelt.cxx : confidence belt calculator
 *  - poleconst.cxx : prints the full construction
 *    This code prints the likelihood ratios obtained in the $s_{hyp}$:$N_{obs}$ plane.
 *
 * @author Fredrik Tegenfeldt (fredrik.tegenfeldt@cern.ch)
 */
//  e^{\frac{-(b-b')^2}{2\sigma_b^2}}\;\;
//  e^{\frac{- (1 -\epsilon')^2}{2\sigma_{\epsilon}^2}}

namespace PDF {
  class Base;
};

namespace LIMITS {

  class Pole;

  struct PoleData {
    const Pole      *polePtr;
    const PDF::Base *pdfEff;
    const PDF::Base *pdfBkg;
    const PDF::Base *pdfObs;
    int              effIndex;
    int              bkgIndex;
    int              nobs;
    double           signal;
    double           effObs;
    double           deffObs;
    double           bkgObs;
    double           dbkgObs;
  };

  class PoleIntegrator {
  public:
    inline PoleIntegrator();
    inline ~PoleIntegrator();

    inline void setPole( const Pole *pole );

    inline void setParameters( std::vector<double> & pars );

    inline const Integrator *getIntegrator() const;
    inline Integrator       *integrator();
    //
    inline void   go();
    inline double result() const;

    inline int    getEffIndex()  const;
    inline double getEffIntMin() const;
    inline double getEffIntMax() const;
    inline int    getBkgIndex()  const;
    inline double getBkgIntMin() const;
    inline double getBkgIntMax() const; 
  private:
    struct PoleData m_poleData;
    IntegratorVegas m_integrator;
  };


  enum RLMETHOD {
    RL_NONE=0,
    RL_FHC2,
    RL_MBT
  };

  class Pole {
  public:
    /*! @name Constructors/destructor/initialisor */
    //@{
    //! main constructor
    Pole() { initDefault(); }

    //! destructor
    ~Pole() {}

    //! initialise to default
    void initDefault();
    //@}

    /*! @name Main interfaces */
    //@{
    //! do the limit calculation or coverage determination
    bool analyseExperiment();

    //! generate a random observation (observable + nuisance parameters)
    void generatePseudoExperiment() {
      m_measurement.generatePseudoExperiment();
      m_validBestMu = false;
    }

    //! running with the current setup
    void execute();

    //! execute one event
    void exeEvent(bool first);

    //! execute events from an input file
    void exeFromFile();
    //@}


    /*! @name Set parameters concerning the measurement */

    //! Set input file
    void setInputFile( const char *s ) { m_inputFile = s; }
    //! Set the number of lines to read from the input file
    void setInputFileLines( int nmax ) { m_inputFileLines = nmax; }

    //! Set measurement
    void setMeasurement( const MEAS::MeasPoisEB & m ) { m_measurement.copy(m); }

    //! set the confidence level
    void setCL(double cl)    { m_cl = cl; if ((cl>1.0)||(cl<0.0)) m_cl=0.9;}

    //! set the CL method to be used
    void setMethod( RLMETHOD m ) { m_method = m; }
    //! set the Cl method to be used (int input)
    void setMethod( int m ) { m_method = RLMETHOD(m); }

    //! set the number of observed events
    void setNObserved(int nobs) { m_nBeltUsed = nobs; m_measurement.setObsVal(nobs); }

    //! set pdf of efficiency
    void setEffPdf(double mean,double sigma, PDF::DISTYPE dist=PDF::DIST_GAUS) {
      m_measurement.setEffPdf(mean,sigma,dist);
      m_validBestMu = false;
    }
    //! set pdf of background
    void setBkgPdf(double mean,double sigma, PDF::DISTYPE dist=PDF::DIST_GAUS) {
      m_measurement.setBkgPdf(mean,sigma,dist);
      m_validBestMu = false;
    }
    //! set pdf mean of efficiency
    void setEffPdfScale( double s=1.0 ) {
      m_measurement.setEffScale(s);
      m_validBestMu = false;
    }
    //! set pdf mean of efficiency
    void setEffPdfMean( double m ) {
      m_measurement.setEffPdfMean(m);
      m_validBestMu = false;
    }
    //! set pdf sigma of efficiency
    void setEffPdfSigma( double s ) {
      m_measurement.setEffPdfSigma(s);
      m_validBestMu = false;
    }
    //! set pdf mean of background
    void setBkgPdfScale( double s=1.0 ) {
      m_measurement.setBkgScale(s);
      m_validBestMu = false;
    }
    //! set pdf mean of background
    void setBkgPdfMean( double m ) {
      m_measurement.setBkgPdfMean(m);
      m_validBestMu = false;
    }
    //! set pdf sigma of background
    void setBkgPdfSigma( double s ) {
      m_measurement.setBkgPdfSigma(s);
      m_validBestMu = false;
    }
    //! set the observed efficiency
    void setEffObs(double mean) {
      m_measurement.setEffObs(mean);
      m_validBestMu = false;
    }
    //! set the observed background
    void setBkgObs(double mean) {
      m_measurement.setBkgObs(mean);
      m_validBestMu = false;
    }
    //! set the observed efficiency using the pdf mean
    void setEffObs() {
      m_measurement.setEffObs();
      m_validBestMu = false;
    }
    //! set the observed background using the pdf mean
    void setBkgObs() {
      m_measurement.setBkgObs();
      m_validBestMu = false;
    }
    //! set eff,bkg correlation...
    void setEffBkgPdfCorr(double corr)    { m_measurement.setBEcorr(corr); }

    //! set the integral range for efficiency
    /*!
      integral range is calculated in initIntegral()
      \param nsigma is the number of sigmas to cover
    */
    void setIntEffNSigma( double nsigma ) { m_effIntNSigma = nsigma; }
    //! set the integral range for background
    /*!
      integral range is calculated in initIntegral()
      \param nsigma is the number of sigmas to cover
    */
    void setIntBkgNSigma( double nsigma ) { m_bkgIntNSigma = nsigma; }

    //! set the number of calls used by GSL integrator
    void setIntGslNCalls(int n) { m_gslIntNCalls = n; }

    //! set range in signal for tabulated integral
    void setIntSigRange( double smin, double smax, int nsteps ) { m_intTabSRange.setRange( smin, smax, 0.0, nsteps ); }

    //! set range in N(obs) for tabulated integral
    void setIntNobsRange( int nmin, int nmax )  { m_intTabNRange.setRange( nmin, nmax, 1 ); }

    //! set tabulation flag
    void setTabulateIntegral( bool f ) { m_tabulateIntegral = f; }

    //! set a scaling number to scale the output limits
    void setScaleLimit(double s) { m_scaleLimit = (s>0.0 ? s:1.0); }
    //@}

    /*! @name Set parameters concerning the precision of the calculations */
    //@{
    //! scan step size for finding s_best
    void setBestMuStep(double dmus,double stepmin=0.001) { m_bestMuStep = (dmus > 0.0 ? dmus:stepmin); }
    //! maximum number of points allowed searching for s_best
    void setBestMuNmax(int n)                            { m_bestMuNmax = n; }

    //! set the cutoff probability for the tails in calcLhRatio()
    /*!
      If the given limit is < 0, then it's calculated as floor(log10(1-CL))-2.0.
    */
    // TODO: need to check this against precision in hypothesis testing
    void setMinMuProb(double m=-1) {
      if (m<0.0) {
        double e=floor(log10(1-m_cl))-2.0;
        m_minMuProb = pow(10.0,e);
      } else {
        m_minMuProb = m;
      }
    }
    //! binary search stopping threshold
    /*!
      The binary search in calcLimit() will stop when the change in mutest is below this value.
      Default value is 0.001.
    */
    void setBSThreshold(double step=0.001)    { m_thresholdBS    = (step>0 ? step:0.001); }
    //! set the threshold defining when a CL is close enough
    /*!
      This value defines when the probability for a signal is close enough to the CL.
      The cutoff is $x = 1.0 - (\alpha_{est}/\alpha_{req})$.
      Default value is 0.01.
    */
    void setPrecThreshold( double da=0.01)   { m_thresholdPrec = (da>0 ? da:0.01); }

    //! set range for mutest in calcBelt() etc (NOT used in the limit calculation)
    void setHypTestRange(double low, double high, double step);

    //! set maximum difference from unity in normOK()
    void setNormMaxDiff(double dpmax=0.001) { m_normMaxDiff=dpmax; }
    //@}


    /*! @name Coverage specific setup */
    //@{
    //! set the true signal
    void setTrueSignal(double s)   { m_measurement.setTrueSignal(s); }
    //! set the 'use-coverage' flag
    /*!
      If this flag is set true, the calcLimitCoverage() is called instead of calcLimit().
    */
    void setUseCoverage(bool flag) {m_coverage = flag;}
    //@}


    /*! @name Consistency checks - THESE NEED WORK - NOT YET USED */
    //@{
    //! check that the distributions are ok
    bool checkEffBkgDists();
    //! check that all parameters are ok
    bool checkParams() { return true; }
    //@}

    /*! @name Debugging */
    //@{
    //! set verbosity - 0 is no verbosity
    void setVerbose(int v=0) { m_verbose=v; }
    //@}

    /*! @name PDF instances */
    //@{
    void setPoisson(  const PDF::Poisson   *pdf)  {m_poisson=pdf;}
    void setGauss(    const PDF::Gauss     *pdf)  {m_gauss=pdf;}
    void setGauss2D(  const PDF::Gauss2D   *pdf)  {m_gauss2d=pdf;}
    void setLogNormal(const PDF::LogNormal *pdf)  {m_logNorm=pdf;}
    //@}

    /*! @name Initialising the calculation */
    //@{
    //! initialises integrals etc, calls initIntArrays() among others
    void initAnalysis();
    //! will initialise belt arrays
    void initBeltArrays();
    //! init pole integrator
    void initIntegral();
    //! init tabulated pole integral
    void initTabIntegral();
    //@}


    /*! @name Calculation of belt, limits, likelihood ratios etc */
    //@{
    //! calculate probability P(N(obs) | signal) using table
    inline double calcProb( int n, double s );
    //! finds all mu_best - only used if the method is RLMETHOD::RL_FHC2
    void findAllBestMu();   // dito for all n (loop n=0; n<m_nMuUsed)
    //! finds the mu_best for the given N
    void findBestMu(int n);
    //! calculate the construct
    void calcConstruct();
    //! calculate the construct for the given signal
    void calcConstruct(double s, bool title);
    //! calculate the confidence belt
    void calcBelt();
    //! calculate the confidence belt for the given signal
    double calcBelt(double s, int & n1, int & n2,bool verb, bool title);
    //! scan for lower limit
    bool scanLowerLimit( double mustart, double p0 );
    //! scan for upper limit
    bool scanUpperLimit( double mustart, double p0 );
    //! calculate the confidence limits
    bool calcLimit();
    //! calculate the confidence limit probability for the given signal
    int  calcLimit(double s, double & prec);
    //! calculate the confidence limit probability for the given signal
    int  calcLimit(double s) { double prec; return calcLimit(s,prec); }
    //! check if s(true) lies within the limit
    bool calcCoverageLimit();
    //! calculate the likelihood ratio
    double calcLhRatio(double s, int & nb1, int & nb2);
    //! calculate the power
    void calcPower();
    //! calculate minimum N(obs) that will reject s=0
    void calcNMin();
    //! reset all variables pertaining to the limit calculation
    void resetCalcLimit();
    //! checks if the limits are OK
    bool limitsOK();
    //@}

    /*! @name Output */
    //@{
    //! set the print style - not very elaborate at the moment - can be 0 or 1
    void setPrintLimitStyle( int style ) { m_printLimStyle = style; }
    //! print the result of the limit calculation
    void printLimit(bool doTitle=false);
    //! print setup
    void printSetup();
    //! print failure message
    void printFailureMsg();
    //@}

    /*! @name Accessor functions */
    //@{
    const double getCL()                  const { return m_cl; }
    const double getObservedSignal()      const { return m_measurement.getSignal(); }
    const double getTrueSignal()          const { return m_measurement.getTrueSignal(); }
    const bool   getMethod()              const { return m_method; }
    const bool   usesMBT()                const { return (m_method==RL_MBT); }
    const bool   usesFHC2()               const { return (m_method==RL_FHC2); }

    const MEAS::MeasPoisEB & getMeasurement()   const { return m_measurement; }
    MEAS::MeasPoisEB &       getMeasurement()         { return m_measurement; }
    const int          getNObserved()     const { return m_measurement.getObsVal(); }

    const PDF::Base   *getObsPdf()        const { return m_measurement.getObservable()->getPdf(); }

    const PDF::Base   *getEffPdf()        const { return m_measurement.getEff()->getPdf(); }
    const double       getEffObs()        const { return m_measurement.getEffObs(); }
    const double       getEffPdfMean()    const { return m_measurement.getEffPdfMean(); }
    const double       getEffPdfSigma()   const { return m_measurement.getEffPdfSigma(); }
    const PDF::DISTYPE getEffPdfDist()    const { return m_measurement.getEffPdfDist(); }
    const double       getEffScale()      const { return m_measurement.getEffScale(); }

    const PDF::Base   *getBkgPdf()        const { return m_measurement.getBkg()->getPdf(); }
    const double       getBkgObs()        const { return m_measurement.getBkgObs(); }
    const double       getBkgPdfMean()    const { return m_measurement.getBkgPdfMean(); }
    const double       getBkgPdfSigma()   const { return m_measurement.getBkgPdfSigma(); }
    const PDF::DISTYPE getBkgPdfDist()    const { return m_measurement.getBkgPdfDist(); }
    const double       getBkgScale()      const { return m_measurement.getBkgScale(); }

    const double       getEffPdfBkgCorr() const { return m_measurement.getBEcorr(); }

    const bool         useCoverage()      const { return m_coverage; }
    const bool         truthCovered()     const { return m_coversTruth; }

    const int            getVerbose()       const { return m_verbose; }

    const Range<double> *getHypTest()       const { return &m_hypTest; }

    const double getEffIntMin() const { return m_poleIntegrator.getEffIntMin(); }
    const double getEffIntMax() const { return m_poleIntegrator.getEffIntMax(); }
    const double getBkgIntMin() const { return m_poleIntegrator.getBkgIntMin(); }
    const double getBkgIntMax() const { return m_poleIntegrator.getBkgIntMax(); }
    //    const double getEffIntNorm() const { return m_measurement.getEff()->getIntegral(); }
    //    const double getBkgIntNorm() const { return m_measurement.getBkg()->getIntegral(); }
    //    const double getIntNorm() const { return m_poleIntTable.get(); }
    //
    const Range<double> *getIntSigRange()  const { return &m_intTabSRange; }
    const Range<int>    *getIntNobsRange() const { return &m_intTabNRange; }

    const double  getSbest(int n) const;
    const double  getLsbest(int n) const;
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
    const double getBSThreshold() const { return m_thresholdBS; }
    const double getPrecThreshold() const { return m_thresholdPrec; }
    const double getSumProb() const    { return m_sumProb; }
    const double getLowerLimit() const { return m_lowerLimit; }
    const double getUpperLimit() const { return m_upperLimit; }
    const double getLowerLimitNorm() const { return m_lowerLimitNorm; }
    const double getUpperLimitNorm() const { return m_upperLimitNorm; }
    const double getLowerLimitPrec() const { return m_lowerLimitPrec; }
    const double getUpperLimitPrec() const { return m_upperLimitPrec; }
    const double getRejS0P()  const { return m_rejs0P; }
    const int    getRejS0N1() const { return m_rejs0N1; }
    const int    getRejS0N2() const { return m_rejs0N2; }

    const std::string & getInputFile() const { return m_inputFile; }
    const int getInputFileLines() const { return m_inputFileLines; }

    //@}

    /*! @name Some tools */
    //@{
    const bool isFullyCorrelated() const { return false; }
    const bool isNotCorrelated()   const { return true;  }
    const bool normOK(double p)    const { return (fabs(p-1.0)<m_normMaxDiff); }
    //@}

    static const int      s_intEffInd = 0;
    static const int      s_intBkgInd = 1;
    static const int      s_tabSigInd  = 1;
    static const int      s_tabNobsInd = 0;

  private:

    const PDF::Poisson   *m_poisson;
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
    bool   m_coverage;
    bool   m_coversTruth;
    // Measurement
    MEAS::MeasPoisEB m_measurement;
    // correlation between eff and bkg [-1.0..1.0]
    double  m_beCorr;

    PoleIntegrator            m_poleIntegrator; /**< Pole Integrator wrapper class */
    int                       m_gslIntNCalls;   /**< number of calls used by GSL integrator */
    double                    m_effIntNSigma;   /**< defines the integration range in N(sigmas) */
    double                    m_bkgIntNSigma;   /**< for bkg */

    Tabulator<PoleIntegrator> m_poleIntTable;   /**< the table of poleIntegrator   */
    Range<double>             m_intTabSRange;   /**< tabulated signal range */
    Range<int>                m_intTabNRange;   /**< tabulated N(obs) range */
    bool                      m_tabulateIntegral; /**< flag; if true, tabulate the integral */

    ////////////////////////////////////////////////////
    //
    // Test range for belt construction (calcBelt()), power calculation (calcPower()) and construction (calcConstruct())
    Range<double>  m_hypTest;
    //
    // Arrays of best fit and limits
    //
    // POLE
    std::vector<double> m_calcProbBuf;
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
    double m_thresholdBS; // binary search threshold
    double m_thresholdPrec; // threshold for accepting a CL in calcLimit(s)
    double m_sumProb;    // sum of probs for conf.belt construction - set by calcLimit()
    double m_scanBeltNorm; // sum(p) for all n used in belt at the current s - idem
    bool   m_lowerLimitFound; // true if lower limit is found
    bool   m_upperLimitFound; // true if an upper limit is found
    double m_lowerLimit; // lowerlimit obtained from ordering (4)
    double m_upperLimit; // upperlimit
    double m_lowerLimitNorm; // sum of the probabilities in the belt given by the lower limit
    double m_upperLimitNorm; // dito upper limit (useful to check whether nbelt was enough or not)
    double m_lowerLimitPrec; // precision for lower limit
    double m_upperLimitPrec; // precision for upper limit
    //
    double m_maxNorm;  // max probability sum (should be very near 1)
    double m_normMaxDiff; // max(norm-1.0) allowed before giving a warning
    double m_rejs0P;  // probability for belt at s=0
    int    m_rejs0N1; // min N at s=0
    int    m_rejs0N2; // ditto max N

    double m_scaleLimit; // a factor to scale the resulting limit (in printLimit() )

    std::string m_inputFile; // input file with data
    int         m_inputFileLines; // number of lines to read; if < 1 => read ALL lines
    //
  };

  PoleIntegrator::PoleIntegrator() {
    m_poleData.polePtr = 0;
  }

  PoleIntegrator::~PoleIntegrator() {
  }

  void PoleIntegrator::setPole( const Pole *pole ) {
    m_poleData.polePtr = pole;
    m_poleData.pdfObs  = pole->getObsPdf();
    m_poleData.pdfEff  = pole->getEffPdf();
    m_poleData.pdfBkg  = pole->getBkgPdf();
    m_poleData.effObs  = pole->getEffObs();
    m_poleData.deffObs = pole->getEffPdfSigma();
    m_poleData.bkgObs  = pole->getBkgObs();
    m_poleData.dbkgObs = pole->getBkgPdfSigma();
    bool effConst = PDF::isConstant(pole->getEffPdfDist());
    bool bkgConst = PDF::isConstant(pole->getBkgPdfDist());
    if (effConst) {
      m_poleData.effIndex=-1;
      if (!bkgConst) m_poleData.effIndex=0;
    }
    if (bkgConst) {
      m_poleData.bkgIndex=-1;
      if (!effConst) m_poleData.effIndex=0;
    }
    if ((!effConst) && (!bkgConst)) {
      m_poleData.effIndex=pole->s_intEffInd;
      m_poleData.bkgIndex=pole->s_intBkgInd;
    }
    
    m_integrator.setFunctionParams( &m_poleData );
  }
  void PoleIntegrator::setParameters( std::vector<double> & pars ) {
    // parameter [s_tabNobsInd] = N(obs)
    // parameter [s_tabSigInd]  = signal
    m_poleData.nobs    = static_cast<int>(pars[m_poleData.polePtr->s_tabNobsInd]+0.5);
    m_poleData.signal  = pars[m_poleData.polePtr->s_tabSigInd];
  }

  const Integrator *PoleIntegrator::getIntegrator() const { return & m_integrator; }
  Integrator       *PoleIntegrator::integrator()          { return & m_integrator; }

  int    PoleIntegrator::getEffIndex()  const { return m_poleData.effIndex; }
  double PoleIntegrator::getEffIntMin() const { return (m_poleData.effIndex<0 ? m_poleData.effObs : m_integrator.getIntXmin( m_poleData.effIndex )); }
  double PoleIntegrator::getEffIntMax() const { return (m_poleData.effIndex<0 ? m_poleData.effObs : m_integrator.getIntXmax( m_poleData.effIndex )); }
  int    PoleIntegrator::getBkgIndex()  const { return m_poleData.bkgIndex; }
  double PoleIntegrator::getBkgIntMin() const { return (m_poleData.bkgIndex<0 ? m_poleData.bkgObs : m_integrator.getIntXmin( m_poleData.bkgIndex )); }
  double PoleIntegrator::getBkgIntMax() const { return (m_poleData.bkgIndex<0 ? m_poleData.bkgObs : m_integrator.getIntXmax( m_poleData.bkgIndex )); }
  //
  void   PoleIntegrator::go()           { m_integrator.go(); }
  double PoleIntegrator::result() const { return m_integrator.result(); }
  
  inline const double Pole::getSbest(int n) const {
    double rval = 0.0;
    if (usesMBT()) {
      if (n>m_measurement.getBkgObs()*m_measurement.getBkgScale()) {
        rval = static_cast<double>(n);
      } else {
        rval = m_measurement.getBkgObs()*m_measurement.getBkgScale();
      }
    } else {
      rval = m_bestMu[n];
    }
    return rval;
  }

  /*!
    Calculates the likelihood of the best fit signal to the given number of observations.
    For RLMETHOD::RL_FHC2 it will use the result from findBestMu().
    For RLMETHOD::RL_MBT it will use Poisson(n,g) where g = (n>bkg ? n : bkg ).
  */
  inline const double Pole::getLsbest(int n) const {
    double rval = 0.0;
    if (usesMBT()) {
      double g;
      if (n>m_measurement.getBkgObs()*m_measurement.getBkgScale()) {
        g = static_cast<double>(n);
      } else {
        g = m_measurement.getBkgObs()*m_measurement.getBkgScale();
      }
      rval = m_poisson->getVal(n,g);
    } else {
      rval = m_bestMuProb[n];
    }
    return rval;
  }

};

inline double LIMITS::Pole::calcProb( int n, double s ) {
  m_calcProbBuf[s_tabSigInd] = s;
  m_calcProbBuf[s_tabNobsInd] = static_cast<double>(n);
  return m_poleIntTable.getValue( m_calcProbBuf );
}

// template<>
// inline double Tabulator<LIMITS::PoleIntegrator>::getValue(double n, double s) {
//   size_t tabind;
//   double nmin = m_tabMin[s_tabNobsInd];
//   double nmax = m_tabMax[s_tabNobsInd];
//   double smin = m_tabMin[s_tabSigInd];
//   double smax = m_tabMax[s_tabSigInd];
//   double nval = static_cast<double>(n);
//   if ((n<nmin) || (n>nmax) || (s<smin) || (s>smax)) {
    
//   } else {
    
//   }
//   return 0;
// }

template<>
inline double Tabulator<LIMITS::PoleIntegrator>::calcValue() {
  // Tabulator parameters in std::vector<double> m_parameters
  // Pole::s_tabSigInd  = index for s
  // Pole::s_tabNobsInd = index for N(obs)
  m_function->setParameters( m_parameters );
  m_function->go();
  return m_function->result();
}


template<>
inline double Tabulator<LIMITS::PoleIntegrator>::interpolate( size_t ind ) const {
  double df    = deriv( ind, LIMITS::Pole::s_tabSigInd );  // derivative wrt S
  double d2f   = deriv2( ind, LIMITS::Pole::s_tabSigInd ); // derivative2 wrt S
  double f0    = m_tabValues[ind];                         // f() at discretized mean
  double x0;
  int    ix0   = calcParIndex(ind,LIMITS::Pole::s_tabSigInd);
  if ( ix0 < 0 ) {
    std::cout << "ERROR: calcParIndex return bad index!" << std::endl;
    return f0;
  }
  x0 = m_tabMin[LIMITS::Pole::s_tabSigInd] + static_cast<double>(ix0*m_tabStep[LIMITS::Pole::s_tabSigInd]);
  double x     = m_parameters[LIMITS::Pole::s_tabSigInd];
  double dx    = x-x0;
  double corr1 = df*dx;
  double corr2 = d2f*dx*dx/2.0;
  //  std::cout << "pole interp: " << ix0 << " , " << x << " , " << x0 << " -> " << f0 << " + " << corr1 << " + " << corr2 << std::endl;
  return  f0 + corr1 + corr2;
}
#endif

