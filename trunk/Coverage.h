#ifndef COVERAGE_H
#define COVERAGE_H
//
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <ctime>
#include <string>

#include "Range.h"
#include "Tabulated.h"
#include "Random.h"
#include "Pole.h"
class Coverage {
public:
  Coverage();
  virtual ~Coverage();
  //
  // Setup basic parameters
  //
  void setPole(Pole *pole) { m_pole = pole;}
  void setSeed(unsigned int r=0);
  void setNloops(int n) {m_nLoops = (n<1 ? 1:n);}
  // call these first
  void setSTrue(  double smin, double smax, double step);
  void setEffTrue(double emin, double emax, double step);
  void setBkgTrue(double bmin, double bmax, double step);
  //
  void setFixedSig(bool flag)  { m_fixedSig  = flag;}
  //
  // Sets the distributions for efficiency and background
  // CHECKPARAMS!!!
  void setEffDist(double mean,double sigma, DISTYPE dist=DIST_GAUS) { m_effMean=mean; m_effSigma=sigma; m_effDist=dist;}
  void setBkgDist(double mean,double sigma, DISTYPE dist=DIST_GAUS) { m_bkgMean=mean; m_bkgSigma=sigma; m_bkgDist=dist;}  //
  void setEffBkgCorr(double coef);
  bool checkEffBkgDists();
  //
  void initTabs();
  //
  void startTimer();
  void stopTimer();
  bool checkTimer(int dt=5);
  void printEstimatedTime(int nest);
  void endofRunTime();
  //
  void printSetup();
  void generateExperiment();   // creates a set of observables (eff,bkg,s)
  void doLoop();               // loops over all requested 'experiments'
  void doExpTest();            // loops over all requested 'experiments', no limit calc
  //
  void updateCoverage();	// Update coverage counters
  void resetCoverage();		// Reset dito
  //
  bool makeDumpName(std::string base, std::string & name);
  void setDumpBase(const char *base) {m_dumpFileNameBase = (base!=0 ? base:"");}
  const char *getDumpBase() {return m_dumpFileNameBase.c_str();}
  const char *getDumpName() {return m_dumpFileName.c_str();}
  //
  void collectStats(bool flag) { m_collectStats = flag; } // if true, then collect statistics
  void pushMeas();
  void pushLimits();
  void updateStatistics();	// Update statistics on limits
  void resetStatistics();	// Reset dito
  void calcStatistics();	// calculate collected stats
  void printStatistics();	// print calculated stuff
  //  void dumpExperiments(const char *name=0);
  void dumpExperiments(std::string name, bool dumpLimits=true);
  void dumpExperiments();
  void calcCoverage();		// Calculate coverage
  virtual void outputCoverageResult(const int flag=0);	// Output coverage
  //
  void setVerbose(int v=0) { m_verbose = v; }
  //
  double getMeasEff()  const {return m_measEff;}
  double getMeasBkg()  const {return m_measBkg;}
  double getMeasNobs() const {return m_measNobs;}

private:
  void calcStats(std::vector<double> & vec, double & average, double & variance);
  double calcStatsCorr(std::vector<double> & x, std::vector<double> & y);
  //
  int    m_verbose;
  //
  unsigned int m_rndSeed;
  Random m_rnd;
  // pointer to Pole class
  Pole  *m_pole;
  // loop ranges
  Range  m_sTrue;
  Range  m_effTrue;
  Range  m_bkgTrue;
  int    m_nLoops;   // number of loops
  // flags if true, the corresponding observable is kept fixed when generating experiments
  bool   m_fixedEff;
  bool   m_fixedBkg;
  bool   m_fixedSig;
  // Efficiency PDF
  double m_effMean;
  double m_effSigma;
  DISTYPE m_effDist;
  // Background PDF
  double m_bkgMean;
  double m_bkgSigma;
  DISTYPE m_bkgDist;
  // Correlation between bkg and eff (...pole)
  double m_beCorr;    // = r = correlation coeff (-1..1)
  double m_beCorrInv; // = sqrt(1-r*r);
  //
  // Signal True mean, poisson
  //
  double m_sTrueMean;
  //
  // Variables for current experiment
  //
  double m_measEff;     // "measured" efficiency (Gauss(m_effMean, m_effSigma))
  double m_measBkg;     // "measured" background (Gauss(m_bkgMean, m_bkgSigma))
  int    m_measNobs;    // "measured" signal (Poisson with mean m_sTrueMean)
  //
  bool   m_isInside;
  int    m_insideCount;
  int    m_totalCount;
  double m_coverage;
  double m_errCoverage;
  //
  // Some statistics
  //
  std::string m_dumpFileNameBase;
  std::string m_dumpFileName;
  bool m_collectStats;
  std::vector<double> m_UL;
  std::vector<double> m_LL;
  std::vector<double> m_effStat;
  std::vector<double> m_bkgStat;
  std::vector<double> m_nobsStat;
  double m_aveUL;
  double m_varUL;
  double m_aveLL;
  double m_varLL;
  double m_aveEff;
  double m_varEff;
  double m_aveBkg;
  double m_varBkg;
  double m_corrEffBkg;
  double m_aveNobs;
  double m_varNobs;

  // some timing stuff
  time_t m_startTime;  // start time
  time_t m_stopTime;   // stop time
  time_t m_estTime;    // estimated time for completion
};
#endif
