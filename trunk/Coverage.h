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
  // call these first - sets the scan ranges of the true s, eff and bkg.
  void setSTrue(  double smin, double smax, double step);
  void setEffTrue(double emin, double emax, double step);
  void setBkgTrue(double bmin, double bmax, double step);
  void setFixedSig(bool flag)  { m_fixedSig  = flag;}
  //
  void printSetup();
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
  void pushMeas(bool ok=true);
  void pushLimits();
  void updateStatistics(bool ok=true);	// Update statistics on limits
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
private:
  void calcStats(std::vector<double> & vec, double & average, double & variance);
  double calcStatsCorr(std::vector<double> & x, std::vector<double> & y);
  //
  int    m_verbose;
  //
  unsigned int m_rndSeed;
  RND::Random m_rnd;
  // pointer to Pole class
  Pole  *m_pole;
  // loop ranges
  Range<double>  m_sTrue;
  Range<double>  m_effTrue;
  Range<double>  m_bkgTrue;
  int    m_nLoops;   // number of loops
  // flags if true, the corresponding observable is kept fixed when generating experiments
  bool   m_fixedEff;
  bool   m_fixedBkg;
  bool   m_fixedSig;
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
  std::vector<double> m_effFrac;
  std::vector<double> m_bkgFrac;
  std::vector<double> m_sumProb;
  std::vector<double> m_status;
  //
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
  double m_aveStatus;
  double m_varStatus;

  // some timing stuff
  time_t m_startTime;  // start time
  time_t m_stopTime;   // stop time
  time_t m_estTime;    // estimated time for completion

  clock_t m_startClock;
  clock_t m_stopClock;
};
#endif
