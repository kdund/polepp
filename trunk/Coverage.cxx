//
// This code contains classes for calculating the coverage of a CL
// belt algorithm.
//
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <ctime>

#include "Coverage.h"

/////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////
inline char *yesNo(bool var) {
  static char yes[4] = {'Y','e','s',char(0)};
  static char no[3] = {'N','o',char(0)};
  return (var ? &yes[0]:&no[0]);
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

Coverage::Coverage() {
  m_saveExperiments = false;
  m_pole = 0;
  m_nLoops = 1;
  m_fixedSig = false;

  // set various pointers to 0
  resetCoverage();
  resetStatistics();
}
Coverage::~Coverage() {
}

void Coverage::setSeed(unsigned int r) {
  m_rndSeed = r;
  m_rnd.setSeed(r);
}

void Coverage::setSTrue(double low, double high, double step) {
  m_sTrue.setRange(low,high,step);
}

void Coverage::setEffTrue(double low, double high, double step) {
  m_effTrue.setRange(low,high,step);
}

void Coverage::setBkgTrue(double low, double high, double step) {
  m_bkgTrue.setRange(low,high,step);
}

void Coverage::setEffBkgCorr(double coef) {
  m_beCorr = (coef>=-1.0 ? (coef<=1.0 ? coef:0.0):0.0); // return 0 if bad
  m_beCorrInv = sqrt(1.0-m_beCorr*m_beCorr);
}

// void Coverage::setUseCorr(bool flag) {
//   m_useCorr = flag;
//   if (flag) {
//     if (m_fixedEff) std::cout << "WARNING: Using correlation - forcing random efficiency." << std::endl;
//     if (m_fixedBkg) std::cout << "WARNING: Using correlation - forcing random background." << std::endl;
//     m_fixedEff  = m_fixedBkg  = false;
//   }
// }

// void Coverage::setFixedEff(bool flag) {
//   m_fixedEff = flag;
//   //  if (flag) {
//   //    if (m_useCorr) std::cout << "WARNING: Efficiency fixed - will not use correlation." << std::endl;
//   //    m_useCorr = false;
//   //  }
// }
// void Coverage::setFixedBkg(bool flag) {
//   m_fixedBkg = flag;
//   //  if (flag) {
//   //    if (m_useCorr) std::cout << "WARNING: Background fixed - will not use correlation." << std::endl;
//   //    m_useCorr = false;
//   //  }
// }

bool Coverage::checkEffBkgDists() {
  bool change=false;
  // If ONLY one is DIST_GAUSCORR, make bkgDist == effDist
  if ( ((m_effDist == DIST_GAUSCORR) && (m_bkgDist != DIST_GAUSCORR)) ||
       ((m_effDist != DIST_GAUSCORR) && (m_bkgDist == DIST_GAUSCORR)) ) {
    m_bkgDist = m_effDist;
    change = true;
  }
  if (m_effDist != DIST_GAUSCORR) {
    m_beCorr = 0;
    change = true;
  }
  //
  if (m_effDist==DIST_NONE) {
    m_effSigma = 0.0;
    change = true;
  }
  if (m_bkgDist==DIST_NONE) {
    m_bkgSigma = 0.0;
    change = true;
  }
  return change;
}

//

void Coverage::generateExperiment() {
  //
  // Assumes:
  //
  //  Parameter	 | Truth      | Measured
  //-------------------------------------------------------
  //  efficiency | m_effMean  | Gauss(m_effmean,m_effSigma)
  //  background | m_bkgMean  | Gauss(m_bkgmean,m_bkgSigma)
  //  n obs.     | lambda     | Poisson(lambda)
  //-------------------------------------------------------
  // with lambda = m_effMean*m_sTrueMean + m_bkgMean
  //
  // If any of the above are fixed, they are set to the 'true' value.
  //
  // All parameters are forced positive
  //
  bool notOk = true;
  //
  switch (m_effDist) {
  case DIST_NONE:
    m_measEff = m_effMean;
    break;
  case DIST_GAUSCORR:
    notOk = true;
    while (notOk) {
      double z1 = m_rnd.gaus(0.0,1.0);
      double z2 = m_rnd.gaus(0.0,1.0);
      m_measEff = m_effMean + z1*m_effSigma;
      m_measBkg = m_bkgMean + m_bkgSigma*(m_beCorr*z1 + m_beCorrInv*z2);
      notOk = ((m_measEff<0.0)||(m_measBkg<0.0));
    }
    break;
  case DIST_LOGN:
    m_measEff = m_rnd.logNormal(m_effMean,m_effSigma);
    break;
  case DIST_FLAT:
    m_measEff = m_rnd.flat(m_effMean, m_effSigma);
    break;
  case DIST_GAUS:
    notOk = true;
    while (notOk) {
      m_measEff    = m_rnd.gaus(m_effMean,m_effSigma); // measured efficiency
      notOk = (m_measEff<0.0);
    }
    break;
  default: // ERROR STATE
    m_measEff = m_effMean;
    break;
  }
  //
  switch (m_bkgDist) {
  case DIST_NONE:
    m_measBkg = m_bkgMean;
    break;
  case DIST_LOGN:
    m_measBkg = m_rnd.logNormal(m_bkgMean,m_bkgSigma);
    break;
  case DIST_FLAT:
    m_measBkg = m_rnd.flat(m_bkgMean, m_bkgSigma);
    break;
  case DIST_GAUS:
    notOk = true;
    while (notOk) {
      m_measBkg    = m_rnd.gaus(m_bkgMean,m_bkgSigma); // measured background
      notOk = (m_measBkg<0.0);
    }
  case DIST_GAUSCORR: // already taken care of above
    break;
  default: // ERROR STATE
    m_measBkg = m_bkgMean;
    break;
  }
  //
  if (m_fixedSig) {
    m_measNobs = static_cast<int>(m_effMean*m_sTrueMean+m_bkgMean+0.5);
  } else {
    m_measNobs   = m_rnd.poisson(m_effMean*m_sTrueMean+m_bkgMean);
  }
  if (m_saveExperiments) {
    m_allNobs.push_back(m_measNobs);
    m_allEff.push_back(m_measEff);
    m_allBkg.push_back(m_measBkg);
  }
  //
}

void Coverage::updateCoverage() {
  m_totalCount++;
  double ll = m_pole->getLowerLimit();
  double ul = m_pole->getUpperLimit();
  double s  = m_sTrueMean;
  if (s>=ll ? (s<=ul ? true:false) : false) {
    m_insideCount++;
    m_isInside = true;
  } else {
    m_isInside = false;
  }
  if (m_verbose>2) {
    std::cout << "LIMIT: " << yesNo(m_isInside) << "\t";
    m_pole->printLimit();
  }
}

void Coverage::calcCoverage() {
  m_coverage = 0;
  m_errCoverage = -1.0;
  //
  if (m_totalCount>0) {
    m_coverage = static_cast<double>(m_insideCount)/static_cast<double>(m_totalCount);
    if (m_insideCount>0) {
      m_errCoverage = m_coverage*sqrt((1.0-m_pole->getCL())/static_cast<double>(m_insideCount));
    }
  }
}

void Coverage::outputCoverageResult(const int flag) { //output coverage result
  static bool firstCall=true;
  std::string header;
  if (firstCall) {
    std::cout << "       Signal \tEfficiency\t\t\tBackground\t\tCorrelation\tCoverage\t\tLoops\tMax loops" << std::endl;
    std::cout << "      --------------------------------------------------------------------------------------------" << std::endl;
    firstCall = false;
  }
  if (flag==1) { // status (intermediate result)
    header = "STATUS: ";
  } else {       // final result
    header = "DATA: ";
  }
  std::cout << header.c_str();
  coutFixed(6,m_sTrueMean); std::cout << "\t";
  coutFixed(6,m_effMean); std::cout << "\t";
  coutFixed(6,m_effSigma); std::cout << "\t";
  coutFixed(6,m_bkgMean); std::cout << "\t";
  coutFixed(6,m_bkgSigma); std::cout << "\t";
  coutFixed(6,m_beCorr); std::cout << "\t";
  coutFixed(6,m_coverage); std::cout << "\t";
  coutFixed(6,m_errCoverage); std::cout << "\t";
  coutFixed(6,m_totalCount); std::cout << "\t";
  coutFixed(6,m_nLoops); std::cout << std::endl;
}


void Coverage::resetCoverage() {
  m_coverage = 0;
  m_errCoverage = 0;
  m_insideCount = 0;
  m_totalCount = 0;
  m_isInside = false;
}

void Coverage::calcStats(std::vector<double> & vec, double & average, double & variance) {
  unsigned int size = vec.size();
  average = 0;
  variance = 0;
  if (size>0) {
    double sum =0;
    double sum2=0;
    for (unsigned int i=0; i<size; i++) {
      sum  += vec[i];
      sum2 += vec[i]*vec[i];
    }
    double n = static_cast<double>(size);
    average = sum/n;
    variance = (sum2 - (sum*sum)/n)/(n-1.0);
  }
}

double Coverage::calcStatsCorr(std::vector<double> & x, std::vector<double> & y) {
  unsigned int size = x.size();
  double rval = 0;
  if (size>0) {
    double sumxy =0;
    double sumx=0;
    double sumy=0;
    for (unsigned int i=0; i<size; i++) {
      sumxy += x[i]*y[i];
      sumx  += x[i];
      sumy  += y[i];
    }
    double n = static_cast<double>(size);
    rval = (sumxy - ((sumx*sumy)/n))/(n-1.0);
  }
  return rval;
}

void Coverage::pushLimits() {
  m_UL.push_back(m_pole->getUpperLimit());
  m_LL.push_back(m_pole->getLowerLimit());
}

void Coverage::pushMeas() {
  m_effStat.push_back(m_measEff);
  m_bkgStat.push_back(m_measBkg);
  m_nobsStat.push_back(m_measNobs);
}

void Coverage::updateStatistics() {
  if (m_collectStats) {
    pushLimits();
    pushMeas();
  }
}

void Coverage::resetStatistics() {
  m_aveUL = 0;
  m_varUL = 0;
  m_aveLL = 0;
  m_varLL = 0;
  m_aveEff = 0;
  m_varEff = 0;
  m_aveBkg = 0;
  m_varBkg = 0;
  m_UL.clear();
  m_LL.clear();
  m_effStat.clear();
  m_bkgStat.clear();
  m_nobsStat.clear();
}

void Coverage::calcStatistics() {
  if (m_collectStats) {
    calcStats(m_UL,m_aveUL,m_varUL);
    calcStats(m_LL,m_aveLL,m_varLL);
    calcStats(m_effStat,m_aveEff,m_varEff);
    calcStats(m_bkgStat,m_aveBkg,m_varBkg);
    //
    m_corrEffBkg = calcStatsCorr(m_effStat,m_bkgStat);
    if ((m_effDist!=DIST_NONE) && (m_varEff>0)) {
      m_corrEffBkg = m_corrEffBkg / sqrt(m_varEff);
    }
    if ((m_bkgDist!=DIST_NONE) && (m_varBkg>0)) {
      m_corrEffBkg = m_corrEffBkg / sqrt(m_varBkg);
    }
    //
    calcStats(m_nobsStat,m_aveNobs,m_varNobs);
  }
}


void Coverage::printStatistics() {
  if (m_collectStats) {
    std::cout << "====== Statistics ======" << std::endl;
    std::cout << " N observed (mu,sig) =\t";
    coutFixed(6,m_aveNobs); std::cout << "\t";
    coutFixed(6,sqrt(m_varNobs)); std::cout << std::endl;
    std::cout << " efficiency (mu,sig) =\t";
    coutFixed(6,m_aveEff); std::cout << "\t";
    coutFixed(6,sqrt(m_varEff)); std::cout << std::endl;
    std::cout << " background (mu,sig) =\t";
    coutFixed(6,m_aveBkg); std::cout << "\t";
    coutFixed(6,sqrt(m_varBkg)); std::cout << std::endl;
    std::cout << " correlation coeff   =\t";
    coutFixed(6,m_corrEffBkg); std::cout << std::endl;
    std::cout << " lower lim. (mu,sig) =\t";
    coutFixed(6,m_aveLL); std::cout << "\t";
    coutFixed(6,sqrt(m_varLL)); std::cout << std::endl;
    std::cout << " upper lim. (mu,sig) =\t";
    coutFixed(6,m_aveUL); std::cout << "\t";
    coutFixed(6,sqrt(m_varUL)); std::cout << std::endl;
    std::cout << "========================" << std::endl;
  }
}

void Coverage::dumpExperiments(std::string name) {
  unsigned i,sz;
  sz = m_effStat.size();
  if (sz==0) return;
  std::ostream  *os=0;
  bool osdel=false;
  //
  if (name.size()>0) {
    std::cout << "Dumping experiments to file: " << name << std::endl;
    os = new std::ofstream(name.c_str());
    osdel = true;
  } else {
    os = &std::cout;
  }
  *os << "#" << std::endl;
  *os << "# N         = " << m_nLoops << std::endl;
  *os << "# s_true    = " << m_sTrue.min() << std::endl;
  *os << "# eff       = " << m_effTrue.min() << std::endl;
  *os << "# eff sigma = " << m_effSigma << std::endl;
  *os << "# bkg       = " << m_bkgTrue.min() << std::endl;
  *os << "# bkg sigma = " << m_bkgSigma << std::endl;
  *os << "#" << std::endl;
  *os << "# N_obs   Efficiency     Background" << std::endl;
  *os << "#-----------------------------------" << std::endl;
#if __GNUC__ > 2
  for (i=0; i<sz; i++) {
    *os << std::fixed
	<< std::setprecision(0)
	<< m_nobsStat[i] << '\t'
	<< std::setprecision(6)
	<< m_effStat[i] << '\t'
        << m_bkgStat[i] << std::endl;
  }
#else
  for (i=0; i<sz; i++) {
    *os << m_nobsStat[i] << '\t'
	<< m_effStat[i] << '\t'
        << m_bkgStat[i] << std::endl;
  }
#endif
  *os << "#EOF" << std::endl;
  if (osdel) delete os;
}

void Coverage::startTimer() {
  //
  static char tstamp[32];
  m_stopTime = 0;
  m_estTime  = 0;
  //
  time(&m_startTime);
  struct tm *timeStruct = localtime(&m_startTime);
  strftime(tstamp,32,"%d/%m/%Y %H:%M:%S",timeStruct);
  std::cout << "Start of run: " << tstamp << std::endl;
}

bool Coverage::checkTimer(int dt) {
  time_t chk;
  int delta;
  time(&chk);
  delta = int(chk - m_startTime);
  return (delta>dt); // true if full time has passed
}

void Coverage::stopTimer() {
  time(&m_stopTime);
}

void Coverage::endofRunTime() {
  static char tstamp[32];
  time_t eorTime;
  //
  time(&eorTime);
  struct tm *timeStruct = localtime(&eorTime);
  strftime(tstamp,32,"%d/%m/%Y %H:%M:%S",timeStruct);
  std::cout << "\nEnd of run: " << tstamp << std::endl;
}

void Coverage::printEstimatedTime(int nest) {
  static char tstamp[32];
  if (m_startTime<=m_stopTime) {
    time_t loopTime  = m_stopTime - m_startTime;
    time_t deltaTime = (m_nLoops*m_sTrue.n()*m_effTrue.n()*m_bkgTrue.n()*loopTime)/nest;
    //
    m_estTime = deltaTime + m_startTime;
    struct tm *timeStruct;
    timeStruct = localtime(&m_estTime);
    strftime(tstamp,32,"%d/%m/%Y %H:%M:%S",timeStruct);
    timeStruct = localtime(&deltaTime);
    int hours,mins,secs;
    hours = static_cast<int>(deltaTime/3600);
    mins = (deltaTime - hours*3600)/60;
    secs =  deltaTime - hours*3600-mins*60;
    std::cout << "Estimated end of run: " << tstamp;
    std::cout << " ( " << hours << "h " << mins << "m " << secs << "s" << " )\n" << std::endl;
  }
}

void Coverage::printSetup() {
  std::cout << "\n";
  std::cout << "==============C O V E R A G E=================\n";
  std::cout << " Random seed        : " << m_rndSeed << std::endl;
  std::cout << " Number of loops    : " << m_nLoops << std::endl;
  std::cout << " Collect statistics : " << yesNo(m_collectStats) << std::endl;
  std::cout << "----------------------------------------------\n";
  std::cout << " Signal min         : " << m_sTrue.min() << std::endl;
  std::cout << " Signal max         : " << m_sTrue.max() << std::endl;
  std::cout << " Signal step        : " << m_sTrue.step() << std::endl;
  std::cout << " Signal N           : " << m_sTrue.n() << std::endl;
  std::cout << " Signal fixed       : " << yesNo(m_fixedSig) << std::endl;
  std::cout << "----------------------------------------------\n";
  std::cout << " Efficiency min     : " << m_effTrue.min() << std::endl;
  std::cout << " Efficiency max     : " << m_effTrue.max() << std::endl;
  std::cout << " Efficiency step    : " << m_effTrue.step() << std::endl;
  std::cout << " Efficiency sigma   : " << m_effSigma << std::endl;
  std::cout << " Efficiency dist    : " << distTypeStr(m_effDist) << std::endl;
  std::cout << "----------------------------------------------\n";
  std::cout << " Background min     : " << m_bkgTrue.min() << std::endl;
  std::cout << " Background max     : " << m_bkgTrue.max() << std::endl;
  std::cout << " Background step    : " << m_bkgTrue.step() << std::endl;
  std::cout << " Background sigma   : " << m_bkgSigma << std::endl;
  std::cout << " Background dist    : " << distTypeStr(m_bkgDist) << std::endl;
  std::cout << "----------------------------------------------\n";
  std::cout << " Correlated bkg,eff : " << yesNo((m_bkgDist==DIST_GAUSCORR)) << std::endl;
  std::cout << " Correlation coef.  : " << m_beCorr << std::endl;
  std::cout << "==============================================\n";
  //
  if (m_collectStats) {
    std::cout << "\nNOTE: collecting statistics -> not optimum for coverage studies\n" << std::endl;
  }
}

void Coverage::doLoop() {
  if (m_pole==0) return;
  //
  // Number of loops to time for an estimate of the total time
  //  int nest = (m_nLoops>100 ? 100:(m_nLoops>10 ? 10:1));
  int nest = 0; // new - times a given period and not a given number of loops
  //
  int is,ie,ib,j;
  bool first = true;
  bool timingDone = false;
  //
  m_pole->setCoverage(true);
  //
  for (is=0; is<m_sTrue.n(); is++) { // loop over all s_true
    m_sTrueMean = m_sTrue.getVal(is);
    m_pole->setTrueSignal(m_sTrueMean);
    for (ie=0; ie<m_effTrue.n(); ie++) { // loop over eff true
      m_effMean = m_effTrue.getVal(ie);
      for (ib=0; ib<m_bkgTrue.n(); ib++) { // loop over bkg true
	m_bkgMean = m_bkgTrue.getVal(ib);

	resetCoverage();   // for each s_true, reset coverage
	resetStatistics(); // dito, statistics
    //
	for (j=0; j<m_nLoops; j++) {  // get stats
	  if (first) {
	    startTimer();       // used for estimating the time for the run
	    first = false;
	  }
	  generateExperiment(); // generate pseudoexperiment using given ditributions of signal,bkg and eff.
	  m_pole->setNobserved(m_measNobs); // set values from generated experiment
	  m_pole->setEffMeas(m_measEff,m_effSigma,m_effDist);
	  m_pole->setBkgMeas(m_measBkg,m_bkgSigma,m_bkgDist);
	  m_pole->setEffBkgCorr(m_beCorr); // always the same...
	  m_pole->setEffInt();         // reset the integral ranges
	  m_pole->setBkgInt();
	  m_pole->initIntArrays();     // will update arrays if nescesarry
	  m_pole->initBeltArrays(true);
	  m_pole->analyseExperiment(); // calculate the limit of the given experiment
	  updateCoverage();            // update the coverage
	  updateStatistics();          // statistics (only if activated)
	  if (!timingDone) {
	    nest++;
	    if (checkTimer(5)) {
	      timingDone = true;
	      stopTimer();
	      printEstimatedTime(nest);
	    }
	  }
	  //	  if ((!timingDone) && (j==nest-1)) { // calculate the estimated run time.
	  //	    stopTimer();
	  //	    printEstimatedTime(nest);
	  //	    timingDone = true;
	  //	  }
	}
	calcCoverage();         // calculate coverage and its uncertainty
	outputCoverageResult(); // print the result
	calcStatistics();       // dito for the statistics...
	printStatistics();
      }
    }
  }
  endofRunTime();
}

void Coverage::doExpTest() {
  //
  int is,ie,ib,j;
  //
  for (is=0; is<m_sTrue.n(); is++) { // loop over all s_true
    m_sTrueMean = m_sTrue.getVal(is);
    for (ie=0; ie<m_effTrue.n(); ie++) { // loop over eff true
      m_effMean = m_effTrue.getVal(ie);
      for (ib=0; ib<m_bkgTrue.n(); ib++) { // loop over bkg true
	m_bkgMean = m_bkgTrue.getVal(ib);
	//
	resetStatistics(); // reset statistics
	//
	for (j=0; j<m_nLoops; j++) {  // get stats
	  generateExperiment(); // generate pseudoexperiment using given ditributions of signal,bkg and eff.
	  pushMeas();           // save 'measured' data
	}
	calcStatistics();       // 
	printStatistics();
      }
    }
  }
}
