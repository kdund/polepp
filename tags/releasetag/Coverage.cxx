//
// This code contains classes for calculating the coverage of a CL
// belt algorithm.
//
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <sys/times.h>

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
  m_pole = 0;
  m_nLoops = 1;
  m_fixedEff = false;
  m_fixedBkg = false;
  m_fixedSig = false;

  // set various pointers to 0
  resetCoverage();
  resetStatistics();
}
Coverage::~Coverage() {
}

void Coverage::setSeed(unsigned int r) {
  m_rndSeed = r;
  m_rnd.SetSeed(r);
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

void Coverage::setCorr(double coef) {
  m_beCorr = (coef>=-1.0 ? (coef<=1.0 ? coef:0.0):0.0); // return 0 if bad
  m_beCorrInv = sqrt(1.0-m_beCorr*m_beCorr);
}

void Coverage::setUseCorr(bool flag) {
  m_useCorr = flag;
  if (flag) {
    if (m_fixedEff) std::cout << "WARNING: Using correlation - forcing random efficiency." << std::endl;
    if (m_fixedBkg) std::cout << "WARNING: Using correlation - forcing random background." << std::endl;
    m_fixedEff  = m_fixedBkg  = false;
  }
}

void Coverage::setFixedEff(bool flag) {
  m_fixedEff = flag;
  if (flag) {
    if (m_useCorr) std::cout << "WARNING: Efficiency fixed - will not use correlation." << std::endl;
    m_useCorr = false;
  }
}
void Coverage::setFixedBkg(bool flag) {
  m_fixedBkg = flag;
  if (flag) {
    if (m_useCorr) std::cout << "WARNING: Background fixed - will not use correlation." << std::endl;
    m_useCorr = false;
  }
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
  if (m_useCorr) {
    notOk = true;
    while (notOk) {
      double z1 = m_rnd.Gaus(0.0,1.0);
      double z2 = m_rnd.Gaus(0.0,1.0);
      m_measEff = m_effMean + z1*m_effSigma;
      m_measBkg = m_bkgMean + m_bkgSigma*(m_beCorr*z1 + m_beCorrInv*z2);
      notOk = ((m_measEff<0.0)||(m_measBkg<0.0));
    }
  } else {
    if (m_fixedEff) {
      m_measEff    = m_effMean;
    } else {
      if (m_useLogNormal) {
	m_measEff    = m_rnd.LogNormal(m_effMean,m_effSigma); //TMP
      } else {
	notOk = true;
	while (notOk) {
	  m_measEff    = m_rnd.Gaus(m_effMean,m_effSigma); // measured efficiency
	  notOk = (m_measEff<0.0);
	}
      }
    }
    if (m_fixedBkg) {
      m_measBkg    = m_bkgMean;
    } else {
      if (m_useLogNormal) {
	m_measBkg    = m_rnd.LogNormal(m_bkgMean,m_bkgSigma); //TMP
      } else {
	notOk = true;
	while (notOk) {
	  m_measBkg    = m_rnd.Gaus(m_bkgMean,m_bkgSigma); // measured background
	  notOk = (m_measBkg<0.0);
	}
      }
    }
  }
  if (m_fixedSig) {
    m_measNobs = static_cast<int>(m_effMean*m_sTrueMean+m_bkgMean+0.5);
  } else {
    m_measNobs   = m_rnd.Poisson(m_effMean*m_sTrueMean+m_bkgMean);
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
  m_errCoverage = 0;
  //
  if (m_totalCount>0) {
    m_coverage = static_cast<double>(m_insideCount)/static_cast<double>(m_totalCount);
    m_errCoverage = m_coverage*sqrt((1.0-m_pole->getCL())/static_cast<double>(m_insideCount));
  }
}

void Coverage::outputCoverageResult() { //output coverage result
  static bool firstCall=true;
  if (firstCall) {
    std::cout << "       Signal \tEfficiency\tBackground\tCoverage\tCov.err." << std::endl;
    std::cout << "       -----------------------------------------------------------------" << std::endl;
    firstCall = false;
  }
  std::cout << "DATA: " << std::fixed << std::setprecision(6)
	    << m_sTrueMean << "\t"
	    << m_effMean << "\t"
	    << m_bkgMean << "\t"
	    << m_coverage << "\t"
	    << m_errCoverage << std::endl;
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
void Coverage::updateStatistics() {
  if (m_collectStats) {
    m_UL.push_back(m_pole->getUpperLimit());
    m_LL.push_back(m_pole->getLowerLimit());
    m_effStat.push_back(m_measEff);
    m_bkgStat.push_back(m_measBkg);
    m_nobsStat.push_back(m_measNobs);
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
    calcStats(m_nobsStat,m_aveNobs,m_varNobs);
  }
}


void Coverage::printStatistics() {
  if (m_collectStats) {
    std::cout << "====== Statistics ======" << std::endl;
    std::cout << " N observed (mu,sig) =\t" << m_aveNobs << "\t" << sqrt(m_varNobs) << std::endl;
    std::cout << " efficiency (mu,sig) =\t" << m_aveEff  << "\t" << sqrt(m_varEff) << std::endl;
    std::cout << " background (mu,sig) =\t" << m_aveBkg  << "\t" << sqrt(m_varBkg) << std::endl;
    std::cout << " lower lim. (mu,sig) =\t" << m_aveLL   << "\t" << sqrt(m_varLL)  << std::endl;
    std::cout << " upper lim. (mu,sig) =\t" << m_aveUL   << "\t" << sqrt(m_varUL)  << std::endl;
    std::cout << "========================" << std::endl;
  }
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
  std::cout << " Efficiency sigma   : " << m_pole->getEffSigma() << std::endl;
  std::cout << " Efficiency fixed   : " << yesNo(m_fixedEff) << std::endl;
  //  std::cout << " Efficiency no dist : " << yesNo(m_noDistEff) << std::endl;
  std::cout << "----------------------------------------------\n";
  std::cout << " Background min     : " << m_bkgTrue.min() << std::endl;
  std::cout << " Background max     : " << m_bkgTrue.max() << std::endl;
  std::cout << " Background step    : " << m_bkgTrue.step() << std::endl;
  std::cout << " Background sigma   : " << m_pole->getBkgSigma() << std::endl;
  std::cout << " Background fixed   : " << yesNo(m_fixedBkg) << std::endl;
  //  std::cout << " Background no dist : " << yesNo(m_noDistBkg) << std::endl;
  std::cout << "----------------------------------------------\n";
  std::cout << " Correlated bkg,eff : " << yesNo(m_useCorr) << std::endl;
  std::cout << " Correlation coef.  : " << m_beCorr << std::endl;
  std::cout << "==============================================\n";
  //
  if (m_useLogNormal) {
    std::cout << "\nNOTE: USING LOG-NORMAL!!\n" << std::endl;
  }
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
  m_effSigma  = m_pole->getEffSigma();
  m_effNoDist = m_pole->isEffNoDist();
  m_bkgSigma  = m_pole->getBkgSigma();
  m_bkgNoDist = m_pole->isBkgNoDist();
  //
  m_pole->useLogNormal(m_useLogNormal); //TMP
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
	  m_pole->setNobserved(m_measNobs);
	  m_pole->setEffDist(m_measEff,m_effSigma,m_effNoDist);
	  m_pole->setBkgDist(m_measBkg,m_bkgSigma,m_bkgNoDist);
	  m_pole->setEffInt();         // will update if nescessary
	  m_pole->setBkgInt();
	  m_pole->initIntArrays();   // will update arrays if nescesarry
	  m_pole->initBeltArrays(true);
	  m_pole->analyseExperiment();  // calculate the limit of the given experiment
	  updateCoverage();     // update the coverage
	  updateStatistics();   // statistics (only if activated)
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