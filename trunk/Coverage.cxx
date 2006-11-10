//
// This code contains classes for calculating the coverage of a CL
// belt algorithm.
//
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <sstream>

#include "Tools.h"
#include "Coverage.h"

Coverage::Coverage() {
  m_doneOneLoop = false;
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

// bool Coverage::checkEffBkgDists() {
//   bool change=false;
//   // If ONLY one is PDF::DIST_GAUS2D, make bkgDist == effDist
//   if ( ((m_effDist == PDF::DIST_GAUS2D) && (m_bkgDist != PDF::DIST_GAUS2D)) ||
//        ((m_effDist != PDF::DIST_GAUS2D) && (m_bkgDist == PDF::DIST_GAUS2D)) ) {
//     m_bkgDist = m_effDist;
//     change = true;
//   }
//   if (m_effDist != PDF::DIST_GAUS2D) {
//     m_beCorr = 0;
//     change = true;
//   }
//   //
//   if (m_effDist==PDF::DIST_NONE) {
//     m_effSigma = 0.0;
//     change = true;
//   }
//   if (m_bkgDist==PDF::DIST_NONE) {
//     m_bkgSigma = 0.0;
//     change = true;
//   }
//   return change;
// }

//

// void Coverage::generateExperiment() {
//   //
//   // Assumes:
//   //
//   //  Parameter	 | Truth      | Measured
//   //-------------------------------------------------------
//   //  efficiency | m_effMean  | Gauss(m_effmean,m_effSigma)
//   //  background | m_bkgMean  | Gauss(m_bkgmean,m_bkgSigma)
//   //  n obs.     | lambda     | Poisson(lambda)
//   //-------------------------------------------------------
//   // with lambda = m_effMean*m_sTrueMean + m_bkgMean
//   //
//   // If any of the above are fixed, they are set to the 'true' value.
//   //
//   // All parameters are forced positive
//   //
//   bool notOk = true;
//   //
//   switch (m_effDist) {
//   case PDF::DIST_NONE:
//     m_measEff = m_effMean;
//     break;
//   case PDF::DIST_GAUS2D:
//     notOk = true;
//     while (notOk) {
//       double z1 = m_rnd.gauss(0.0,1.0);
//       double z2 = m_rnd.gauss(0.0,1.0);
//       m_measEff = m_effMean + z1*m_effSigma;
//       m_measBkg = m_bkgMean + m_bkgSigma*(m_beCorr*z1 + m_beCorrInv*z2);
//       notOk = ((m_measEff<0.0)||(m_measBkg<0.0));
//     }
//     break;
//   case PDF::DIST_LOGN:
//     m_measEff = m_rnd.logNormal(m_effMean,m_effSigma);
//     break;
//   case PDF::DIST_FLAT:
//     m_measEff = m_rnd.flat(m_effMean, m_effSigma);
//     break;
//   case PDF::DIST_GAUS:
//     notOk = true;
//     while (notOk) {
//       m_measEff    = m_rnd.gauss(m_effMean,m_effSigma); // measured efficiency
//       notOk = (m_measEff<0.0);
//     }
//     break;
//   default: // ERROR STATE
//     m_measEff = m_effMean;
//     break;
//   }
//   //
//   switch (m_bkgDist) {
//   case PDF::DIST_NONE:
//     m_measBkg = m_bkgMean;
//     break;
//   case PDF::DIST_LOGN:
//     m_measBkg = m_rnd.logNormal(m_bkgMean,m_bkgSigma);
//     break;
//   case PDF::DIST_FLAT:
//     m_measBkg = m_rnd.flat(m_bkgMean, m_bkgSigma);
//     break;
//   case PDF::DIST_GAUS:
//     notOk = true;
//     while (notOk) {
//       m_measBkg    = m_rnd.gauss(m_bkgMean,m_bkgSigma); // measured background
//       notOk = (m_measBkg<0.0);
//     }
//   case PDF::DIST_GAUS2D: // already taken care of above
//     break;
//   default: // ERROR STATE
//     m_measBkg = m_bkgMean;
//     break;
//   }
//   //
//   if (m_fixedSig) {
//     m_measNobs = static_cast<int>(m_effMean*m_sTrueMean+m_bkgMean+0.5);
//   } else {
//     m_measNobs   = m_rnd.poisson(m_effMean*m_sTrueMean+m_bkgMean);
//   }
//   //
// }

void Coverage::updateCoverage() {
  m_totalCount++;
  if (m_pole->truthCovered()) {
    m_insideCount++;
    m_isInside = true;
  } else {
    m_isInside = false;
  }
  if (m_verbose>2) {
    std::cout << "LIMIT: " << TOOLS::yesNo(m_isInside) << "\t";
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
      //      m_errCoverage = m_coverage*( (sqrt((1.0-m_pole->getCL())/static_cast<double>(m_insideCount))));
      m_errCoverage = sqrt(m_coverage*(1.0-m_coverage)/static_cast<double>(m_totalCount));
    }
  }
}

void Coverage::outputCoverageResult(const int flag) { //output coverage result
  static bool firstCall=true;
  std::string header;
  if (firstCall) {
    //
    std::cout << "#==================================================================================================================" << std::endl;
    std::cout << "#         Signal    |     Efficiency       |      Background       |      Coverage         |    Loops    |    Time " << std::endl;
    std::cout << "#                   | mean        sigma    |  mean        sigma    |  mean        sigma    | done    max |    [ms] " << std::endl;
    std::cout << "#==================================================================================================================" << std::endl;
    firstCall = false;
  }
  if (flag==1) { // status (intermediate result)
    header = "STATUS: ";
  } else {       // final result
    header = "COVERAGE: ";
  }
  std::cout << header.c_str();
  TOOLS::coutFixed(6,m_pole->getTrueSignal()); std::cout << "    ";
  TOOLS::coutFixed(6,m_pole->getEffPdfMean()); std::cout << "    ";
  TOOLS::coutFixed(6,m_pole->getEffPdfSigma()); std::cout << "    ";
  TOOLS::coutFixed(6,m_pole->getBkgPdfMean()); std::cout << "    ";
  TOOLS::coutFixed(6,m_pole->getBkgPdfSigma()); std::cout << "    ";
  //TOOLS::coutFixed(6,m_beCorr); std::cout << "    ";
  TOOLS::coutFixed(6,m_coverage); std::cout << "    ";
  TOOLS::coutFixed(6,m_errCoverage); std::cout << "    ";
  TOOLS::coutFixed(6,m_totalCount); std::cout << "    ";
  TOOLS::coutFixed(6,m_nLoops); std::cout << "      ";
  TOOLS::coutFixed(2,m_timer.getUsedClock()); std::cout << std::endl;
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

void Coverage::pushLimits(bool ok) {
  m_UL.push_back(m_pole->getUpperLimit());
  m_LL.push_back(m_pole->getLowerLimit());
  m_sumProb.push_back(m_pole->getSumProb());
  m_status.push_back((ok ? 1.0:0.0));
}

void Coverage::pushMeas() {
  m_effStat.push_back(m_pole->getEffObs());
  m_bkgStat.push_back(m_pole->getBkgObs());
  m_nobsStat.push_back(m_pole->getNObserved());
//   if ((m_pole->getEffPdfSigma()>0) && (m_pole->getEffPdfDist() != PDF::DIST_NONE)) {
//     m_effFrac.push_back(m_pole->getEffPdfSigma()/m_pole->getEffObs());
//   } else {
//     m_effFrac.push_back(-1.0);
//   }
//   if ((m_pole->getBkgPdfSigma()>0) && (m_pole->getBkgPdfDist() != PDF::DIST_NONE)) {
//     m_bkgFrac.push_back(m_pole->getBkgPdfSigma()/m_pole->getBkgObs());
//   } else {
//     m_bkgFrac.push_back(-1.0);
//   }
}

void Coverage::updateStatistics(bool ok) {
  if (m_collectStats) {
    pushLimits(ok);
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
  m_effFrac.clear();
  m_bkgFrac.clear();
  m_nobsStat.clear();
  m_sumProb.clear();
  m_status.clear();
}

void Coverage::calcStatistics() {
  if (m_collectStats) {
    calcStats(m_UL,m_aveUL,m_varUL);
    calcStats(m_LL,m_aveLL,m_varLL);
    calcStats(m_effStat,m_aveEff,m_varEff);
    calcStats(m_bkgStat,m_aveBkg,m_varBkg);
    calcStats(m_status,m_aveStatus, m_varStatus);
    //
    m_corrEffBkg = calcStatsCorr(m_effStat,m_bkgStat);
    if ((m_pole->getEffPdfDist()!=PDF::DIST_NONE) && (m_varEff>0)) {
      m_corrEffBkg = m_corrEffBkg / sqrt(m_varEff);
    }
    if ((m_pole->getBkgPdfDist()!=PDF::DIST_NONE) && (m_varBkg>0)) {
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
    TOOLS::coutFixed(6,m_aveNobs); std::cout << "\t";
    TOOLS::coutFixed(6,sqrt(m_varNobs)); std::cout << std::endl;
    std::cout << " efficiency (mu,sig) =\t";
    TOOLS::coutFixed(6,m_aveEff); std::cout << "\t";
    TOOLS::coutFixed(6,sqrt(m_varEff)); std::cout << std::endl;
    std::cout << " background (mu,sig) =\t";
    TOOLS::coutFixed(6,m_aveBkg); std::cout << "\t";
    TOOLS::coutFixed(6,sqrt(m_varBkg)); std::cout << std::endl;
    std::cout << " correlation coeff   =\t";
    TOOLS::coutFixed(6,m_corrEffBkg); std::cout << std::endl;
    std::cout << " lower lim. (mu,sig) =\t";
    TOOLS::coutFixed(6,m_aveLL); std::cout << "\t";
    TOOLS::coutFixed(6,sqrt(m_varLL)); std::cout << std::endl;
    std::cout << " upper lim. (mu,sig) =\t";
    TOOLS::coutFixed(6,m_aveUL); std::cout << "\t";
    TOOLS::coutFixed(6,sqrt(m_varUL)); std::cout << std::endl;
    std::cout << " Success rate        =\t";
    TOOLS::coutFixed(6,m_aveStatus); std::cout << "\t";
    TOOLS::coutFixed(6,sqrt(m_varStatus)); std::cout << std::endl;
    std::cout << "========================" << std::endl;
  }
}

bool Coverage::makeDumpName(std::string base, std::string & name) {
  bool rval=false;
  if (base=="") return rval;
  int eind  = static_cast<int>(100.0*m_pole->getEffObs());
  int deind = static_cast<int>(100.0*m_pole->getEffPdfSigma());
  int bind  = static_cast<int>(100.0*m_pole->getBkgObs());
  int dbind = static_cast<int>(100.0*m_pole->getBkgPdfSigma());
  int sind  = static_cast<int>(100.0*m_pole->getTrueSignal());
  int we = (eind>999 ? 4:3);
  int wb = (bind>999 ? 4:3);
  int ws = (sind>999 ? 4:3);
  std::ostringstream ostr;
  ostr << std::fixed << std::setfill('0') << "_s" << std::setw(ws) << sind
       << "_e" << PDF::distTypeStr(m_pole->getEffPdfDist())
       << "-"  << std::setw(we) << eind
       << "-"  << std::setw(we) << deind
       << "_b" << PDF::distTypeStr(m_pole->getBkgPdfDist())
       << "-"  << std::setw(wb) << bind
       << "-"  << std::setw(wb) << dbind;
  name = base+ostr.str()+".dat";
  rval=true;
  return rval;
}

void Coverage::dumpExperiments(bool limits) {
  if (makeDumpName(m_dumpFileNameBase,m_dumpFileName)) {
    std::cout << "Dumping data to file : " << m_dumpFileName << std::endl;
    dumpExperiments(m_dumpFileName,true);
  }
}

void Coverage::dumpExperiments(std::string name, bool limits) {
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
  bool dumpLimits=false;
  if (limits) {
    dumpLimits = (m_UL.size() == sz);
  }
  //
  *os << "#" << std::endl;
  *os << "# N             = " << m_nLoops << std::endl;
  *os << "# s_true        = " << m_sTrue.min() << std::endl;
  *os << "# eff           = " << m_effTrue.min() << std::endl;
  *os << "# eff sigma     = " << m_pole->getEffPdfSigma() << std::endl;
  *os << "# eff dist      = " << distTypeStr(m_pole->getEffPdfDist()) << std::endl;
  *os << "# bkg           = " << m_bkgTrue.min() << std::endl;
  *os << "# bkg sigma     = " << m_pole->getBkgPdfSigma() << std::endl;
  *os << "# bkg dist      = " << distTypeStr(m_pole->getBkgPdfDist()) << std::endl;
  *os << "# lhRatio       = " << (m_pole->usesMBT() ? "MBT":"FHC2") << std::endl;
  //  *os << "# corr.         = " << m_beCorr << std::endl;
  *os << "# coverage      = " << m_coverage << std::endl;
  *os << "# coverage unc. = " << m_errCoverage << std::endl;
  *os << "#" << std::endl;
  *os << "# N_obs   Efficiency     Background" << std::endl;
  *os << "#-----------------------------------" << std::endl;

#if __GNUC__ > 2
  for (i=0; i<sz; i++) {

    *os << "COVSTAT: "
        << std::fixed
	<< std::setprecision(0)
	<< m_nobsStat[i] << '\t'
	<< m_pole->getEffPdfDist() << '\t'
	<< std::setprecision(6)
	<< m_effStat[i] << '\t'
	<< m_pole->getEffPdfSigma() << '\t'
	<< std::setprecision(0)
	<< m_pole->getBkgPdfDist() << '\t'
	<< std::setprecision(6)
        << m_bkgStat[i] << '\t'
        << m_pole->getBkgPdfSigma() << '\t';
    if (dumpLimits) {
      *os << std::setprecision(0)
          << m_status[i] << '\t'
          << std::setprecision(2)
          << m_LL[i] << '\t'
          << m_UL[i] << '\t'
          << m_sumProb[i] << '\t';
    }
    *os << std::endl;
  }
#else
  for (i=0; i<sz; i++) {
    *os << "COVSTAT: "
        << m_nobsStat[i] << '\t'
	<< m_effStat[i] << '\t'
        << m_bkgStat[i] <<
    if (dumpLimits) {

      *os << '\t'
	  << m_LL[i] << '\t'
	  << m_UL[i];
    }
    *os << m_std::endl;
  }
#endif
  *os << "#EOF" << std::endl;
  if (osdel) delete os;
}

void Coverage::printSetup() {
  std::cout << "\n";
  std::cout << "==============C O V E R A G E=================\n";
  std::cout << " Random seed        : " << m_rndSeed << std::endl;
  std::cout << " Number of loops    : " << m_nLoops << std::endl;
  std::cout << " Collect statistics : " << TOOLS::yesNo(m_collectStats) << std::endl;
  std::cout << "----------------------------------------------\n";
  std::cout << " Signal min         : " << m_sTrue.min() << std::endl;
  std::cout << " Signal max         : " << m_sTrue.max() << std::endl;
  std::cout << " Signal step        : " << m_sTrue.step() << std::endl;
  std::cout << " Signal N           : " << m_sTrue.n() << std::endl;
  std::cout << " Signal fixed       : " << TOOLS::yesNo(m_fixedSig) << std::endl;
  std::cout << "----------------------------------------------\n";
  std::cout << " Efficiency min     : " << m_effTrue.min() << std::endl;
  std::cout << " Efficiency max     : " << m_effTrue.max() << std::endl;
  std::cout << " Efficiency step    : " << m_effTrue.step() << std::endl;
  std::cout << " Efficiency sigma   : " << m_pole->getEffPdfSigma() << std::endl;
  std::cout << " Efficiency dist    : " << PDF::distTypeStr(m_pole->getEffPdfDist()) << std::endl;
  std::cout << "----------------------------------------------\n";
  std::cout << " Background min     : " << m_bkgTrue.min() << std::endl;
  std::cout << " Background max     : " << m_bkgTrue.max() << std::endl;
  std::cout << " Background step    : " << m_bkgTrue.step() << std::endl;
  std::cout << " Background sigma   : " << m_pole->getBkgPdfSigma() << std::endl;
  std::cout << " Background dist    : " << PDF::distTypeStr(m_pole->getBkgPdfDist()) << std::endl;
  std::cout << "----------------------------------------------\n";
  std::cout << " Correlated bkg,eff : " << TOOLS::yesNo((m_pole->getBkgPdfDist()==PDF::DIST_GAUS2D)) << std::endl;
  //  std::cout << " Correlation coef.  : " << m_beCorr << std::endl;
  std::cout << "==============================================\n";
  //
  if (m_collectStats) {
    std::cout << "\nNOTE: collecting statistics -> not optimum for coverage studies\n" << std::endl;
  }
}

void Coverage::doLoop() {
  m_doneOneLoop = false;
  if (m_pole==0) return;
  //
  //
  // Number of loops to time for an estimate of the total time
  //  int nest = (m_nLoops>100 ? 100:(m_nLoops>10 ? 10:1));
  //
  int nest = 0; // new - times a given period and not a given number of loops
  //
  int is,ie,ib,j;
  bool first = true;
  bool timingDone = false;
  int nWarnings=0;
  int nTotal=m_sTrue.n()*m_effTrue.n()*m_bkgTrue.n()*m_nLoops;
  const int maxWarnings=10;
  //
  // If we want to collect statistics, the full limit is required.
  // For just finding the coverage, a faster method is used.
  //
  m_pole->setUseCoverage(!m_collectStats); // make full limits when collecting statistics
  //
  // Start timer
  //
  m_timer.startClock();
  //
  for (is=0; is<m_sTrue.n(); is++) { // loop over all s_true
    m_pole->setTrueSignal( m_sTrue.getVal(is) );
    for (ie=0; ie<m_effTrue.n(); ie++) { // loop over eff true
      m_pole->setEffPdfMean( m_effTrue.getVal(ie) );
      for (ib=0; ib<m_bkgTrue.n(); ib++) { // loop over bkg true
        m_pole->setEffPdfMean( m_effTrue.getVal(ie) );
        m_pole->setBkgPdfMean( m_bkgTrue.getVal(ib) );
        //
        // Reset various counters
        //
	resetCoverage();   // for each s_true, reset coverage
	resetStatistics(); // dito, statistics
        //
        // For each 'truth', make nLoops pseudoexperiments and check if the truth is with the limits
        // 
        m_timer.startClock();
	for (j=0; j<m_nLoops; j++) {  // get stats
	  if (first) {
            m_timer.startTimer();       // used for estimating the time for the run
	    first = false;
	  }
          // generate pseudoexperiment using given ditributions of signal,bkg and eff.
	  m_pole->generatePseudoExperiment();
          // set the truth to the observed value since the construction done on the assumed truth
          // this will be reset after the loop
          m_pole->setEffPdfMean( m_pole->getEffObs() );
          m_pole->setBkgPdfMean( m_pole->getBkgObs() );
          m_pole->initAnalysis();
          // calculate the limit of the given experiment
	  if (!m_pole->analyseExperiment()) {
	    updateStatistics(false);          // statistics (only if activated)
	    if (nWarnings<maxWarnings) {
	      m_pole->printFailureMsg();
	      m_pole->getMeasurement().dump();
	      m_pole->printSetup();
	      std::cout << "s(true) = " << m_pole->getTrueSignal() << std::endl;
	      std::cout << "Pseudoexperiment will be ignored!" << std::endl;
	      if (nWarnings==maxWarnings) {
		std::cout << "WARNING: previous message will not be repeated." << std::endl;
	      }
	    }
	    nWarnings++;
	  } else {
            m_doneOneLoop = true;
	    updateCoverage();            // update the coverage
	    updateStatistics(true);          // statistics (only if activated)
	    if (!timingDone) {
	      nest++;
	      if (m_timer.checkTimer(5)) {
		timingDone = true;
                m_timer.stopTimer();
                m_timer.printEstimatedTime(nTotal, nest);
	      }
	    }
	  }
          m_pole->setEffPdfMean( m_effTrue.getVal(ie) );
          m_pole->setBkgPdfMean( m_bkgTrue.getVal(ib) );
	}
	calcCoverage();         // calculate coverage and its uncertainty
        m_timer.stopClock();
	outputCoverageResult(); // print the result
	calcStatistics();       // dito for the statistics...
	printStatistics();
	dumpExperiments();
      }
    }
  }
  double frate = (nTotal>0 ? double(nWarnings)/double(nTotal):0.0);
  std::cout << ">>>Limit calculation failure rate: " << frate << std::endl;
  if (frate>0.01) {
    std::cout << "WARNING: The failure rate in the limit calculations is large (>0.01)." << std::endl;
    std::cout << "         Possible cures:" << std::endl;
    std::cout << "         1. Increase N_belt  ( Pole::setBelt(N) )" << std::endl;
    std::cout << "         2. Increase hypothesis range ( Pole::setTestHyp() )" << std::endl;
    std::cout << "         3. Increase integration precision (Pole::setEffInt(),setBkgInt() )" << std::endl;
  }
  m_timer.printCurrentTime("\nEnd of run: ");
}

// void Coverage::doLoopOLD() {
//   if (m_pole==0) return;
//   //
//   // Number of loops to time for an estimate of the total time
//   //  int nest = (m_nLoops>100 ? 100:(m_nLoops>10 ? 10:1));
//   //
//   int nest = 0; // new - times a given period and not a given number of loops
//   //
//   int is,ie,ib,j;
//   bool first = true;
//   bool timingDone = false;
//   int nWarnings=0;
//   int nTotal=0;
//   const int maxWarnings=10;
//   //
//   m_pole->setCoverage(!m_collectStats); // make full limits when collecting statistics
//   //
//   startClock();
//   //
//   for (is=0; is<m_sTrue.n(); is++) { // loop over all s_true
//     m_sTrueMean = m_sTrue.getVal(is);
//     m_pole->setTrueSignal(m_sTrueMean);
//     for (ie=0; ie<m_effTrue.n(); ie++) { // loop over eff true
//       m_effMean = m_effTrue.getVal(ie);
//       for (ib=0; ib<m_bkgTrue.n(); ib++) { // loop over bkg true
// 	m_bkgMean = m_bkgTrue.getVal(ib);

// 	resetCoverage();   // for each s_true, reset coverage
// 	resetStatistics(); // dito, statistics
//     //
// 	startClock();
// 	for (j=0; j<m_nLoops; j++) {  // get stats
// 	  if (first) {
// 	    startTimer();       // used for estimating the time for the run
// 	    first = false;
// 	  }
// 	  nTotal++;
// 	  generateExperiment(); // generate pseudoexperiment using given ditributions of signal,bkg and eff.
// 	  m_pole->setNObserved(m_measNobs); // set values from generated experiment
// 	  m_pole->setEffMeas(m_measEff,m_effSigma,m_effDist);
// 	  m_pole->setBkgMeas(m_measBkg,m_bkgSigma,m_bkgDist);
// 	  m_pole->setEffBkgCorr(m_beCorr); // always the same...
// 	  m_pole->setEffInt();         // reset the integral ranges
// 	  m_pole->setBkgInt();
// 	  m_pole->setTestHyp(m_hypRMin,m_hypRMax,m_hypRStep);        // recalculate hypothesis range
// 	  //	  m_pole->initIntArrays(); // DONE in analyseExperiment()
// 	  //	  m_pole->initBeltArrays();
// 	  if (!m_pole->analyseExperiment()) { // calculate the limit of the given experiment
// 	    updateStatistics(false);          // statistics (only if activated)
// 	    if (nWarnings<maxWarnings) {
// 	      m_pole->printFailureMsg();
// 	      m_pole->getMeasurement().dump();
// 	      m_pole->printSetup();
// 	      std::cout << "s(true) = " << m_sTrueMean << std::endl;
// 	      std::cout << "Pseudoexperiment will be ignored!" << std::endl;
// 	      if (nWarnings==maxWarnings) {
// 		std::cout << "WARNING: previous message will not be repeated." << std::endl;
// 	      }
// 	    }
// 	    nWarnings++;
// 	  } else {
// 	    updateCoverage();            // update the coverage
// 	    updateStatistics(true);          // statistics (only if activated)
// 	    if (!timingDone) {
// 	      nest++;
// 	      if (checkTimer(5)) {
// 		timingDone = true;
// 		stopTimer();
// 		printEstimatedTime(nest);
// 	      }
// 	    }
// 	    //	  if ((!timingDone) && (j==nest-1)) { // calculate the estimated run time.
// 	    //	    stopTimer();
// 	    //	    printEstimatedTime(nest);
// 	    //	    timingDone = true;
// 	    //	  }
// 	  }
// 	}
// 	calcCoverage();         // calculate coverage and its uncertainty
// 	stopClock();
// 	outputCoverageResult(); // print the result
// 	calcStatistics();       // dito for the statistics...
// 	printStatistics();
// 	dumpExperiments();
//       }
//     }
//   }
//   double frate = (nTotal>0 ? double(nWarnings)/double(nTotal):0.0);
//   std::cout << ">>>Limit calculation failure rate: " << frate << std::endl;
//   if (frate>0.01) {
//     std::cout << "WARNING: The failure rate in the limit calculations is large (>0.01)." << std::endl;
//     std::cout << "         Possible cures:" << std::endl;
//     std::cout << "         1. Increase N_belt  ( Pole::setBelt(N) )" << std::endl;
//     std::cout << "         2. Increase hypothesis range ( Pole::setTestHyp() )" << std::endl;
//     std::cout << "         3. Increase integration precision (Pole::setEffInt(),setBkgInt() )" << std::endl;
//   }
//   endofRunTime();
//   //  printClockUsage(m_nLoops*m_sTrue.n()*m_effTrue.n()*m_bkgTrue.n());
// }

void Coverage::doExpTest() {
  //
  if (m_pole==0) return;
  int is,ie,ib,j;
  //
  for (is=0; is<m_sTrue.n(); is++) { // loop over all s_true
    m_pole->setTrueSignal( m_sTrue.getVal(is) );
    for (ie=0; ie<m_effTrue.n(); ie++) { // loop over eff true
      m_pole->setEffPdfMean( m_effTrue.getVal(ie) );
      for (ib=0; ib<m_bkgTrue.n(); ib++) { // loop over bkg true
        m_pole->setBkgPdfMean( m_bkgTrue.getVal(ib) );
	//
	resetStatistics(); // reset statistics
	//
	for (j=0; j<m_nLoops; j++) {  // get stats
	  m_pole->generatePseudoExperiment(); // generate pseudoexperiment using given ditributions of signal,bkg and eff.
	  pushMeas();           // save 'measured' data
	}
	calcStatistics();       // 
	printStatistics();
      }
    }
  }
}
