#include <iomanip>
#include <algorithm>
#include <iterator>
#include "Pole.h"
#include "Power.h"
#include "PseudoExperiment.h"


Power::Power() {
  m_pole = 0;
  m_power = 0;
  m_sHyp = 0.0;
  m_nLoops = 1;
}

Power::Power(Pole *pole) {
  m_pole = pole;
  m_power = 0;
  m_sHyp = 0.0;
  m_nLoops = 1;
  if (pole) m_experiment.setMeasurement(pole->getMeasurement());
}


Power::~Power() {
}

bool Power::calculate(double strue) {
  if (m_pole==0) return false;
  //
//   if (strue<m_sHyp) {
//     s0 = strue;
//     s1 = m_sHyp;
//   } else {
//     s1 = strue;
//     s0 = m_sHyp;
//   }
  //  m_pole->setTestHyp();
  //  m_pole->initAnalysis();
  if (!m_pole->usesNLR()) {
    if (m_verbose>1) std::cout << "Finding all best mu" << std::endl;
    m_pole->findAllBestMu();
  }
  //  const Range<double> *hyprng = m_pole->getHypTest();
  //  std::cout << "Hyp test range min/max: " << hyprng->min() << " : " << hyprng->max() << std::endl;
  //  std::cout << "H0 = " << m_sHyp << " , H1 = " << strue << std::endl;
  if (m_verbose>1) std::cout << "Calculating belt" << std::endl;
  m_probHyp  = m_pole->calcBelt(m_sHyp,m_n1hyp, m_n2hyp,false,-1.0);
  m_probTrue = m_pole->calcBelt(strue,m_n1true,m_n2true,false,-1.0); // muProb is now for H1 (true)
//   std::cout << std::endl;
//   std::cout << "Nobs = " << m_pole->getNObserved() << std::endl;
//   std::cout << "H0: s = " << m_sHyp  << " , [n1,n2] = " << m_n1hyp  << " , " << m_n2hyp  << std::endl;
//   std::cout << "H1: s = " << strue   << " , [n1,n2] = " << m_n1true << " , " << m_n2true << std::endl;
  //
  return true;
}

void Power::resetPower() {
  m_nOutside = 0;
  m_nTotal = 0;
  m_power = 0;
  m_sumPOutside = 0;
  m_sumP = 0;
}

void Power::updatePower() {
  int nobs = m_pole->getNObserved();
  double ph1 = m_pole->getMuProb(nobs);
  if ((nobs<m_n1hyp) || (nobs>m_n2hyp)) {
    m_nOutside++;
    m_sumPOutside += ph1;
  }
  m_nTotal++;
  m_sumP += ph1;
}

void Power::outputResult() {
  std::cout << "POWER:\t" << m_sHyp << "\t" << m_pole->getSTrue() << "\t" << m_power << "\t" << m_powerUnc << "\t" << m_sumP
    //	    << "\t" << m_sumPOutside << "\t" << m_sumP
    //	    << "\t" << m_nOutside << "\t" << m_nTotal
	    << std::endl;
}

void Power::doLoop() {
  if (m_pole==0) return;
  //
  // Number of loops to time for an estimate of the total time
  //  int nest = (m_nLoops>100 ? 100:(m_nLoops>10 ? 10:1));
  int nest = 0; // new - times a given period and not a given number of loops
  //
  int is,ie,ib,j;
  bool first = true;
  bool timingDone = false;
  int nWarnings=0;
  int nTotal=0;
  const int maxWarnings=10;
  double smin,smax,strue;
  //
  //  startClock();
  //
  //
  for (is=0; is<m_sRange.n(); is++) {
    strue = m_sRange.getVal(is);
    if (m_verbose>1) std::cout << "POWERLOOP: S(true) = " << strue << std::endl;
    if (strue>m_sHyp) {
      smin = m_sHyp;
      smax = strue;
    } else {
      smax = m_sHyp;
      smin = strue;
    }
    m_pole->setTrueSignal(strue);
    m_pole->setTestHyp(smin,smax,m_sRange.step());        // hyp range [s0,...s1]
    m_experiment.setTrueSignal(strue);
    //
    resetPower();   // for each s_true, reset power
    //
    //startClock();
    for (j=0; j<m_nLoops; j++) {  // get stats
//       if (first) {
// 		startTimer();       // used for estimating the time for the run
// 	first = false;
//       }
      nTotal++;
      if (m_verbose>1) std::cout << "Generating experiment " << j << std::endl;
      m_experiment.generateMeasurement(m_measurement); // generate pseudoexperiment using given ditributions of signal,bkg and eff.
      
      if (m_verbose>1) std::cout << "Set measurement+init arrays" << std::endl;
      m_pole->setMeasurement(m_measurement);
      m_pole->setEffInt();         // reset the integral ranges
      m_pole->setBkgInt();
      m_pole->initIntArrays();
      m_pole->initBeltArrays();
      if (m_verbose>1) std::cout << "Init integral" << std::endl;
      m_pole->initIntegral();
      if (m_verbose>1) std::cout << "Calculate..." << std::endl;
      if (!calculate(strue)) {
	if (nWarnings<maxWarnings) {
	  m_pole->printFailureMsg();
	  m_pole->getMeasurement().dump();
	  m_pole->printSetup();
	  std::cout << "s(true) = " << strue << std::endl;
	  std::cout << "Pseudoexperiment will be ignored!" << std::endl;
	  if (nWarnings==maxWarnings) {
	    std::cout << "WARNING: previous message will not be repeated." << std::endl;
	  }
	}
	nWarnings++;
      } else {
	updatePower();            // update the coverage
// 	if (!timingDone) {
// 	  nest++;
// 	  if (checkTimer(5)) {
// 	    timingDone = true;
// 	    stopTimer();
// 	    printEstimatedTime(nest);
// 	  }
// 	}
      }
    }
    calcPower();         // calculate coverage and its uncertainty
    //    stopClock();
    outputResult(); // print the result
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
  //  endofRunTime();
}
