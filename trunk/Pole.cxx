//
// Issues:: Number of points in the integral construction -> check! setBkgInt
//
//
//
// This code contains classes for calculating the coverage of a CL
// belt algorithm.
//
#include <iostream>
#include <fstream>
#include <iomanip>

#include "Pole.h"

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

inline char *yesNo(bool var) {
  static char yes[4] = {'Y','e','s',char(0)};
  static char no[3] = {'N','o',char(0)};
  return (var ? &yes[0]:&no[0]);
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
/*!
  Main constructor.
 */
// HERE: Should have a truly empty with NO initiation of variables -> SPEED!
Pole::Pole() {
  m_verbose=0;
  m_coverage = false;
  m_nUppLim = 10;
  m_normMaxDiff = 0.01;
  m_lowerLimitNorm = -1.0;
  m_upperLimitNorm = -1.0;
  m_normInt = 1.0;
  m_normBkg = 1.0;
  m_normEff = 1.0;
  //
  m_sTrue = 1;
  //
  // Default observation
  //
  m_measurement.setNObserved(1);
  m_measurement.setEff(1.0,0.1,DIST_GAUS);
  m_measurement.setBkg(0.0,0.0,DIST_NONE);
  m_measurement.setBEcorr(0.0);
  //
  m_stepMin = 0.001;
  //
  m_dmus = 0.01;
  //
  m_effIntScale = 5.0;
  m_bkgIntScale = 5.0;
  //
  setEffInt(m_effIntScale,21);
  setBkgInt(m_bkgIntScale,21);
  //
  setTestHyp(0.01);
  //
  m_validInt  = false;
  m_nInt      = 0; // should be equal to the size of Int vectors below
  //
  m_validBestMu = false;
  m_nBelt       = 50;
  m_nBeltMinUsed = m_nBelt;
  m_nBeltMaxUsed = 0;
  m_suggestBelt = (m_nBelt<1);
  //
  // init list of suggested nBelt
  //
//   m_nBeltList.push_back(18); //0
//   m_nBeltList.push_back(20);
//   m_nBeltList.push_back(20);
//   m_nBeltList.push_back(22);
//   m_nBeltList.push_back(22);
//   m_nBeltList.push_back(23); //5
//   m_nBeltList.push_back(25);
//   m_nBeltList.push_back(32);
//   m_nBeltList.push_back(38);
//   m_nBeltList.push_back(38);
//   m_nBeltList.push_back(40); //10
  //
  m_useNLR = false;
}

Pole::~Pole() {
}

bool Pole::checkEffBkgDists() {
  bool change=false;
  // If ONLY one is DIST_GAUSCORR, make bkgDist == effDist
  if ( ((m_measurement.getEffDist() == DIST_GAUSCORR) && (m_measurement.getBkgDist() != DIST_GAUSCORR)) ||
       ((m_measurement.getEffDist() != DIST_GAUSCORR) && (m_measurement.getBkgDist() == DIST_GAUSCORR)) ) {
    m_measurement.setBkgDist(m_measurement.getEffDist());
    change = true;
  }
  if (m_measurement.getEffDist() != DIST_GAUSCORR) {
    m_measurement.setBEcorr(0);
    change = true;
  }
  //
  if (m_measurement.getEffDist()==DIST_NONE) {
    m_measurement.setEffSigma(0.0);
    change = true;
  }
  if (m_measurement.getBkgDist()==DIST_NONE) {
    m_measurement.setBkgSigma(0.0);
    change = true;
  }
  return change;
}
//
// Calculates the range for integration.
//
void Pole::setInt(double & low, double & high, double scale, double mean, double sigma, DISTYPE dist) {
  //
  double dx;
  double nsigma;
  double nmean;
  //
  if (dist==DIST_NONE) {
    low  = mean;
    high = mean;
  } else {
    switch (dist) {
    case DIST_GAUSCORR:
    case DIST_GAUS:
      low  = mean - scale*sigma;
      high = mean + scale*sigma;
      break;
    case DIST_FLAT:
      dx=sigma*1.73205081; // == sqrt(12)*0.5; ignore scale - always use full range
      low  = mean-dx;
      high = mean+dx;
      break;
    case DIST_LOGN:
      nmean  = m_gauss->getLNMean(mean,sigma);
      nsigma = m_gauss->getLNSigma(mean,sigma);
      low  = nmean - scale*nsigma;
      high = nmean + scale*nsigma;
      break;
    default: // ERROR STATE
      break;
    }
    if (dist!=DIST_LOGN) {
      if (low<0) { // lower limit below 0 - shift the range
	high-=low;
	low=0;
      }
    }
  }
}

void Pole::setEffInt(double scale,int n) {
  if (n<1) n = m_effRangeInt.n();
  if (n<2) n=2;
  if (scale<0) scale = m_effIntScale;
  m_validInt = false;
  m_validBestMu = false;
  //
  double low,high;
  if (isFullyCorrelated() && (m_measurement.getEffDist() == DIST_GAUSCORR)) { // Fully correlated - integrate only bkg with sigma/sqrt(2)
    low = m_measurement.getEffMeas();
    high = m_measurement.getEffMeas();
    n = 1;
  } else {
    setInt(low,high,scale,m_measurement.getEffMeas(),m_measurement.getEffSigma(),m_measurement.getEffDist());  
  }
  //
  m_effRangeInt.setRange(low,high,n);
}

// void Pole::setEffInt(double scale,double step) {
//   m_validInt = false;
//   m_validBestMu = false;
//   //
//   double low,high;
//   if (isFullyCorrelated() && (m_measurement.getEffDist() == DIST_GAUSCORR)) { // Fully correlated - integrate only bkg with sigma/sqrt(2)
//     low = m_measurement.getEffMeas();
//     high = m_measurement.getEffMeas();
//     step = 1.0;
//   } else {
//     if (step<=0) step = m_effRangeInt.step(); // use previous step size
//     if (scale<=0) { // use previous scale
//       scale = m_effIntScale;
//     } else { // set new scale
//       m_effIntScale = scale;
//     }
//     setInt(low,high,step,scale,m_measurement.getEffMeas(),m_measurement.getEffSigma(),m_measurement.getEffDist());
//   }
//   //
//   m_effRangeInt.setRange(low,high,step);
// }

void Pole::setBkgInt(double scale,int n) {
  if (n<1) n = m_bkgRangeInt.n();
  if (n<2) n=2;
  if (scale<0) scale = m_bkgIntScale;
  m_validInt = false;
  m_validBestMu = false;
  //
  double low,high;
  if (isFullyCorrelated() && (m_measurement.getBkgDist() == DIST_GAUSCORR)) { // Fully correlated - integrate only bkg with sigma/sqrt(2)
    low = m_measurement.getBkgMeas();
    high = m_measurement.getBkgMeas();
    n = 1;
  } else {
    setInt(low,high,scale,m_measurement.getBkgMeas(),m_measurement.getBkgSigma(),m_measurement.getBkgDist());  
  }
  //
  m_bkgRangeInt.setRange(low,high,n);
}
// void Pole::setBkgInt(double scale,double step) {
//   m_validInt = false;
//   m_validBestMu = false;
//   //
//   double low,high;
//   if (step<=0)  step  = m_bkgRangeInt.step();
//   if (scale<=0) {
//     scale = m_bkgIntScale;
//   } else {
//     m_bkgIntScale = scale;
//   }
//   setInt(low,high,step,scale,m_measurement.getBkgMeas(),m_measurement.getBkgSigma(),m_measurement.getBkgDist());
//   //
//   m_bkgRangeInt.setRange(low,high,step);
// }

void Pole::setTestHyp(double step) {
  //
  // Find hypothesis test range based on the input measurement
  // MUST be called after the measurement is set.
  //
  if (step<=0) {
    step = m_hypTest.step();
    if (step<=0) step = 0.01;
  }
  
  double low = BeltEstimator::getSigLow(getNObserved(),
					getEffDist(), getEffMeas(), getEffSigma(),
					getBkgDist(), getBkgMeas(), getBkgSigma(), m_normInt);
  double up  = BeltEstimator::getSigUp( getNObserved(),
					getEffDist(), getEffMeas(), getEffSigma(),
					getBkgDist(), getBkgMeas(), getBkgSigma(), m_normInt);
  m_hypTest.setRange(low,up,step);
}

void Pole::setTestHyp(double low, double high, double step) {
  //
  // Set explicitely the test range.
  // * step<=0  => step = (high-low)/1000
  // * high<low => call setTestHyp(step) [i.e, estimate the test range needed]
  //
  if (high<low) {
    setTestHyp(step);
  } else {
    if (step<=0) step = (high-low)/1000.0; // default 1000 pts
    m_hypTest.setRange(low,high,step);
  }
  
}
//
// MAYBE REMOVE THIS // HERE
bool Pole::checkParams() {
  bool rval=true;
  std::cout << "<<Pole::CheckParams() is disabled - Manually make sure that the integration limits are ok>>" << std::endl;
  return rval;
  /////////////////////////////////////////
  // check efficiency distribution
  // check background distribution
  // check true signal range
  // check mu_test range
  // check efficiency integration
  std::cout << "Checking efficiency integration - ";
  double dsLow, dsHigh;
  // remember, RangeInt for bkg and eff, LOGN, are in lnx, hence exp() is the true eff or bkg
  dsLow  = (m_measurement.getEffMeas() - m_effRangeInt.min())/m_measurement.getEffSigma();
  dsHigh = (m_effRangeInt.max()  - m_measurement.getEffMeas())/m_measurement.getEffSigma();
  if ( (m_measurement.getEffDist()!=DIST_NONE) && ( ((dsLow<4.0)&&(m_effRangeInt.min()>0)) ||
                                   (dsHigh<4.0) ) ) {
    std::cout << "Not OK" << std::endl;
    std::cout << "  Efficiency range for integration does not cover 4 sigma around the true efficiency." << std::endl;
    std::cout << "  Change range or efficiency distribution." << std::endl;
    rval = false;
  } else {
    std::cout << "OK" << std::endl;
  }
  // check background integration
  std::cout << "Checking background integration - ";
  dsLow  = (m_measurement.getBkgMeas() - getBkgIntMin())/m_measurement.getBkgSigma();
  dsHigh = (getBkgIntMax()  - m_measurement.getBkgMeas())/m_measurement.getBkgSigma();
  if ( (m_measurement.getBkgDist()!=DIST_NONE) && ( ((dsLow<4.0)&&(getBkgIntMin()>0)) ||
			    (dsHigh<4.0) ) ) {
    std::cout << "Not OK" << std::endl;
    std::cout << "  Background range for integration does not cover 4 sigma around the true background." << std::endl;
    std::cout << "  Change range or background distribution." << std::endl;
    rval = false;
  } else {
    std::cout << "OK" << std::endl;
  }
  return rval;
}

//

void Pole::initPoisson(int nlambda, int nn, double lmbmax) {
  //  if (m_poisson) m_poisson->init(nlambda, nn, lmbmax);
}

void Pole::initGauss(int ndata, double mumax) {
  //  if (m_gauss) m_gauss->init(ndata, mumax);
}

//
// Must be called after setEffInt, setBkgInt and setTestHyp
//
void Pole::initIntArrays() {
  m_nInt = m_effRangeInt.n()*m_bkgRangeInt.n();
  //  if (m_weightInt.size()<static_cast<unsigned int>(m_nInt)) {
  m_weightInt.resize(m_nInt,0.0);
  m_effInt.resize(m_nInt,0.0);
  m_bkgInt.resize(m_nInt,0.0);
  //  }
}

// int Pole::suggestBelt() {
//   int rval=50; // default value
//   int nbelt = static_cast<int>(m_nBeltList.size());
//   if (m_measurement.getNObserved()>=0) {
//     if (m_measurement.getNObserved()<nbelt) {
//       rval = m_nBeltList[m_measurement.getNObserved()];
//     } else {
//       rval = m_measurement.getNObserved()*4;
//     }
//   }
//   return rval;
// }


int Pole::suggestBelt() {
  int rval=50;
  if (m_measurement.getNObserved()>=0) {
    rval = BeltEstimator::getBeltMin(m_measurement.getNObserved(),
				     m_measurement.getEffDist(),
				     m_measurement.getEffMeas(), m_measurement.getEffSigma(),
				     m_measurement.getBkgDist(),
				     m_measurement.getBkgMeas(), m_measurement.getBkgSigma(),
				     m_normInt) + 5;
    if (rval<20) rval=20; // becomes less reliable for small n
  }
  if (m_verbose>1) std::cout << "Using max N(belt) = " << m_nBelt << std::endl;
  return rval;
}

void Pole::initBeltArrays() {
  if (m_suggestBelt) m_nBelt = suggestBelt();
  //
  // should not be initiated here
  //
  m_nBeltMinUsed = m_nBelt;
  m_nBeltMaxUsed = 0;

  //  unsigned int nbs = static_cast<unsigned int>(m_nBelt);
  //  if (m_muProb.size()!=nbs) {
  m_muProb.resize(m_nBelt,0.0);
  m_bestMuProb.resize(m_nBelt,0.0);
  m_bestMu.resize(m_nBelt,0.0);
  m_lhRatio.resize(m_nBelt,0.0);
    //  }
}

void Pole::initIntegral() {
  if (m_validInt) return;
  //
  static bool firstEffWarn=true;
  static bool firstBkgWarn=true;
  double effs;
  double bkgs;
  double eff_prob, bkg_prob, norm_prob;
  double sum_eff_prob=0;
  double sum_bkg_prob=0;
  unsigned int count=0;
  norm_prob = 1.0;
  bool full2d = false;
  /////////////////////////////////
  double sdetC=0, detC=0, seff=0, sbkg=0, vceff=0;
  double effMean,effSigma;
  double bkgMean,bkgSigma;
  if (m_measurement.getEffDist()==DIST_LOGN) {
    effMean  = m_gauss->getLNMean(m_measurement.getEffMeas(),m_measurement.getEffSigma());
    effSigma = m_gauss->getLNSigma(m_measurement.getEffMeas(),m_measurement.getEffSigma());
  } else {
    effMean  = m_measurement.getEffMeas();
    effSigma = m_measurement.getEffSigma();
  }
  if (m_measurement.getBkgDist()==DIST_LOGN) {
    bkgMean  = m_gauss->getLNMean(m_measurement.getBkgMeas(),m_measurement.getBkgSigma());
    bkgSigma = m_gauss->getLNSigma(m_measurement.getBkgMeas(),m_measurement.getBkgSigma());
  } else {
    bkgMean  = m_measurement.getBkgMeas();
    bkgSigma = m_measurement.getBkgSigma();
  }
  if (m_measurement.getEffDist()==DIST_GAUSCORR) {
    if (isFullyCorrelated() || isNotCorrelated()) {
      full2d=false;
    } else {
      full2d=true;
      detC = m_gauss->getDetC(m_measurement.getEffSigma(),m_measurement.getBkgSigma(),m_measurement.getBEcorr());
      sdetC = sqrt(detC);
      vceff = m_gauss->getVeffCorr(detC,m_measurement.getEffSigma(),m_measurement.getBkgSigma(),m_measurement.getBEcorr());
      seff  = sqrt(m_gauss->getVeff(detC,m_measurement.getBkgSigma())); // effective sigma for efficiency v1eff = detC/(s2*s2)
      sbkg  = sqrt(m_gauss->getVeff(detC,m_measurement.getEffSigma())); // dito for background            v2eff = detC/(s1*s1)
      if (m_verbose>2) {
	std::cout << "DetC = " << detC << std::endl;
	std::cout << "sDetC = " << sdetC << std::endl;
	std::cout << "Vceff = " << vceff << std::endl;
	std::cout << "seff = " << seff << std::endl;
	std::cout << "sbkg = " << sbkg << std::endl;
	std::cout << "corr = " << m_measurement.getBEcorr() << std::endl;
      }
    }
  }
  //////////////////////////////////
  //
  bool dosumbkg = false;
  double de=1.0;
  double db=1.0;
  if (m_effRangeInt.n() != 1) de=m_effRangeInt.step();
  if (m_bkgRangeInt.n() != 1) db=m_bkgRangeInt.step();
  double dedb  = de*db;
//   if ((m_effRangeInt.n() != 1) && (m_bkgRangeInt.n() != 1)) { // both eff and bkg vary
//     dedb = m_bkgRangeInt.step()*m_effRangeInt.step();
//   } else {
//     if ((m_effRangeInt.n() != 1) && (m_bkgRangeInt.n() == 1)) { // bkg const
//       dedb = m_effRangeInt.step();
//     } else {
//       if ((m_effRangeInt.n() == 1) && (m_bkgRangeInt.n() != 1)) {  // eff const
// 	dedb = m_bkgRangeInt.step();
//       }
//     }
//   }
  if (m_verbose>1) {
    std::cout << "InitMatrix: eff dist (mu,s)= " << m_measurement.getEffMeas() << "\t"
	      << m_measurement.getEffSigma() << std::endl;
    std::cout << "InitMatrix: bkg dist (mu,s)= " << m_measurement.getBkgMeas() << "\t"
	      << m_measurement.getBkgSigma() << std::endl;
    std::cout << "InitMatrix: dedb           = " << dedb << std::endl;
  }
  //
  for (int i=0;i<m_effRangeInt.n();i++) { // Loop over all efficiency points
    effs =  i*m_effRangeInt.step() + m_effRangeInt.min();
    if ( m_effRangeInt.n() == 1 ) {
      eff_prob = 1.0;
    } else {
      switch(m_measurement.getEffDist()) {
      case DIST_GAUS:
	eff_prob = m_gauss->getVal(effs,m_measurement.getEffMeas(),m_measurement.getEffSigma());
	//	std::cout << "GAUSSEFF: " << effs << " " << m_measurement.getEffMeas() << " " << m_measurement.getEffSigma() << " -> " << eff_prob << std::endl;
	break;
      case DIST_LOGN:
	effs = exp(effs);
	eff_prob = m_gauss->getValLogN(effs,effMean,effSigma)*effs;
	break;
      case DIST_GAUSCORR:
	if (isNotCorrelated()) {
	  eff_prob = m_gauss->getVal(effs,m_measurement.getEffMeas(),m_measurement.getEffSigma());
	} else {
	  eff_prob = 1.0; // will be set in the bkg loop
	}
	break;
      case DIST_FLAT:
	if (m_measurement.getEffSigma()>0) {
	  eff_prob = 1.0/(m_measurement.getEffSigma()*3.46410161513775);
	} else {
	  eff_prob = 1.0;
	}
	break;
      case DIST_NONE:
	eff_prob = 1.0;
	break;
      default:
	eff_prob = 1.0;
	break;
      }
    }
    sum_eff_prob += eff_prob;
    //
    for (int j=0;j<m_bkgRangeInt.n();j++) { // Loop over all background points
      bkgs =  j*m_bkgRangeInt.step() + m_bkgRangeInt.min();
      if ( m_bkgRangeInt.n() == 1 ) {
	bkg_prob = 1.0;
	dosumbkg = (i==0);
      } else {
	switch(m_measurement.getBkgDist()) {
	case DIST_GAUS:
	  bkg_prob = m_gauss->getVal(bkgs,m_measurement.getBkgMeas(),m_measurement.getBkgSigma());
	  dosumbkg = (i==0);
	  //	  std::cout << "GAUSSBKG: " << bkgs << " " << m_measurement.getBkgMeas() << " " << m_measurement.getBkgSigma() << " -> " << bkg_prob << std::endl;
	  break;
	case DIST_LOGN:
	  bkgs = exp(bkgs);
	  bkg_prob = m_gauss->getValLogN(bkgs,bkgMean,bkgSigma)*bkgs;
	  dosumbkg = (i==0);
	  break;
	case DIST_GAUSCORR:
	  if (isFullyCorrelated()) {
	    bkg_prob = m_gauss->getVal(bkgs,m_measurement.getBkgMeas(),m_measurement.getBkgSigma()/sqrt(2.0));
	    dosumbkg = (i==0);
	  } else {
	    if (isNotCorrelated()) {
	      bkg_prob = m_gauss->getVal(bkgs,m_measurement.getBkgMeas(),m_measurement.getBkgSigma());
	      dosumbkg = (i==0);
	    } else {
	      bkg_prob = m_gauss->getVal2D(effs,m_measurement.getEffMeas(),bkgs,m_measurement.getBkgMeas(),sdetC,seff,sbkg,vceff);
	      dosumbkg = true;
	    }
	  }
	  break;
	case DIST_FLAT:
	  if (m_measurement.getBkgSigma()>0) {
	    bkg_prob = 1.0/(m_measurement.getBkgSigma()*3.46410161513775);
	  } else {
	    bkg_prob = 1.0;
	  }
	  dosumbkg = (i==0);
	  break;
	case DIST_NONE:
	  bkg_prob = 1.0;
	  dosumbkg = (i==0);
	  break;
	default:
	  bkg_prob = 1.0;
	  dosumbkg = (i==0);
	  break;
	}
      }
      if (dosumbkg) sum_bkg_prob += bkg_prob; // only for first loop
      norm_prob = eff_prob*bkg_prob*dedb;
      //
      // Fill the arrays (integral)
      //
      if( (m_effRangeInt.n() == 1) && (m_bkgRangeInt.n() == 1) ) {
	m_bkgInt[count]    = m_measurement.getBkgMeas();
	m_effInt[count]    = m_measurement.getEffMeas();
	m_weightInt[count] = 1.0;
      } else {
	m_bkgInt[count]    = bkgs;
	m_effInt[count]    = effs;
	m_weightInt[count] = norm_prob;
      }
      //
      //      if (m_verbose>1) std::cout << "Matrix: (" << i << ", " <<  j << ") -> (" << m_effInt[count] << ", " << m_bkgInt[count] << ") " << m_weightInt[count] << "\n";
      //
      count++;
      if (count>m_nInt) {
	std::cout << "WARNING: Array index overflow. count = " << count
		  << ", n_e = " << m_effRangeInt.n()
		  << ", n_b = " << m_bkgRangeInt.n() << std::endl;
      }
    }
  }
  double norm_bkg = sum_bkg_prob*db;
  double norm_eff = sum_eff_prob*de;
  if (full2d) {
    norm_bkg *= de;  // if the distribution is a 2d gauss, the sum is over all area elements
    norm_eff  = 1.0; // per construction
  }
  double norm = norm_bkg*norm_eff;
  //
  for (unsigned int i=0; i<count; i++) {
    m_weightInt[i] = m_weightInt[i]/norm;
  }
  m_normEff = norm_eff;
  m_normBkg = norm_bkg;
  m_normInt = norm;
  if (m_verbose>2) {
    std::cout << "s_V  = " << getSVar() << std::endl;
    std::cout << "norm = " << m_normInt << std::endl;
  }
  if (firstEffWarn && (!full2d) && (!normOK(norm_eff))) {
    std::cout << std::endl;
    std::cout << "WARNING: normalisation of EFF is off: " << norm_eff << std::endl;
    std::cout << "         This might be OK if the intention is to use a bounded pdf." << std::endl;
    std::cout << "         The pdf is reweighted such that it is normalised." << std::endl;
    std::cout << "         Int. range = [" << m_effRangeInt.min() << ":" <<  m_effRangeInt.max() << "], de = " << de << std::endl;
    m_measurement.dump();
    std::cout << "--- MESSAGE IS NOT REPEATED! ---" << std::endl;
    std::cout << std::endl;
    

    firstEffWarn = false;
  }
  if (firstBkgWarn && (!normOK(norm_bkg))) {
    std::cout << std::endl;
    std::cout << "WARNING: normalisation of BKG is off: " << norm_bkg << std::endl;
    std::cout << "         This might be OK if the intention is to use a bounded pdf." << std::endl;
    std::cout << "         The pdf is reweighted such that it is normalised." << std::endl;
    std::cout << "         Int. range = [" << m_bkgRangeInt.min() << ":" <<  m_bkgRangeInt.max() << "]" << std::endl;
    m_measurement.dump();
    std::cout << "--- MESSAGE IS NOT REPEATED! ---" << std::endl;
    std::cout << std::endl;
    firstBkgWarn = false;
  }
  if (m_verbose>1) {
    std::cout << "First 10 from double integral (bkg,eff,w):" << std::endl;
    std::cout << "Eff: " << m_effRangeInt.n()   << "\t" << m_effRangeInt.step()
	      << "\t"    << m_effRangeInt.min() << "\t" << m_effRangeInt.max() << std::endl;
    std::cout << "Bkg: " << m_bkgRangeInt.n()   << "\t" << m_bkgRangeInt.step()
	      << "\t"    << m_bkgRangeInt.min() << "\t" << m_bkgRangeInt.max() << std::endl;
    for (int i=0; i<10; i++) {
      std::cout << m_bkgInt[i] << "\t" << m_effInt[i] << "\t"
		<< m_weightInt[i] << std::endl;
    }
  }
  m_validInt = true;
}



void Pole::findBestMu(int n) {
  // finds the best fit (mu=s) for a given n. Fills m_bestMu[n]
  // and m_bestMuProb[n].
  int i;
  double mu_test,lh_test,mu_best,mu_s_max,mu_s_min;
  //,dmu_s;
  double lh_max = 0.0;
  //
  mu_best = 0;
  if(n<m_measurement.getBkgMeas()) {
    m_bestMu[n] = 0; // best mu is 0
    m_bestMuProb[n] = calcProb(n,0);
  } else {
    //    mu_s_max = double(n)-m_measurement.getBkgMeas(); // OLD version
    mu_s_max = (double(n) - getBkgIntMin())/m_measurement.getEffMeas();
    if(mu_s_max<0) {mu_s_max = 0.0;}
    //    mu_s_min = (double(n) - m_bkgRangeInt.max())/m_effRangeInt.max();
    //    mu_s_min = mu_s_max/2.0; //
    mu_s_min = (double(n) - getBkgIntMax())/getEffIntMax();
    if(mu_s_min<0) {mu_s_min = 0.0;}
    //    dmu_s = 0.01; // HARDCODED:: Change!
    int ntst = 1+int((mu_s_max-mu_s_min)/m_dmus);
    //// TEMPORARY CODE - REMOVE /////
    //    ntst = 1000;
    //    m_dmus = (mu_s_max-mu_s_min)/double(ntst);
    //////////////////////////////////
    if (m_verbose>1) std::cout << "FindBestMu range: " << " I " << getBkgIntMax() << " " << getEffIntMax() << " "
			       << n << " " << m_measurement.getBkgMeas() << " " << ntst << " [" << mu_s_min << "," << mu_s_max << "] => ";
    for (i=0;i<ntst;i++) {
      mu_test = mu_s_min + i*m_dmus;
      lh_test = calcProb(n,mu_test);
      if(lh_test > lh_max) {
	lh_max = lh_test;
	mu_best = mu_test;
      }
    }  
    if (m_verbose>1) std::cout <<"s_best = " << mu_best << ", LH = " << lh_max << std::endl;
    m_bestMu[n] = mu_best; m_bestMuProb[n] = lh_max;  
  }
}

void Pole::findAllBestMu() {
  if (m_validBestMu) return;
// fills m_bestMuProb and m_bestMu (L(s_best + b)[n])
  for (int n=0; n<m_nBelt; n++) {
    findBestMu(n);
  }
  if (m_verbose>1) {
    std::cout << "First 10 from best fit (mean,prob):" << std::endl;
    std::cout << m_bestMu.size() << ":" << m_bestMuProb.size() << std::endl;
    for (int i=0; i<10; i++) {
      std::cout << m_bestMu[i] << "\t" << m_bestMuProb[i] << std::endl;
    }
  }
  m_validBestMu = true;
}

void Pole::calcLh(double s) {
  //  double norm_p=0.0;
  for (int n=0; n<m_nBelt; n++) {
    m_muProb[n] = calcProb(n, s);
    //    norm_p += m_muProb[n]; // needs to be renormalised - NO!! 
  }
}

double Pole::calcLhRatio(double s) {
  double norm_p = 0;
  if (m_useNLR) { // use method by Gary Hill
    double pbf;
    for (int n=0; n<m_nBelt; n++) {
      m_muProb[n] =  calcProb(n, s);
      if (n>m_measurement.getBkgMeas()) {
	pbf = static_cast<double>(n);
      } else {
	pbf = m_measurement.getBkgMeas();
      }
      pbf = m_poisson->getVal(n,pbf);
      m_lhRatio[n]  = m_muProb[n]/pbf;
      norm_p += m_muProb[n]; // needs to be renormalised
    }
  } else {
    for (int n=0; n<m_nBelt; n++) {
      m_muProb[n] =  calcProb(n, s);
      m_lhRatio[n]  = m_muProb[n]/m_bestMuProb[n];
      norm_p += m_muProb[n]; // needs to be renormalised
    }
  }
  std::cout << "BUG:: NOT YET OK (calcLhRatio())" << std::endl;
  return 0.0;
}

double Pole::calcLimit(double s) {
  int k,i;
  //
  double norm_p = 0;
  m_sumProb = 0;
  bool lowNfound=false;
  bool upNfound=false;
  int nInBelt=0;

  int nBeltMaxUsed = m_nBelt;
  int nBeltMinUsed = 0;
  const double minMuProb = 1e-10;
  //
  //  std::cout << "calcLimit for " << s << std::endl;
  if (m_useNLR) { // use method by Gary Hill
    double g,pbf;
    int n=0;
    while ((n<m_nBelt) && (!upNfound)) {
      //    for (int n=0; n<m_nBelt; n++) {
      if (!upNfound) {
	m_muProb[n] =  calcProb(n, s);
	if ((!lowNfound) && (m_muProb[n]>minMuProb)) {
	  lowNfound=true;
	  nBeltMinUsed = n;
	} else {
	  if ((nInBelt>1) && lowNfound && (m_muProb[n]<minMuProb)) {
	    upNfound = true;
	    nBeltMaxUsed = n-1;
	  }
	}
	if (lowNfound && (!upNfound)) nInBelt++;
      } else {
	m_muProb[n] =  0.0;
	m_lhRatio[n] = 0.0;
      }
      //      std::cout << "m_muProb[" << n << "] = " << m_muProb[n] << ", s = " << s << std::endl;
      if (!upNfound) {
	if (n>m_measurement.getBkgMeas()) {
	  g = static_cast<double>(n);
	} else {
	  g = m_measurement.getBkgMeas();
	}
	pbf = m_poisson->getVal(n,g);
	if (g==0) {
	  pbf=1.0;
	}
	m_lhRatio[n]  = m_muProb[n]/pbf;
	norm_p += m_muProb[n]; // check norm
      }
      n++;
    }
    if (m_verbose>2) {
      if ((norm_p>1.5) || (norm_p<0.5)) {
	std::cout << "Normalisation off (" << norm_p << ") for s= " << s << std::endl;
	for (int n=0; n<nBeltMaxUsed; n++) {
	  std::cout << "muProb[" << n << "] = " << m_muProb[n] << std::endl;
	}
      }
    }
  } else {
    int n=0;
    while ((n<m_nBelt) && (!upNfound)) {
      //    for (int n=0; n<m_nBelt; n++) {
      if (!upNfound) {
	m_muProb[n] =  calcProb(n, s);
	if ((!lowNfound) && (m_muProb[n]>minMuProb)) {
	  lowNfound=true;
	  nBeltMinUsed = n;
	} else {
	  if ((nInBelt>1) && lowNfound && (m_muProb[n]<minMuProb)) {
	    upNfound = true;
	    nBeltMaxUsed = n-1;
	  }
	}
	if (lowNfound && (!upNfound)) nInBelt++;
      } else {
	m_muProb[n] =  0.0;
	m_lhRatio[n] = 0.0;
      }
      if (!upNfound) {
	m_lhRatio[n]  = m_muProb[n]/m_bestMuProb[n];
	norm_p += m_muProb[n];
      }
      n++;
    }
  }
  if (norm_p>m_maxNorm) m_maxNorm=norm_p;
  //
  if (nBeltMinUsed<m_nBeltMinUsed) m_nBeltMinUsed = nBeltMinUsed;
  if (nBeltMaxUsed>m_nBeltMaxUsed) m_nBeltMaxUsed = nBeltMaxUsed;
  //
  //  if (m_verbose>1) std::cout << "Used max NBelt = " << m_nBeltMaxUsed << " ( " << m_nBelt << " )" << std::endl;
  k = m_measurement.getNObserved();

  if ((k>nBeltMaxUsed) || (k<=nBeltMinUsed)) m_lhRatio[k] = 0.0;

  if (k>=m_nBelt) {
    k=m_nBelt; // WARNING::
    std::cout << "WARNING:: n_observed is larger than the maximum n used for R(n,s)!!" << std::endl;
    std::cout << "          -> increase nbelt such that it is more than n_obs = " << m_measurement.getNObserved() << std::endl;
  }									\
  if (m_verbose>2) std::cout << "Got nBelt range: " << nBeltMinUsed << ":" << nBeltMaxUsed << "( max = " << m_nBelt << " )" << std::endl;
  // Calculate the probability for all n and the given s.
  // The Feldman-Cousins method dictates that for each n a
  // likelihood ratio (R) is calculated. The n's are ranked according
  // to this ratio. Values of n are included starting with that giving
  // the highest R and continuing with decreasing R until the total probability
  // matches the searched CL.
  // Below, the loop sums the probabilities for a given s and for all n with R>R0.
  // R0 is the likelihood ratio for n_observed.
  i=nBeltMinUsed;
  bool done=false;
  //  std::cout << "Norm_p = " << norm_p << std::endl;
  while (!done) {
    //  for(i=0;i<m_nBelt;i++) {
    //    m_muProb[i] = m_muProb[i]/norm_p; DO NOT NORMALISE - NOT NEEDED AS THEY ARE ALL 1
    //    for(k=0;k<m_nBelt;k++) {
    if(i != k) { 
      //      std::cout << "LHratio: s= " << s << "   i:k " << i << ":" << k << "    RL(i:k) = " << m_lhRatio[i] << ":" << m_lhRatio[k]
      //		<< "   prob = " << m_sumProb << std::endl;
      //    }
      if(m_lhRatio[i] > m_lhRatio[k])  {
	m_sumProb  +=  m_muProb[i];
      }
      if (m_verbose>9) {
	std::cout << "RL[" << i << "] = " << m_lhRatio[i] << ", RLmax[" << k << "] = " << m_lhRatio[k] << ", sumP = " << m_sumProb << std::endl;
      }
      //    }
    }
    i++;
    done = ((i>nBeltMaxUsed) || m_sumProb>m_cl); // CHANGE 11/8
  }
  //
  // Check if limit is reached.
  // For a given n_observed, we should have
  //   Sum(p) > cl for s < s_low
  //   Sum(p) < cl for s_low < s < s_upp
  //   Sum(p) > cl for s_upp < s
  //
  if (m_sumProb<m_cl) {
    if (m_foundLower) {
      m_upperLimit = s;
      m_upperLimitNorm = norm_p;
      //      m_foundUpper = true;
    } else {
      m_lowerLimit = s;
      m_lowerLimitNorm = norm_p;
      m_foundLower = true;
      m_foundUpper = false;
    }
  } else {
    if (m_foundLower) {
      m_foundUpper = true;
    }
  }
  return norm_p;
}

bool Pole::limitsOK() {
  bool rval=false;
  if (m_foundLower && m_foundUpper) {
    rval = (normOK(m_lowerLimitNorm) && normOK(m_upperLimitNorm));
  }
  return rval;
}

//*********************************************************************//
//*********************************************************************//
//*********************************************************************//

void Pole::calcConstruct(double s) {
  int i;
  //
  double norm_p = 0;
  double p;
  std::vector<double> muProb;
  std::vector<double> lhRatio;
  //
  //  std::cout << "calcLimit for " << s << std::endl;
  if (m_useNLR) { // use method by Gary Hill
    double pbf;
    for (int n=0; n<m_nBelt; n++) {
      p =  calcProb(n, s);
      if (n>m_measurement.getBkgMeas()) {
	pbf = static_cast<double>(n);
      } else {
	pbf = 0;
      }
      pbf = m_poisson->getVal(n,pbf);
      lhRatio.push_back(p/pbf);
      muProb.push_back(p);
      norm_p += p; // needs to be renormalised

    }
  } else {
    for (int n=0; n<m_nBelt; n++) {
      p =  calcProb(n, s);
      lhRatio.push_back(p/m_bestMuProb[n]);
      muProb.push_back(p);
      norm_p += p; // needs to be renormalised
    }
  }
  for (i=0; i<m_nBelt; i++) {
    muProb[i] = muProb[i]/norm_p;
  }
  //
  for (i=0; i<m_nBelt; i++) {
    std::cout << "CONSTRUCT: " << s << "\t" << i << "\t" << lhRatio[i] << "\t" << muProb[i] << std::endl;
  }
}

// NOTE: Simple sorting - should be put somewhere else...

void sort_index(std::vector<double> & input, std::vector<int> & index, bool reverse=false) {
  int ndata = input.size();
  if (ndata<=0) return;
  //
  int i;
  std::list< std::pair<double,int> > dl;
  //
  for (i=0; i<ndata; i++) {
    dl.push_back(std::pair<double,int>(input[i],i));
  }
  dl.sort();
  //
  if (!reverse) {
    std::list< std::pair<double,int> >::iterator dliter;
    for (dliter = dl.begin(); dliter != dl.end(); dliter++) {
      index.push_back(dliter->second);
    }
  } else {
    std::list< std::pair<double,int> >::reverse_iterator dliter;
    for (dliter = dl.rbegin(); dliter != dl.rend(); dliter++) {
      index.push_back(dliter->second);
    }
  }
}

double Pole::calcBelt(double s, int & n1, int & n2) {
  int i;
  //
  double norm_p = 0;
  double sumProb = 0;
  double p;
  std::vector<double> muProb;
  std::vector<double> lhRatio;
  std::vector<int> index;
  //
  //  std::cout << "calcLimit for " << s << std::endl;
  if (m_useNLR) { // use method by Gary Hill
    double pbf;
    for (int n=0; n<m_nBelt; n++) {
      p =  calcProb(n, s);
      if (n>m_measurement.getBkgMeas()) {
	pbf = static_cast<double>(n);
      } else {
	pbf = 0;
      }
      pbf = m_poisson->getVal(n,pbf);
      lhRatio.push_back(p/pbf);
      muProb.push_back(p);
      norm_p += p; // needs to be renormalised

    }
  } else {
    for (int n=0; n<m_nBelt; n++) {
      p =  calcProb(n, s);
      lhRatio.push_back(p/m_bestMuProb[n]);
      muProb.push_back(p);
      norm_p += p; // needs to be renormalised
    }
  }
  sort_index(lhRatio,index,true); // reverse sort
  //
  // Calculate the probability for all n and the given s.
  // The Feldman-Cousins method dictates that for each n a
  // likelihood ratio (R) is calculated. The n's are ranked according
  // to this ratio. Values of n are included starting with that giving
  // the highest R and continuing with decreasing R until the total probability
  // matches the searched CL.
  // Below, the loop sums the probabilities for a given s and for all n with R>R0.
  // R0 is the likelihood ratio for n_observed.
  for (i=0; i<m_nBelt; i++) {
    muProb[i] = muProb[i]/norm_p;
  }
  i=0;
  bool done=false;
  int nmin=-1;
  int nmax=-1;
  int n;

  while (!done) {
    n = index[i];
    p = muProb[i];
    sumProb +=p;
    if ((n<nmin)||(nmin<0)) nmin=n;
    if ((n>nmax)||(nmax<0)) nmax=n;
    //
    i++;
    done = ((i==m_nBelt) || sumProb>m_cl);
  }
  if ((nmin<0) || ((nmin==0)&&(nmax==0))) {
    nmin=0;
    nmax=1;
    sumProb=1.0;
  }
  n1 = nmin;
  n2 = nmax;
  std::cout << "CONFBELT: " << s << "\t" << n1 << "\t" << n2 << "\t" << sumProb << "\t"
	    << lhRatio[n1] << "\t" << lhRatio[n2] << "\t" << norm_p << "\t"
	    << index[0] << "\t" << lhRatio[index[0]] << std::endl;
  return sumProb;
}

//*********************************************************************//
//*********************************************************************//
//*********************************************************************//

void Pole::findConstruct() {
  double mu_test;
  int i = 0;
  bool done = (i==m_hypTest.n());
  //
  while (!done) {
    mu_test = m_hypTest.min() + i*m_hypTest.step();
    calcConstruct(mu_test);
    i++;
    done = (i==m_hypTest.n()); // must loop over all hypothesis
  }
}

void Pole::findBelt() {
  double mu_test;
  int i = 0;
  bool done = (i==m_hypTest.n());
  //
  int n1,n2;
  while (!done) {
    mu_test = m_hypTest.min() + i*m_hypTest.step();
    calcBelt(mu_test,n1,n2);
    i++;
    done = (i==m_hypTest.n()); // must loop over all hypothesis
  }
}

bool Pole::findLimits() {
  //  static bool heavyDBG=false;
  m_maxNorm = -1.0;
  m_foundLower = false;
  m_foundUpper = false;
  m_lowerLimit = 0;
  m_upperLimit = 0;
  m_lowerLimitNorm = 0;
  m_upperLimitNorm = 0;
  //
  if (m_measurement.getNObserved()>=m_nBelt) return false;
  double mu_test;
  int i = 0;
  bool done = (i==m_hypTest.n());
  bool nuppOK = false;
  int nupp=0;
  //
  //
  double p;
  //  int n1,n2;
  bool firstLower=true;
  int nAfterFirstLow=0;
  int nInBelt=0;
  while (!done) {
    mu_test = m_hypTest.min() + i*m_hypTest.step();
//     if (heavyDBG) {
//       savedVerb = m_verbose;
//       m_verbose=10;
//     }
    p=calcLimit(mu_test);
    //    calcBelt(mu_test,n1,n2);
    if (m_verbose>2) std::cout << "findLimits(): " << mu_test << " norm = " << p << " sump = " << m_sumProb << std::endl;
    if (m_foundLower) {
      if (firstLower) {
	firstLower = false;
      }
      if (!m_foundUpper) nAfterFirstLow++;
    }
    if (m_foundUpper) {
      if (nAfterFirstLow==5) {
	m_foundLower=false;
	m_lowerLimit=0;
	m_lowerLimitNorm=0;
	m_foundUpper=false;
	m_upperLimit=0;
	m_upperLimitNorm=0;
	//
	nAfterFirstLow=0;
	firstLower = true;
      }
    }
    i++;
    if (m_foundUpper) { // sample a few afterwards just to make sure we get the last one
//       nupp++;
//       nuppOK = ((nupp>=m_nUppLim) || (m_nUppLim<1));
      nuppOK = true; // SKIP THE SCAN...
    }
    done = ((i==m_hypTest.n()) || // Done if last hypothesis reached
	    (nuppOK)           || //
	    (!normOK(p)));        // normalisation not OK
  }
  if (m_verbose>1) {
    if (limitsOK()) {
      // && (!heavyDBG)) {
      std::cout << "LIMITS(N,e,b,l,u): ";
      coutFixed(4,m_measurement.getNObserved());  std::cout << "\t";
      coutFixed(4,m_measurement.getEffMeas());    std::cout << "\t";
      coutFixed(4,m_measurement.getBkgMeas());    std::cout << "\t";
      coutFixed(4,m_lowerLimit); std::cout << "\t";
      coutFixed(4,m_upperLimit); std::cout << std::endl;
      //      heavyDBG = false;
      //      findLimits();
      //      heavyDBG = false;
    } else {
      std::cout << "Limits NOT OK!" << std::endl;
    }
  }
  if (limitsOK()) {
    if (m_verbose>1) {
      if (m_upperLimit-m_lowerLimit<1.0) { // DEBUG! Sometimes limits are accepted for large s when the range is too small
	std::cout << "\n<************* WARNING **************>" << std::endl;
	std::cout << "n points in belt = " << nAfterFirstLow << std::endl;
	printFailureMsg();
	printSetup();
	i=0;
	while (i<m_hypTest.n()) {
	  mu_test = m_hypTest.min() + i*m_hypTest.step();
	  p=calcLimit(mu_test);
	  std::cout << "findLimitsDEBUG(): " << mu_test << " norm = " << p << " sump = " << m_sumProb << std::endl;
	  std::cout << std::endl;
	  i++;
	}
      }
    }
  }
  return limitsOK();
}

bool Pole::findCoverageLimits() {
  //
  // Same as findLimits() but it will stop the scan as soon it is clear
  // that the true s (m_sTrueMean) is inside or outside the range.
  //
  m_maxNorm = -1.0;
  m_foundLower = false;
  m_foundUpper = false;
  m_lowerLimit = 0;
  m_upperLimit = 0;
  //
  if (m_measurement.getNObserved()>=m_nBelt) return false;
  //
  double mu_test;
  int i = 0;
  bool done = (i==m_hypTest.n());
  bool decided = false;
  bool nuppOK = false;
  int  nupp=0;
  //
  double p;
  while (!done) {
    mu_test = m_hypTest.min() + i*m_hypTest.step();
    p=calcLimit(mu_test);
    if (m_verbose>2) std::cout << "findCoverageLimits():"
			       << " low,up = " << m_foundLower << ":" << m_foundUpper
			       << " s_true = " << m_sTrue
			       << " s_test = " << mu_test
			       << " norm = " << p
			       << " sump = " << m_sumProb << std::endl;
    i++;
    if (mu_test>m_sTrue) {
      //    if ((!m_foundLower) && (mu_test>m_sTrue)) { // for sure outside
      if (!m_foundLower) { // outside, below
	decided = true;
	m_lowerLimit = mu_test;     // create a fake set of limits
	m_upperLimit = mu_test+0.1;
	m_lowerLimitNorm = p; // not really needed
	m_upperLimitNorm = p; // not exactly true - only the lower limit norm is relevant
	m_foundLower = true;
	m_foundUpper = true;
      } else if (m_foundLower && (!m_foundUpper)) { // inside
	decided = true;
	m_upperLimit = mu_test+0.1;
      }
    } else if (m_foundUpper) {
      if (m_upperLimit>m_sTrue) { // for sure inside
	decided = true; // ...skip limitsOK()
      } else { // m_sTrue is above - wait a bit before deciding
	// SKIP THIS WAITING AS IT CAUSES PROBLEMS!
// 	nupp++;
// 	nuppOK = ((nupp>=m_nUppLim) || (m_nUppLim<1));
// 	if (nuppOK) decided = true;
	decided = true;
      }
    }
    done = ((decided) ||
	    (i==m_hypTest.n()));
  }
  if (m_verbose>1) {
    if (limitsOK()) {
      std::cout << "COVLIMITS(N,s,e,b,l,u): ";
      coutFixed(4,m_measurement.getNObserved());  std::cout << "\t";
      coutFixed(4,m_sTrue);      std::cout << "\t";
      coutFixed(4,m_measurement.getEffMeas());    std::cout << "\t";
      coutFixed(4,m_measurement.getBkgMeas());    std::cout << "\t";
      coutFixed(4,m_lowerLimit); std::cout << "\t";
      coutFixed(4,m_upperLimit); std::cout << std::endl;;
    } else {
      std::cout << "Limits NOT OK!" << std::endl;
    }
  }
  return decided;
}

void Pole::initAnalysis() {
  if (m_verbose>0) std::cout << "Initialise arrays" << std::endl;
  initIntArrays();
  initBeltArrays();
  if (m_verbose>0) std::cout << "Constructing integral" << std::endl;
  initIntegral();  // generates predefined parts of double integral ( eq.7)
}

bool Pole::analyseExperiment() {
  bool rval=false;
  initAnalysis();
  if (!m_useNLR) {
    if (m_verbose>0) std::cout << "Finding s_best" << std::endl;
    findAllBestMu(); // loops
  }
  if (m_verbose>3 && m_verbose<10) {
    findBelt();
  }
  if (m_verbose>0) std::cout << "Calculating limit" << std::endl;
  if (m_coverage) {
    rval=findCoverageLimits();
  } else {
    rval=findLimits();
  }
  // Should not do this - if probability is OK then the belt is also OK...?
  // The max N(Belt) is defined by a cutoff in probability (very small)
  //  if (m_nBeltMaxUsed==m_nBelt) rval=false; // reject limit if the full belt is used
  return rval;
}

void Pole::printLimit(bool doTitle) {
  if (doTitle) {
    std::cout << "-------------------------------------------------" << std::endl;
    std::cout << " Max N(belt) : " << m_nBeltMaxUsed << std::endl;
    std::cout << "-------------------------------------------------" << std::endl;
    std::cout << " Nobs  \t  Eff   \t Bkg" << std::endl;
    std::cout << "-------------------------------------------------" << std::endl;
  }
  coutFixed(4,m_measurement.getNObserved()); std::cout << "\t";
  coutFixed(6,m_measurement.getEffMeas()); std::cout << "\t";
  coutFixed(6,m_measurement.getEffSigma()); std::cout << "\t";
  coutFixed(6,m_measurement.getBkgMeas()); std::cout << "\t";
  coutFixed(6,m_measurement.getBkgSigma()); std::cout << "\t";
  std::cout << "[ ";
  coutFixed(2,m_lowerLimit); std::cout << ", ";
  coutFixed(2,m_upperLimit); std::cout << " ]";
  std::cout << std::endl;
}

void Pole::printSetup() {
  std::cout << "\n";
  std::cout << "================ P O L E ==================\n";
  std::cout << " Confidence level   : " << m_cl << std::endl;
  std::cout << " N observed         : " << m_measurement.getNObserved() << std::endl;
  std::cout << "----------------------------------------------\n";
  std::cout << " Coverage friendly  : " << yesNo(m_coverage) << std::endl;
  std::cout << " True signal        : " << m_sTrue << std::endl;
  std::cout << "----------------------------------------------\n";
  std::cout << " Efficiency meas    : " << m_measurement.getEffMeas() << std::endl;
  std::cout << " Efficiency sigma   : " << m_measurement.getEffSigma() << std::endl;
  std::cout << " Efficiency dist    : " << distTypeStr(m_measurement.getEffDist()) << std::endl;
  std::cout << "----------------------------------------------\n";
  std::cout << " Background meas    : " << m_measurement.getBkgMeas() << std::endl;
  std::cout << " Background sigma   : " << m_measurement.getBkgSigma() << std::endl;
  std::cout << " Background dist    : " << distTypeStr(m_measurement.getBkgDist()) << std::endl;
  std::cout << "----------------------------------------------\n";
  std::cout << " Bkg-Eff correlation: " << m_measurement.getBEcorr() << std::endl;
  std::cout << "----------------------------------------------\n";
  std::cout << " Int. eff. min      : " << m_effRangeInt.min() << std::endl;
  std::cout << " Int. eff. max      : " << m_effRangeInt.max() << std::endl;
  std::cout << " Int. eff. N pts    : " << m_effRangeInt.n() << std::endl;
  std::cout << "----------------------------------------------\n";
  std::cout << " Int. bkg. min      : " << m_bkgRangeInt.min() << std::endl;
  std::cout << " Int. bkg. max      : " << m_bkgRangeInt.max() << std::endl;
  std::cout << " Int. bkg. N pts    : " << m_bkgRangeInt.n() << std::endl;
  std::cout << "----------------------------------------------\n";
  std::cout << " Test hyp. min      : " << m_hypTest.min() << std::endl;
  std::cout << " Test hyp. max      : " << m_hypTest.max() << std::endl;
  std::cout << " Test hyp. step     : " << m_hypTest.step() << std::endl;
  std::cout << "----------------------------------------------\n";
  if (m_suggestBelt) {
    std::cout << " Belt N max         : variable" << std::endl;
  } else {
    std::cout << " Belt N max         : " << m_nBelt << std::endl;
  }
  std::cout << " Step mu_best       : " << m_dmus << std::endl;
  std::cout << "----------------------------------------------\n";
  std::cout << " Use NLR            : " << yesNo(m_useNLR) << std::endl;
  std::cout << "----------------------------------------------\n";
  std::cout << " Verbosity          : " << m_verbose << std::endl;
  std::cout << "==============================================\n";
  //
}

void Pole::printFailureMsg() {
  std::cout << "ERROR: limit calculation failed. Possible causes:" << std::endl;
  std::cout << "1. nbelt is too small (set nbelt = " << getNBelt() << ", max used = " << getNBeltMaxUsed() << ")" << std::endl;
  std::cout << "2. precision in integrations (eff,bkg) not sufficient" << std::endl;
  std::cout << "3. hypethesis test range too small ( max = " << m_hypTest.max() << " )" << std::endl;
  std::cout << "4. if Poisson table is used, insufficient precision. Symptom: probability of lower/upper limit diverges from unity." << std::endl;
  std::cout << "Input:" << std::endl;
  std::cout << "   N(obs)     = " << getNObserved() << std::endl;
  std::cout << "Results:" << std::endl;
  std::cout << "   probability (should be = CL)  = " << getSumProb() << std::endl;
  std::cout << "   lower lim norm (should be 1)  = " << getLowerLimitNorm() << std::endl;
  std::cout << "   upper lim norm (ditto)        = " << getUpperLimitNorm() << std::endl;
  std::cout << "   lower lim                     = " << getLowerLimit() << std::endl;
  std::cout << "   upper lim                     = " << getUpperLimit() << std::endl;
}
