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
  m_method = RL_FHC2;
  m_nUppLim = 10;
  m_normMaxDiff = 0.01;
  m_lowerLimitNorm = -1.0;
  m_upperLimitNorm = -1.0;
  //
  m_sTrue = 1;
  //
  m_poisson = &PDF::gPoisTab;
  m_gauss   = &PDF::gGauss;
  m_gauss2d = &PDF::gGauss2D;
  m_logNorm = &PDF::gLogNormal;
  //
  // Default observation
  //
  m_measurement.setNObserved(1);
  m_measurement.setEffPdf(1.0,0.1,PDF::DIST_GAUS);
  m_measurement.setEffObs(); // sets the mean from PDF as the observed efficiency
  m_measurement.setBkgPdf(0.0,0.0,PDF::DIST_NONE);
  m_measurement.setBkgObs();
  m_measurement.setBEcorr(0.0);
  //
  m_stepMin = 0.001;
  m_minMuProb = 1e-6;
  //
  m_dmus = 0.01;
  m_nmusMax = 100;
  //
  m_measurement.setEffInt(5.0,21);
  m_measurement.setBkgInt(5.0,21);
  //
  setTestHyp(0.01);
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
}

Pole::~Pole() {
}

bool Pole::checkEffBkgDists() {
  bool change=false;
  // If ONLY one is PDF::DIST_GAUS2D, make bkgDist == effDist
  if ( ((m_measurement.getEffPdfDist() == PDF::DIST_GAUS2D) && (m_measurement.getBkgPdfDist() != PDF::DIST_GAUS2D)) ||
       ((m_measurement.getEffPdfDist() != PDF::DIST_GAUS2D) && (m_measurement.getBkgPdfDist() == PDF::DIST_GAUS2D)) ) {
    m_measurement.setBkgPdf(m_measurement.getEffPdfMean(),m_measurement.getEffPdfSigma(),m_measurement.getEffPdfDist());
    change = true;
  }
  if (m_measurement.getEffPdfDist() != PDF::DIST_GAUS2D) {
    m_measurement.setBEcorr(0);
    change = true;
  }
  //
  if (m_measurement.getEffPdfDist()==PDF::DIST_NONE) {
    m_measurement.setEffPdfSigma(0.0);
    change = true;
  }
  if (m_measurement.getBkgPdfDist()==PDF::DIST_NONE) {
    m_measurement.setBkgPdfSigma(0.0);
    change = true;
  }
  return change;
}

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
					getBkgDist(), getBkgMeas(), getBkgSigma(), m_measurement.getNuisanceIntNorm());
  double up  = BeltEstimator::getSigUp( getNObserved(),
					getEffDist(), getEffMeas(), getEffSigma(),
					getBkgDist(), getBkgMeas(), getBkgSigma(), m_measurement.getNuisanceIntNorm());
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
// bool Pole::checkParams() {
//   bool rval=true;
//   std::cout << "<<Pole::CheckParams() is disabled - Manually make sure that the integration limits are ok>>" << std::endl;
//   return rval;
//   /////////////////////////////////////////
//   // check efficiency distribution
//   // check background distribution
//   // check true signal range
//   // check mu_test range
//   // check efficiency integration
//   std::cout << "Checking efficiency integration - ";
//   double dsLow, dsHigh;
//   // remember, RangeInt for bkg and eff, LOGN, are in lnx, hence exp() is the true eff or bkg
//   dsLow  = (m_measurement.getEffObs() - m_effRangeInt.min())/m_measurement.getEffPdfSigma();
//   dsHigh = (m_effRangeInt.max()  - m_measurement.getEffObs())/m_measurement.getEffPdfSigma();
//   if ( (m_measurement.getEffPdfDist()!=PDF::DIST_NONE) && ( ((dsLow<4.0)&&(m_effRangeInt.min()>0)) ||
//                                    (dsHigh<4.0) ) ) {
//     std::cout << "Not OK" << std::endl;
//     std::cout << "  Efficiency range for integration does not cover 4 sigma around the true efficiency." << std::endl;
//     std::cout << "  Change range or efficiency distribution." << std::endl;
//     rval = false;
//   } else {
//     std::cout << "OK" << std::endl;
//   }
//   // check background integration
//   std::cout << "Checking background integration - ";
//   dsLow  = (m_measurement.getBkgObs() - getBkgIntMin())/m_measurement.getBkgPdfSigma();
//   dsHigh = (getBkgIntMax()  - m_measurement.getBkgObs())/m_measurement.getBkgPdfSigma();
//   if ( (m_measurement.getBkgPdfDist()!=PDF::DIST_NONE) && ( ((dsLow<4.0)&&(getBkgIntMin()>0)) ||
// 			    (dsHigh<4.0) ) ) {
//     std::cout << "Not OK" << std::endl;
//     std::cout << "  Background range for integration does not cover 4 sigma around the true background." << std::endl;
//     std::cout << "  Change range or background distribution." << std::endl;
//     rval = false;
//   } else {
//     std::cout << "OK" << std::endl;
//   }
//   return rval;
// }

//
// OBSOLETE - TODO: Remove
// void Pole::initPoisson(int nlambda, int nn, double lmbmax) {
//   //  if (m_poisson) m_poisson->init(nlambda, nn, lmbmax);
// }

// void Pole::initGauss(int ndata, double mumax) {
//   //  if (m_gauss) m_gauss->init(ndata, mumax);
// }


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
				     m_measurement.getEffPdfDist(),
				     m_measurement.getEffObs(), m_measurement.getEffPdfSigma(),
				     m_measurement.getBkgPdfDist(),
				     m_measurement.getBkgObs(), m_measurement.getBkgPdfSigma(),
				     m_measurement.getNuisanceIntNorm());
    if (rval<20) rval=20; // becomes less reliable for small n
  }
  if (m_verbose>1) {
    std::cout << "Using max N(belt) = " << rval << std::endl;
    std::cout << "t(test var) = " << getSVar() << std::endl;
  }
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


void Pole::findBestMu(int n) {
  // finds the best fit (mu=s) for a given n. Fills m_bestMu[n]
  // and m_bestMuProb[n].
  int i;
  double mu_test,lh_test,mu_best,sMax,sMin;
  //,dmu_s;
  double lh_max = 0.0;
  //
  mu_best = 0;
  if(n<m_measurement.getBkgObs()) {
    m_bestMu[n] = 0; // best mu is 0
    m_bestMuProb[n] = m_measurement.calcProb(n,0);
  } else {
    //    sMax = double(n)-m_measurement.getBkgObs(); // OLD version
    sMax = (double(n) - getBkgIntMin())/m_measurement.getEffObs();
    if(sMax<0) {sMax = 0.0;}
    //    sMin = (double(n) - m_bkgRangeInt.max())/m_effRangeInt.max();
    //    sMin = sMax/2.0; //
    sMin = (double(n) - getBkgIntMax())/getEffIntMax();
    if(sMin<0) {sMin = 0.0;}
    sMin = (sMax-sMin)*0.6 + sMin;
    //    dmu_s = 0.01; // HARDCODED:: Change!
    int ntst = 1+int((sMax-sMin)/m_dmus);
    double dmus=m_dmus;
    if (ntst>m_nmusMax) {
      ntst=m_nmusMax;
      dmus=(sMax-sMin)/double(ntst-1);
    }
    //// TEMPORARY CODE - REMOVE /////
    //    ntst = 1000;
    //    m_dmus = (sMax-sMin)/double(ntst);
    //////////////////////////////////
    if (m_verbose>1) std::cout << "FindBestMu range: " << " I " << getBkgIntMax() << " " << getEffIntMax() << " "
			       << n << " " << m_measurement.getBkgObs() << " " << ntst << " [" << sMin << "," << sMax << "] => ";
    int imax=-10;
    for (i=0;i<ntst;i++) {
      mu_test = sMin + i*dmus;
      lh_test = m_measurement.calcProb(n,mu_test);
      if(lh_test > lh_max) {
	imax = i;
	lh_max = lh_test;
	mu_best = mu_test;
      }
    }
    bool lowEdge = (imax==0)&&(ntst>0);
    bool upEdge  = (imax==ntst-1);
    if (lowEdge) { // lower edge found - just to make sure, scan downwards a few points
      int nnew = ntst/10;
      i=1;
      while ((i<=nnew)) {
	mu_test = sMin - i*dmus;
	lh_test = m_measurement.calcProb(n,mu_test);
	if(lh_test > lh_max) {
	  imax = i;
	  lh_max = lh_test;
	  mu_best = mu_test;
	}
	i++;
      }
      lowEdge = (i==nnew); //redfine lowEdge
    }
    //
    if (upEdge) { // ditto, but upper edge
      int nnew = ntst/10;
      i=0;
      while ((i<nnew)) {
	mu_test = sMin - (i+ntst)*dmus;
	lh_test = m_measurement.calcProb(n,mu_test);
	if(lh_test > lh_max) {
	  imax = i;
	  lh_max = lh_test;
	  mu_best = mu_test;
	}
	i++;
      }
      upEdge = (i==nnew-1); //redfine lowEdge
    }

    if (m_verbose>1) {
      if (upEdge || lowEdge) {
	std::cout << "WARNING: In FindBestMu -> might need to increase scan range in s! imax = " << imax << " ( " << ntst-1 << " )" << std::endl;
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
  if (m_verbose>2) {
    std::cout << "First 10 from best fit (mean,prob):" << std::endl;
    std::cout << m_bestMu.size() << ":" << m_bestMuProb.size() << std::endl;
    for (int i=0; i<10; i++) {
      std::cout << m_bestMu[i] << "\t" << m_bestMuProb[i] << std::endl;
    }
  }
  m_validBestMu = true;
}

void Pole::calcLh(double s) { // TODO: Need this one????
  //  double norm_p=0.0;
  for (int n=0; n<m_nBelt; n++) {
    m_muProb[n] = m_measurement.calcProb(n, s);
    //    norm_p += m_muProb[n]; // needs to be renormalised - NO!! 
  }
}

double Pole::calcLhRatio(double s, int & nbMin, int & nbMax) { //, double minMuProb) {
  double norm_p = 0;
  double lhSbest;
  bool upNfound  = false;
  bool lowNfound = false;
  int  nInBelt   = 0;
  int n=0;
  //
  nbMin = 0;
  nbMax = m_nBelt-1;
  //
  //
  while ((n<m_nBelt) && (!upNfound)) {
    m_muProb[n] =  m_measurement.calcProb(n, s);
    lhSbest = getLsbest(n);
    m_lhRatio[n]  = m_muProb[n]/lhSbest;
    norm_p += m_muProb[n]; // needs to be renormalised
    //
    if ((!lowNfound) && (m_muProb[n]>m_minMuProb)) {
      lowNfound=true;
      nbMin = n;
    } else {
      if ((nInBelt>1) && lowNfound && (m_muProb[n]<m_minMuProb)) {
	upNfound = true;
	nbMax = n-1;
      }
    }
    if (lowNfound && (!upNfound)) nInBelt++;
    n++;
    //    std::cout << "calcLhRatio: " << n << ", p = " << m_muProb[n] << ", lhSbest = " << lhSbest << std::endl;
  }
  //
  if (nbMin<m_nBeltMinUsed) m_nBeltMinUsed = nbMin;
  if (nbMax>m_nBeltMaxUsed) m_nBeltMaxUsed = nbMax;
  //
  return norm_p;
}


double Pole::calcLimit(double s) {
  int k,i;
  int nBeltMinUsed;
  int nBeltMaxUsed;
  //
  // Get RL(n,s)
  //
  double norm_p = calcLhRatio(s,nBeltMinUsed,nBeltMaxUsed);
  if (m_verbose>2) {
    std::cout << "Normalisation over n for s = " << s << " is " << norm_p << std::endl;
    if ((norm_p>1.5) || (norm_p<0.5)) {
      std::cout << "Normalisation off (" << norm_p << ") for s= " << s << std::endl;
      for (int n=0; n<nBeltMaxUsed; n++) {
	std::cout << "muProb[" << n << "] = " << m_muProb[n] << std::endl;
      }
    }
  }
  //
  k = m_measurement.getNObserved();

  if ((k>nBeltMaxUsed) || (k<nBeltMinUsed)) m_lhRatio[k] = 0.0;

  if (k>=m_nBelt) {
    k=m_nBelt; // WARNING::
    std::cout << "WARNING:: n_observed is larger than the maximum n used for R(n,s)!!" << std::endl;
    std::cout << "          -> increase nbelt such that it is more than n_obs = " << m_measurement.getNObserved() << std::endl;
  }									\
  if (m_verbose>2) std::cout << "Got nBelt range: " << nBeltMinUsed << ":" << nBeltMaxUsed << "( max = " << m_nBelt-1 << " )" << std::endl;
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
  m_sumProb = 0;
  if (m_verbose>9) std::cout << "\nSearch s: = " <<  s << std::endl;
  while (!done) {
    if(i != k) { 
      if(m_lhRatio[i] > m_lhRatio[k])  {
	m_sumProb  +=  m_muProb[i];
      }
      if (m_verbose>9) {
	std::cout << "RL[" << i << "] = " << m_lhRatio[i] << ", RLmax[" << k << "] = " << m_lhRatio[k] << ", sumP = " << m_sumProb << std::endl;
      }
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

double Pole::calcLimitOLD(double s) {
  int k,i;
  //
  double norm_p = 0;
  m_sumProb = 0;
  bool lowNfound=false;
  bool upNfound=false;
  int nInBelt=0;

  int nBeltMaxUsed = m_nBelt;
  int nBeltMinUsed = 0;
  //  const double minMuProb = 1e-5;
  //
  //  std::cout << "calcLimit for " << s << std::endl;
  if (usesMBT()) { // use method by Gary Hill
    double g,pbf;
    int n=0;
    while ((n<m_nBelt) && (!upNfound)) {
      //    for (int n=0; n<m_nBelt; n++) {
      if (!upNfound) {
	m_muProb[n] =  m_measurement.calcProb(n, s);
	if ((!lowNfound) && (m_muProb[n]>m_minMuProb)) {
	  lowNfound=true;
	  nBeltMinUsed = n;
	} else {
	  if ((nInBelt>1) && lowNfound && (m_muProb[n]<m_minMuProb)) {
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
	if (n>m_measurement.getBkgObs()) {
	  g = static_cast<double>(n);
	} else {
	  g = m_measurement.getBkgObs();
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
	m_muProb[n] =  m_measurement.calcProb(n, s);
	if ((!lowNfound) && (m_muProb[n]>m_minMuProb)) {
	  lowNfound=true;
	  nBeltMinUsed = n;
	} else {
	  if ((nInBelt>1) && lowNfound && (m_muProb[n]<m_minMuProb)) {
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

void Pole::calcConstruct(double s, bool verb) {
  int i;
  int nb1,nb2;
  //
  double norm_p = calcLhRatio(s,nb1,nb2);
  if (verb) {
    for (i=nb1; i<nb2; i++) {
      std::cout << "CONSTRUCT: " << s << "\t" << i << "\t" << m_lhRatio[i] << "\t" << m_muProb[i] << "\t" << norm_p << std::endl;
    }
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

double Pole::calcBelt(double s, int & n1, int & n2, bool verb) { //, double muMinProb) {
  //
  int nBeltMinUsed;
  int nBeltMaxUsed;
  double p;
  std::vector<int> index;
  //
  // Get RL(n,s)
  //
  //  double norm_p = 
  calcLhRatio(s,nBeltMinUsed,nBeltMaxUsed); //,muMinProb);
  //
  // Sort RL
  //
  sort_index(m_lhRatio,index,true); // reverse sort
  //
  // Calculate the probability for all n and the given s.
  // The Feldman-Cousins method dictates that for each n a
  // likelihood ratio (R) is calculated. The n's are ranked according
  // to this ratio. Values of n are included starting with that giving
  // the highest R and continuing with decreasing R until the total probability
  // matches the searched CL.
  // Below, the loop sums the probabilities for a given s and for all n with R>R0.
  // R0 is the likelihood ratio for n_observed.
  bool done=false;
  int nmin=-1;
  int nmax=-1;
  int n;//imax = nBeltMaxUsed - nBeltMinUsed;
  double sumProb = 0;
  int i=0;
  //  int nn=nBeltMaxUsed-nBeltMinUsed+1;
  //
  while (!done) {
    n = index[i];
    if ((n>=nBeltMinUsed) && (n<=nBeltMaxUsed)) {
      p = m_muProb[n];
      sumProb +=p;
      if ((n<nmin)||(nmin<0)) nmin=n;
      if ((n>nmax)||(nmax<0)) nmax=n;
    }
    //    std::cout << "calcBelt: " << i << " , " << n << " , " << m_lhRatio[n] << std::endl;
    //
    i++;
    done = ((i==nBeltMaxUsed) || sumProb>m_cl);
  }
  if ((nmin<0) || ((nmin==0)&&(nmax==0))) {
    nmin=0;
    nmax=1;
    sumProb=1.0;
  }
  n1 = nmin;
  n2 = nmax;
  if (verb) {
    std::cout << "CONFBELT: " << s << "\t" << n1 << "\t" << n2 << "\t" << sumProb << std::endl;
    //	      << m_lhRatio[n1] << "\t" << m_lhRatio[n2] << "\t" << norm_p << "\t"
    //	      << index[0] << "\t" << m_lhRatio[index[0]] << " , nBelt: " << nBeltMinUsed << " , " << nBeltMaxUsed << std::endl;
  }
  return sumProb;
}

//*********************************************************************//
//*********************************************************************//
//*********************************************************************//

void Pole::findPower() {
  double muTest;
  std::vector< std::vector<double> > fullConstruct;
  std::vector< double > probVec;
  std::vector< int > n1vec;
  std::vector< int > n2vec;
  int n1,n2;
  int nhyp = m_hypTest.n();
  double sumP;
  if (m_verbose>-1) std::cout << "Make full construct" << std::endl;
  for (int i=0; i<nhyp; i++) {
    muTest = m_hypTest.min() + i*m_hypTest.step();
    //    if (m_verbose>-1) std::cout << "Hypothesis: " << muTest << std::endl;
    sumP = calcBelt(muTest,n1,n2,false); //,-1.0);
    fullConstruct.push_back(m_muProb);

    n1vec.push_back(n1);
    n2vec.push_back(n2);
  }
  if (m_verbose>-1) std::cout << "Calculate power(s)" << std::endl;
  // calculate power P(n not accepted by H0 | H1)
  int n01,n02,n;
  //,n11,n12,nA,nB;
  double power,powerm,powerp;
  bool usepp;
  int np,nm;
  //  int npTot;
  int i1 = int(m_sTrue/m_hypTest.step());
  int i2 = i1+1;
  for (int i=i1; i<i2; i++) { // loop over all H0
    //    std::cout << fullConstruct[i].size() << std::endl;
    //    if (m_verbose>-1) std::cout << "H0 index = " << i << " of " << nhyp << std::endl;
    n01 = n1vec[i]; // belt at H0
    n02 = n2vec[i];
    powerm=0.0;
    powerp=0.0;
    nm = 0;
    np = 0;
    for (int j=0; j<nhyp; j++) { // loop over all H1
//       n11 = n1vec[j]; // belt at H1
//       n12 = n2vec[j];
      usepp = (j>i);
      if (usepp) {
	np++;
      } else {
	nm++;
      }
//       if (usepp) { // H1>H0
// 	nA = 0;
// 	nB = n01; // lower lim H0
// 	np++;
//       } else {
// 	nA = n12;     // upper limit for H1
// 	nB = m_nBelt; // man N
// 	nm++;
//       }
      //      if (nA!=nB) std::cout << "N loop = " << nA << "\t" << nB << std::endl;
      powerm=0.0;
      powerp=0.0;
      double pcl=0.0;
      probVec = fullConstruct[j];
      for ( n=0; n<n01; n++) { // loop over all n outside acceptance of H0
	if (usepp) {
	  powerp += probVec[n]; // P(n|H1) H1>H0
	} else {
	  powerm += probVec[n]; // P(n|H1) H1<H0
	}
      }
      for (n=n01; n<=n02; n++) {
	pcl += probVec[n];
      }
      for ( n=n02+1; n<m_nBelt; n++) {
	if (usepp) {
	  powerp += probVec[n]; // P(n|H1) H1>H0
	} else {
	  powerm += probVec[n]; // P(n|H1) H1<H0
	}
      }
      power = powerp+powerm;
      std::cout << "POWER:\t" << m_hypTest.getVal(i) << "\t" << m_hypTest.getVal(j) << "\t" << power << "\t" << "0\t" << pcl << std::endl;
    }
    power = powerm+powerp;
    if (nm+np==0) {
      power = 0.0;
    } else {
      power = power/double(nm+np);
    }
    if (nm==0) {
      powerm = 0.0;
    } else {
      powerm = powerm/double(nm);
    }
    if (np==0) {
      powerp = 0.0;
    } else {
      powerp = powerp/double(np);
    }
  }
}

void Pole::findConstruct() {
  double mu_test;
  int i = 0;
  bool done = (i==m_hypTest.n());
  //
  while (!done) {
    mu_test = m_hypTest.min() + i*m_hypTest.step();
    calcConstruct(mu_test,true);
    i++;
    done = (i==m_hypTest.n()); // must loop over all hypothesis
  }
}

int Pole::findNMin() { // calculates the minimum N rejecting s = 0.0
  int n1,n2;
  double sumP = calcBelt(0.0,n1,n2,true);//,-1.0);
  //
  std::cout << "NMIN0:\t" << n2 << "\t" << sumP << std::endl;
  return n2;
}

void Pole::findBelt() {
  double mu_test;
  int i = 0;
  bool done = (i==m_hypTest.n());
  //
  int n1,n2;
  while (!done) {
    mu_test = m_hypTest.min() + i*m_hypTest.step();
    calcBelt(mu_test,n1,n2,true);
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
  findNMin();
  if (m_measurement.getNObserved()>=m_nBelt) return false;
  double mu_test;
  int i = 0;
  bool done = (i==m_hypTest.n());
  //  bool nuppOK = false;
  //  int nupp=0;
  //
  //
  double p;
  //  int n1,n2;
  bool firstLower=true;
  int nAfterFirstLow=0;
  //
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
//    if (m_foundUpper) { // sample a few afterwards just to make sure we get the last one
//       nupp++;
//       nuppOK = ((nupp>=m_nUppLim) || (m_nUppLim<1));
//      nuppOK = true; // SKIP THE SCAN...
//    }
    done = ((i==m_hypTest.n()) || // Done if last hypothesis reached
	    //	    (nuppOK)           || //
	    (m_foundUpper) ||
	    (!normOK(p)));        // normalisation not OK
  }
  if (m_verbose>1) {
    if (limitsOK()) {
      // && (!heavyDBG)) {
      std::cout << "LIMITS(N,e,b,l,u): ";
      coutFixed(4,m_measurement.getNObserved());  std::cout << "\t";
      coutFixed(4,m_measurement.getEffObs());    std::cout << "\t";
      coutFixed(4,m_measurement.getBkgObs());    std::cout << "\t";
      coutFixed(4,m_lowerLimit); std::cout << "\t";
      coutFixed(4,m_upperLimit); std::cout << std::endl;
      //      heavyDBG = false;
      //      findLimits();
      //      heavyDBG = false;
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
    } else {
      std::cout << "Limits NOT OK!" << std::endl;
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
  //  bool nuppOK = false;
  //  int  nupp=0;
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
      coutFixed(4,m_measurement.getEffObs());    std::cout << "\t";
      coutFixed(4,m_measurement.getBkgObs());    std::cout << "\t";
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
  m_measurement.initIntNuisance();
  initBeltArrays();
  if (m_verbose>0) std::cout << "Constructing integral" << std::endl;
  m_measurement.fillIntNuisance();  // fill f(x)dx per nuisance parameter
  m_measurement.initNuisanceWeights(); // calculate all f(x)dx*g(y)dy... for all combinations of (x,y,...)
}

bool Pole::analyseExperiment() {
  bool rval=false;
  //  initAnalysis();
  if (usesFHC2()) {
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
    std::cout << " Max N(belt) set  : " << m_nBelt << std::endl;
    std::cout << " Max N(belt) used : " << m_nBeltMaxUsed << std::endl;
    std::cout << "-------------------------------------------------" << std::endl;
    std::cout << " Nobs  \t  Eff   \t Bkg" << std::endl;
    std::cout << "-------------------------------------------------" << std::endl;
  }
  coutFixed(4,m_measurement.getNObserved()); std::cout << "\t";
  coutFixed(6,m_measurement.getEffObs()); std::cout << "\t";
  coutFixed(6,m_measurement.getEffPdfSigma()); std::cout << "\t";
  coutFixed(6,m_measurement.getBkgObs()); std::cout << "\t";
  coutFixed(6,m_measurement.getBkgPdfSigma()); std::cout << "\t";
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
  std::cout << " Efficiency meas    : " << m_measurement.getEffObs() << std::endl;
  std::cout << " Efficiency sigma   : " << m_measurement.getEffPdfSigma() << std::endl;
  std::cout << " Efficiency dist    : " << PDF::distTypeStr(m_measurement.getEffPdfDist()) << std::endl;
  std::cout << "----------------------------------------------\n";
  std::cout << " Background meas    : " << m_measurement.getBkgObs() << std::endl;
  std::cout << " Background sigma   : " << m_measurement.getBkgPdfSigma() << std::endl;
  std::cout << " Background dist    : " << PDF::distTypeStr(m_measurement.getBkgPdfDist()) << std::endl;
  std::cout << "----------------------------------------------\n";
  std::cout << " Bkg-Eff correlation: " << m_measurement.getBEcorr() << std::endl;
  std::cout << "----------------------------------------------\n";
  std::cout << " Int. eff. min      : " << getEffIntMin() << std::endl;
  std::cout << " Int. eff. max      : " << getEffIntMax() << std::endl;
  std::cout << " Int. eff. N pts    : " << getEffIntN() << std::endl;
  std::cout << "----------------------------------------------\n";
  std::cout << " Int. bkg. min      : " << getBkgIntMin() << std::endl;
  std::cout << " Int. bkg. max      : " << getBkgIntMax() << std::endl;
  std::cout << " Int. bkg. N pts    : " << getBkgIntN() << std::endl;
  std::cout << "----------------------------------------------\n";
  std::cout << " Test hyp. min      : " << m_hypTest.min() << std::endl;
  std::cout << " Test hyp. max      : " << m_hypTest.max() << std::endl;
  std::cout << " Test hyp. step     : " << m_hypTest.step() << std::endl;
  std::cout << "----------------------------------------------\n";
  std::cout << " Min. probability   : " << m_minMuProb << std::endl;
  std::cout << "----------------------------------------------\n";
  if (m_suggestBelt) {
    std::cout << " Belt N max         : variable" << std::endl;
  } else {
    std::cout << " Belt N max         : " << m_nBelt << std::endl;
  }
  std::cout << " Step mu_best       : " << m_dmus << std::endl;
  std::cout << " Max N, mu_best     : " << m_nmusMax << std::endl;
  std::cout << "----------------------------------------------\n";
  std::cout << " Method             : ";
  switch (m_method) {
  case RL_FHC2:
    std::cout << "FHC2";
    break;
  case RL_MBT:
    std::cout << "MBT";
    break;
  default:
    std::cout << "Unknown!? = " << m_method;
    break;
  }
  std::cout << std::endl;
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
  std::cout << "4. Symptom: probability of lower/upper limit diverges from unity." << std::endl;
  std::cout << "   -> if Poisson table is used, insufficient precision." << std::endl;
  std::cout << "   -> minimum probability too large; min = " << m_minMuProb << std::endl;
  std::cout << "Input:" << std::endl;
  std::cout << "   N(obs)     = " << getNObserved() << std::endl;
  std::cout << "Results:" << std::endl;
  std::cout << "   probability (should be = CL)  = " << getSumProb() << std::endl;
  std::cout << "   lower lim norm (should be 1)  = " << getLowerLimitNorm() << std::endl;
  std::cout << "   upper lim norm (ditto)        = " << getUpperLimitNorm() << std::endl;
  std::cout << "   lower lim                     = " << getLowerLimit() << std::endl;
  std::cout << "   upper lim                     = " << getUpperLimit() << std::endl;
}
