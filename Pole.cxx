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

Pole::Pole() {
  m_useLogNormal = false;
  m_stepMin = 0.001;
  m_s2pi = sqrt(2.0L*M_PIl);
  //
  m_dmus = 0.01;
  //
  m_nObserved = -1;
  m_nBeltMax = 50;
  //
  m_effScaleInt = 5.0;
  m_bkgScaleInt = 5.0;
  //
  m_validInt  = false;
  m_nInt      = 0;
  m_weightInt = 0;
  m_effInt    = 0;
  m_bkgInt    = 0;
  //
  m_validBestMu = false;
  m_nBeltMax = 0;
  m_muProb = 0;
  m_bestMuProb = 0;
  m_bestMu = 0;
  m_lhRatio = 0;
  //
  // init list of suggested nBelt
  //
  m_nBeltList.push_back(18); //0
  m_nBeltList.push_back(20);
  m_nBeltList.push_back(20);
  m_nBeltList.push_back(22);
  m_nBeltList.push_back(22);
  m_nBeltList.push_back(23); //5
  m_nBeltList.push_back(25);
  m_nBeltList.push_back(32);
  m_nBeltList.push_back(38);
  m_nBeltList.push_back(38);
  m_nBeltList.push_back(40); //10
}
Pole::~Pole() {
  // delete arrays
  if (m_weightInt)  delete [] m_weightInt;
  if (m_effInt)     delete [] m_effInt;
  if (m_bkgInt)     delete [] m_bkgInt;
  //
  if (m_bestMuProb) delete [] m_bestMuProb;
  if (m_muProb)     delete [] m_muProb;
  if (m_bestMu)     delete [] m_bestMu;
  if (m_lhRatio)    delete [] m_lhRatio;
}

void Pole::setEffDist(double mean,double sigma, bool nodist) {
  m_effMeas   = mean;
  m_effSigma  = sigma;
  m_effNoDist = nodist;
  m_validInt    = false;
  m_validBestMu = false;
}

void Pole::setBkgDist(double mean,double sigma, bool nodist) {
  m_bkgMeas   = mean;
  m_bkgSigma  = sigma;
  m_bkgNoDist = nodist;
  m_validInt    = false;
  m_validBestMu = false;
}

void Pole::setEffInt(double low, double high, double step) {
  if (step<0) {
    setEffInt();
  } else {
    m_validInt = false;
    m_validBestMu = false;
    //
    // Will set the range for the efficiency integration.
    // If step<0, then set range according to the selected efficiency distribution.
    // If fixed, force it to the fixed mean value.
    //
    if (m_effNoDist) {
      //      std::cout << "Efficiency known exactly, excluding it from integral" << std::endl;
      low  = m_effMeas;
      high = m_effMeas;
      step = 1.0;
    } else {
      if (high>low) {
	if (step<m_stepMin) step=m_stepMin;
      } else {
	high = low;
	step = 1.0;
      }
    }
    if (m_effSigma>0) m_effScaleInt = 0.5*(high-low)/m_effSigma;
    m_effRangeInt.setRange(low,high,step,m_stepMin);
    m_effIntMin = m_effRangeInt.min();
    m_effIntMax = m_effRangeInt.max();
  }
}

void Pole::setEffInt(double scale, double step) {
  double low,high;
  m_validInt = false;
  m_validBestMu = false;
  //
  // Will set the range for the efficiency integration.
  // If step<0, then set range according to the selected efficiency distribution.
  // If fixed, force it to the fixed mean value.
  //
  if (m_effNoDist) {
    //    std::cout << "Efficiency known exactly, excluding it from integral" << std::endl;
    low  = m_effMeas;
    high = m_effMeas;
    step = 1.0;
  } else {
    if (scale<1.0) {
      scale = m_effScaleInt;
    } else {
      m_effScaleInt = scale;
    }
    low  = m_effMeas - scale*m_effSigma;
    high = m_effMeas + scale*m_effSigma;
    if (low<0) low = 0;
    if (step<=0.0) {
      step = (high-low)/20.0;     // minimum 21 pts
      if (step>0.1) step = 0.1; // not too coarse
    }
    if (step<m_stepMin) step = m_stepMin;

  }
  m_effRangeInt.setRange(low,high,step,m_stepMin);
  m_effIntMin = m_effRangeInt.min();
  m_effIntMax = m_effRangeInt.max();
}

void Pole::setBkgInt(double low, double high, double step) {
  if (step<0) {
    setBkgInt();
  } else {
    m_validInt = false;
    m_validBestMu = false;

    if (m_bkgNoDist) {
      low  = m_bkgMeas;
      high = m_bkgMeas;
      step = 1.0;
    } else {
      if (high>low) {
	if (step<m_stepMin) step=m_stepMin;
      } else {
	low = high;
	step = 1.0;
      }
    }
    if (m_bkgSigma>0) m_bkgScaleInt = 0.5*(high-low)/m_bkgSigma;
    m_bkgRangeInt.setRange(low,high,step,m_stepMin);
    m_bkgIntMin = m_bkgRangeInt.min();
    m_bkgIntMax = m_bkgRangeInt.max();
  }
}

void Pole::setBkgInt(double scale, double step) {
  double low,high;
  m_validInt = false;
  m_validBestMu = false;

  if (m_bkgNoDist) {
    low  = m_bkgMeas;
    high = m_bkgMeas;
    step = 1.0;
  } else {
    if (scale<1.0) {
      scale = m_bkgScaleInt;
    } else {
      m_bkgScaleInt = scale;
    }
    low  = m_bkgMeas - scale*m_bkgSigma;
    high = m_bkgMeas + scale*m_bkgSigma;
    if (low<0) low = 0;
    if (step<=0.0) {
      step = (high-low)/20.0;     // minimum 21 points
      if (step>0.1) step = 0.1; // not too coarse
    }
    if (step<m_stepMin) step = m_stepMin;
  }
  m_bkgRangeInt.setRange(low,high,step,m_stepMin);
  m_bkgIntMin = m_bkgRangeInt.min();
  m_bkgIntMax = m_bkgRangeInt.max();
}

void Pole::setTestHyp(double low, double high, double step) {
  int n;
  if (high>low) {
    if (step<m_stepMin) step=m_stepMin;
  } else {
    low = high;
    step = 1.0;
    n = 1;
  }
  m_hypTest.setRange(low,high,step,m_stepMin);
}
//
bool Pole::checkParams() {
  bool rval=true;
  // check efficiency distribution
  // check background distribution
  // check true signal range
  // check mu_test range
  // check efficiency integration
  std::cout << "Checking efficiency integration - ";
  double dsLow, dsHigh;
  dsLow  = (m_effMeas - m_effRangeInt.min())/m_effSigma;
  dsHigh = (m_effRangeInt.max()  - m_effMeas)/m_effSigma;
  if ( (!m_effNoDist) && ( ((dsLow<4.0)&&(m_effRangeInt.min()>0)) ||
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
  dsLow  = (m_bkgMeas - m_bkgRangeInt.min())/m_bkgSigma;
  dsHigh = (m_bkgRangeInt.max()  - m_bkgMeas)/m_bkgSigma;
  if ( (!m_bkgNoDist) && ( ((dsLow<4.0)&&(m_bkgRangeInt.min()>0)) ||
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
  m_poisson.init(nlambda, nn, lmbmax);
}

void Pole::initGauss(int ndata, double mumax) {
  m_gauss.init(ndata, mumax);
}

//
// Must be called after setEffInt, setBkfInt and setTestHyp
//
void Pole::initIntArrays() {
  m_nInt = m_effRangeInt.n()*m_bkgRangeInt.n();
  if ((m_weightInt==0) || (m_nInt>m_nIntMax)) {
    m_nIntMax = m_nInt*2; // make a safe margin
    if (m_weightInt) { // arrays needs resizing
      delete [] m_weightInt;
      delete [] m_effInt;
      delete [] m_bkgInt;
    }
    m_weightInt = new double[m_nIntMax];
    m_effInt    = new double[m_nIntMax];
    m_bkgInt    = new double[m_nIntMax];
  }
}

int Pole::suggestBelt() {
  int rval=50;
  int nbelt = static_cast<int>(m_nBeltList.size());
  if ((m_nObserved>=0) &&(m_nObserved<nbelt)) {
    rval = m_nBeltList[m_nObserved];
  }
  return rval;
}

void Pole::initBeltArrays(bool suggest) {
  if (suggest || (m_nBelt<1)) m_nBelt = suggestBelt();
  if (m_verbose>1) std::cout << "Belt range: " << m_nBelt << std::endl;
  if ((m_muProb==0) || (m_nBelt>m_nBeltMax) || (m_nBelt<1)) {
    m_nBeltMax = m_nBelt*2;
    if (m_muProb) { // arrays needs resizing
      delete [] m_muProb;
      delete [] m_bestMuProb;
      delete [] m_bestMu;
      delete [] m_lhRatio;
    }
    m_muProb     = new double[m_nBeltMax];
    m_bestMuProb = new double[m_nBeltMax];
    m_bestMu     = new double[m_nBeltMax];
    m_lhRatio    = new double[m_nBeltMax];
  }
}

void Pole::initIntegral() {
  if (m_validInt) return;
  //
  int i,j;
  double effs;
  double bg;
  double eff_prob, bkg_prob, norm_prob;
  int count=0;
  norm_prob = 1.0;
  //////////////////////////////////
  double effMean,effSigma;   //TMP
  double bkgMean,bkgSigma;   //TMP
  if (m_useLogNormal) {  //TMP
    effMean  = log(m_effMeas*m_effMeas/sqrt(m_effSigma*m_effSigma + m_effMeas*m_effMeas));
    effSigma = sqrt(log(m_effSigma*m_effSigma/m_effMeas/m_effMeas+1));
    bkgMean  = log(m_bkgMeas*m_bkgMeas/sqrt(m_bkgSigma*m_bkgSigma + m_bkgMeas*m_bkgMeas));
    bkgSigma = sqrt(log(m_bkgSigma*m_bkgSigma/m_bkgMeas/m_bkgMeas+1));
  } else {
    effMean  = m_effMeas;
    effSigma = m_effSigma;
    bkgMean  = m_bkgMeas;
    bkgSigma = m_bkgSigma;
  }
  //////////////////////////////////
  //
  if (m_verbose>1) {
    std::cout << "eff dist (mu,s)= " << m_effMeas << "\t"
	      << m_effSigma << std::endl;
    std::cout << "bkg dist (mu,s)= " << m_bkgMeas << "\t"
	      << m_bkgSigma << std::endl;
  }
  double dedb=1.0;
  if ((m_effRangeInt.n() != 1) && (m_bkgRangeInt.n() != 1)) { // both eff and bkg vary
    dedb = m_bkgRangeInt.step()*m_effRangeInt.step();
  } else {
    if ((m_effRangeInt.n() != 1) && (m_bkgRangeInt.n() == 1)) { // bkg const
      dedb = m_effRangeInt.step();
    } else {
      if ((m_effRangeInt.n() == 1) && (m_bkgRangeInt.n() != 1)) {  // eff const
	dedb = m_bkgRangeInt.step();
      }
    }
  }
  for (i=0;i<m_effRangeInt.n();i++) { // Loop over all efficiency points
    effs =  i*m_effRangeInt.step() + m_effRangeInt.min();
    if ( m_effRangeInt.n() == 1 ) {
      eff_prob = 1.0;
    } else {
      if (m_useLogNormal) {
	eff_prob = logNormal(effs,effMean,effSigma); // TMP
      } else {
	eff_prob = m_gauss.getVal(effs,m_effMeas,m_effSigma);
      }
    }
    for (j=0;j<m_bkgRangeInt.n();j++) { // Loop over all background points
      bg =  j*m_bkgRangeInt.step() + m_bkgRangeInt.min();
      if ( m_bkgRangeInt.n() == 1 ) {
	bkg_prob = 1.0;
      } else {
	if (m_useLogNormal) {
	  bkg_prob = logNormal(bg,bkgMean,bkgSigma);
	} else {
	  bkg_prob = m_gauss.getVal(bg,m_bkgMeas,m_bkgSigma);
	}
      }
      norm_prob = eff_prob*bkg_prob*dedb;
      //
      // Fill the arrays (integral)
      //
      if( (m_effRangeInt.n() == 1) && (m_bkgRangeInt.n() == 1) ) {
	m_bkgInt[count]    = m_bkgMeas;
	m_effInt[count]    = m_effMeas;
	m_weightInt[count] = 1.0;
      } else {
	m_bkgInt[count]    = bg;
	m_effInt[count]    = effs;
	m_weightInt[count] = norm_prob;
      }
      count++;
      if (count>m_nInt) {
	std::cout << "WARNING: Array index overflow. count = " << count
		  << ", n_e = " << m_effRangeInt.n()
		  << ", n_b = " << m_bkgRangeInt.n() << std::endl;
      }
    }
  }
  if (m_verbose>1) {
    std::cout << "First 10 from double integral (bkg,eff,w):" << std::endl;
    std::cout << "Eff: " << m_effRangeInt.n()   << "\t" << m_effRangeInt.step()
	      << "\t"    << m_effRangeInt.min() << "\t" << m_effRangeInt.max() << std::endl;
    std::cout << "Bkg: " << m_bkgRangeInt.n()   << "\t" << m_bkgRangeInt.step()
	      << "\t"    << m_bkgRangeInt.min() << "\t" << m_bkgRangeInt.max() << std::endl;
    for (i=0; i<10; i++) {
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
  if(n<m_bkgMeas) {
    m_bestMu[n] = 0; // best mu is 0
    m_bestMuProb[n] = calcProb(n,0);
  } else {
    mu_s_max = double(n)-m_bkgMeas;
    //    mu_s_max = (double(n) - m_bkgRangeInt.min())/m_effRangeInt.min();
    //    mu_s_min = (double(n) - m_bkgRangeInt.max())/m_effRangeInt.max();
    mu_s_min = (double(n) - m_bkgIntMax)/m_effIntMax;
    if(mu_s_min<0) {mu_s_min = 0.0;}
    //    dmu_s = 0.01; // HARDCODED:: Change!
    int ntst = 1+int((mu_s_max-mu_s_min)/m_dmus);
    if (m_verbose>1) std::cout << "FindBestMu range: " << ntst << " [" << mu_s_min << "," << mu_s_max << "] => ";
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
    for (int i=0; i<10; i++) {
      std::cout << m_bestMu[i] << "\t" << m_bestMuProb[i] << std::endl;
    }
  }
  m_validBestMu = true;
}

void Pole::calcLimit(double s) {
  int k,i;
  //
  double norm_p = 0;
  m_sumProb = 0;
  //
  //  std::cout << "calcLimit for " << s << std::endl;
  for (int n=0; n<m_nBelt; n++) {
    m_muProb[n] =  calcProb(n, s);
    m_lhRatio[n]  = m_muProb[n]/m_bestMuProb[n]; // POLE
    norm_p += m_muProb[n]; // needs to be renormalised
  }
  //
  for(i=0;i<m_nBelt;i++) {
    m_muProb[i] = m_muProb[i]/norm_p;
    //    for(k=0;k<m_nBelt;k++) {
    k = m_nObserved;
    if(i != k) { 
      if(m_lhRatio[i]  > m_lhRatio[k])  {
	m_sumProb  +=  m_muProb[i];
      }                
      //    }
    }
  }
  if (m_sumProb<m_cl) {
    if (m_foundLower) {
      m_upperLimit = s;
      m_foundUpper = true;
    } else {
      m_lowerLimit = s;
      m_foundLower = true;
      m_foundUpper = false;
    }
  }
//  else {
//     if (m_foundUpper) {
//       if (!m_foundUpper) {
// 	//	std::cout << "Last upper found" << std::endl;
// 	m_foundUpper = true;
//       } else {
// 	//	std::cout << "WARNING:: found a new set of upper limits!!!" << std::endl;
//       }
//     }
//   }
}

bool Pole::findLimits() {
  double mu_test;
  int i = 0;
  bool done = (i==m_hypTest.n());
  //
  m_foundLower = false;
  m_foundUpper = false;
  m_lowerLimit = 0;
  m_upperLimit = 0;
  //
  while (!done) {
    mu_test = m_hypTest.min() + i*m_hypTest.step();
    calcLimit(mu_test);
    i++;
    done = (i==m_hypTest.n());
  }
  if (m_verbose>1) {
    if (m_foundUpper) {
      std::cout << std::fixed << std::setprecision(4)
		<< "LIMITS(N,e,b,l,u): " << m_nObserved << "\t"
		<< m_effMeas << "\t" << m_bkgMeas << "\t"
		<< m_lowerLimit << "\t" << m_upperLimit << std::endl;
    } else {
      std::cout << "Limits NOT OK!" << std::endl;
    }
  }
  return (m_foundUpper);
}

bool Pole::findCoverageLimits() {
  //
  // Same as findLimits() but it will stop the scan as soon it is clear
  // that the true s (m_sTrueMean) is inside or outside the range.
  //
  double mu_test;
  int i = 0;
  bool done = (i==m_hypTest.n());
  bool decided = false;
  //
  m_foundLower = false;
  m_foundUpper = false;
  m_lowerLimit = 0;
  m_upperLimit = 0;
  //
  while (!done) {
    mu_test = m_hypTest.min() + i*m_hypTest.step();
    calcLimit(mu_test);
    i++;
    if ((!m_foundLower) && (mu_test>m_sTrue)) { // for sure outside
      decided = true;
      m_lowerLimit = mu_test;     // create a fake set of limits
      m_upperLimit = mu_test+0.1;
      m_foundLower = true;
      m_foundUpper = true;
    } else if ((m_foundUpper) && (m_upperLimit>m_sTrue)) { // for sure inside
      decided = true;
    }
    //   done = ((limitOk) || (i==m_hypTest.n()));
    done = ((decided) || (i==m_hypTest.n())); // loop until limit is found or end of
    //    done = (i==m_hypTest.n());
  }
  if (m_verbose>1) {
    if (m_foundUpper) {
      std::cout << std::fixed << std::setprecision(4)
		<< "COVLIMITS(N,s,e,b,l,u): " << m_nObserved << "\t" << m_sTrue << "\t"
		<< m_effMeas << "\t" << m_bkgMeas << "\t"
		<< m_lowerLimit << "\t" << m_upperLimit << std::endl;
    } else {
      std::cout << "Limits NOT OK!" << std::endl;
    }
  }
  return decided;
}

void Pole::analyseExperiment() {
  if (m_verbose>0) std::cout << "Initialise arrays" << std::endl;
  initIntArrays();
  initBeltArrays();
  if (m_verbose>0) std::cout << "Constructing integral" << std::endl;
  initIntegral();  // generates predefined parts of double integral ( eq.7)
  if (m_verbose>0) std::cout << "Finding s_best" << std::endl;
  findAllBestMu(); // loops
  if (m_verbose>0) std::cout << "Calculating limit" << std::endl;
  if (m_coverage) {
    findCoverageLimits();
  } else {
    findLimits();
  }
}

void Pole::printLimit(bool doTitle) {
  if (doTitle) {
    std::cout << " Nobs  \t  Eff   \t Bkg" << std::endl;
    std::cout << "-------------------------------------------------" << std::endl;
  }
  std::cout << std::fixed << std::setw(4)
	    << m_nObserved << "\t"
	    << std::fixed << std::setprecision(6)
	    << m_effMeas << "\t"
	    << m_effSigma << "\t"
	    << m_bkgMeas << "\t"
	    << m_bkgSigma << "\t[ "
	    << std::fixed << std::setprecision(2)
	    << m_lowerLimit << ", "
	    << m_upperLimit << " ]" << std::endl;
}

void Pole::printSetup() {
  std::cout << "\n";
  std::cout << "================ P O L E ==================\n";
  std::cout << " Confidence level   : " << m_cl << std::endl;
  std::cout << " N observed         : " << m_nObserved << std::endl;
  std::cout << "----------------------------------------------\n";
  std::cout << " Coverage friendly  : " << yesNo(m_coverage) << std::endl;
  std::cout << " True signal        : " << m_sTrue << std::endl;
  std::cout << "----------------------------------------------\n";
  std::cout << " Efficiency meas    : " << m_effMeas << std::endl;
  std::cout << " Efficiency sigma   : " << m_effSigma << std::endl;
  std::cout << " Efficiency no dist : " << yesNo(m_effNoDist) << std::endl;
  std::cout << "----------------------------------------------\n";
  std::cout << " Background meas    : " << m_bkgMeas << std::endl;
  std::cout << " Background sigma   : " << m_bkgSigma << std::endl;
  std::cout << " Background no dist : " << yesNo(m_bkgNoDist) << std::endl;
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
  std::cout << " Test hyp. N pts    : " << m_hypTest.n() << std::endl;
  std::cout << "----------------------------------------------\n";
  std::cout << " Belt N max         : " << m_nBelt << std::endl;
  std::cout << " Step mu_best       : " << m_dmus << std::endl;
  std::cout << "==============================================\n";
  //
  if (m_useLogNormal) {
    std::cout << "\nNOTE: USING LOG-NORMAL!!\n" << std::endl;
  }
}
