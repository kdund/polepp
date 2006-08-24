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
  m_scaleLimit=1.0;
  m_inputFile="";
  m_inputFileLines=0;
  m_printLimStyle=0; // computer readable
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
  // N = 1
  // eff = 1.0, Gaus(1.0,0.1)
  // bkg = 0.0
  //
  setNObserved(1);
  setEffPdf(1.0,0.1,PDF::DIST_GAUS);
  setEffObs();
  //  m_measurement.setEffPdf(1.0,0.1,PDF::DIST_GAUS);
  //  m_measurement.setEffObs(); // sets the mean from PDF as the observed efficiency
  setBkgPdf(0.0,0.0,PDF::DIST_NONE);
  setBkgObs();
  //  m_measurement.setBkgPdf(0.0,0.0,PDF::DIST_NONE);
  //  m_measurement.setBkgObs();
  setEffPdfBkgCorr(0.0);
  //
  m_minMuProb = 1e-6;
  //
  m_bestMuStep = 0.01;
  m_bestMuNmax = 100;
  //
  m_measurement.setEffInt(5.0,21);
  m_measurement.setBkgInt(5.0,21);
  //
  m_limitHypStep = 0.01;
  setTestHyp(0.01);
  //
  m_validBestMu = false;
  m_nBelt        = 50;
  m_nBeltUsed    = 5; // not really needed to be set here - set in initBeltArrays()
  m_nBeltMinUsed = m_nBeltUsed; // idem
  m_nBeltMaxUsed = 0;
}

Pole::~Pole() {
}

void Pole::execute() {
  if (m_inputFile.size()>0) {
    exeFromFile();
  } else {
    exeEvent(true);
  }
}

void Pole::exeEvent(bool first) {
  initAnalysis();
  //
  if (first) printSetup();
  if (analyseExperiment()) {
    printLimit(first);
  } else {
    printFailureMsg();
  }
}

void Pole::exeFromFile() {
  std::ifstream inpf( m_inputFile.c_str() );
  if ( ! inpf.is_open() ) {
    std::cout << "ERROR: File <" << m_inputFile << "> could not be found/opened" << std::endl;
    return;
  }
  std::cout << "--- Reading data from input file: " << m_inputFile << std::endl;
  int nch;
  std::string dataLine;
  int nskipped = 0;
  int nread    = 0;
  int n;
  double em,es,bm,bs;
  double y,z,t,u;
  int ed,bd;
  bool first=true;
  int nlines=0;
  TOOLS::Timer loopTime;
  loopTime.start();
  bool notDone=true;
  //
  while (notDone && ((nch=inpf.peek())>-1)) {
    getline(inpf,dataLine);
    if ((dataLine[0]!='#') && (dataLine.size()>0)) {
      std::istringstream sstr(dataLine);
      //  sstr >> n >> ed >> em >> es >> bd >> bm >> bs;
      sstr >> n >> y >> z >> t >> u;
      if (sstr) {
        if ((m_verbose>0) && (nlines<10)) {
          std::cout << "LINE: ";
          //          std::cout << n << "\t" << ed << "\t" << em << "\t" << es << "\t" << bd << "\t" << bm << "\t" << bs << std::endl;
          std::cout << n << "\t" << y << "\t" << z << "\t" << t << "\t" << u  << std::endl;
        }
        nread++;
        setNObserved(n*t);
        setEffPdf( z, z, PDF::DIST_POIS );
        setEffObs();
        setBkgPdf( y, y, PDF::DIST_POIS );
        setBkgObs();
        setScaleLimit(u/t);

//         setNObserved(n);
//         setEffPdf( em, es, PDF::DISTYPE(ed) );
//         setEffObs();
//         setBkgPdf( bm, bs, PDF::DISTYPE(bd) );
//         setBkgObs();
        //      
        exeEvent(first);
        first=false;
        nlines++;
        if ((m_inputFileLines>0) && (nlines>=m_inputFileLines)) notDone = false; // enough lines read
      } else {
        if (nskipped<10) {
          std::cout << "Skipped line " <<  nskipped+1 << " in input file" << std::endl;
        }
        nskipped++;
      }
    }
  }
  loopTime.stop();
  loopTime.printUsedTime();
  loopTime.printUsedClock(nlines);
  if (first) {
    std::cout << "Failed processing any lines in the given input file." << std::endl;
  }
}

bool Pole::checkEffBkgDists() {
  bool change=false;
  // If ONLY one is PDF::DIST_GAUS2D, make bkgDist == effDist
  if ( ((getEffPdfDist() == PDF::DIST_GAUS2D) && (getBkgPdfDist() != PDF::DIST_GAUS2D)) ||
       ((getEffPdfDist() != PDF::DIST_GAUS2D) && (getBkgPdfDist() == PDF::DIST_GAUS2D)) ) {
    setBkgPdf(getEffPdfMean(),getEffPdfSigma(),getEffPdfDist());
    change = true;
  }
  if (getEffPdfDist() != PDF::DIST_GAUS2D) {
    setEffPdfBkgCorr(0);
    change = true;
  }
  //
  if (getEffPdfDist()==PDF::DIST_NONE) {
    m_measurement.setEffPdfSigma(0.0);
    change = true;
  }
  if (getBkgPdfDist()==PDF::DIST_NONE) {
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
					getEffPdfDist(), getEffObs(), getEffPdfSigma(),
					getBkgPdfDist(), getBkgObs(), getBkgPdfSigma(), m_measurement.getNuisanceIntNorm());
  double up  = BeltEstimator::getSigUp( getNObserved(),
					getEffPdfDist(), getEffObs(), getEffPdfSigma(),
					getBkgPdfDist(), getBkgObs(), getBkgPdfSigma(), m_measurement.getNuisanceIntNorm());
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
//   dsLow  = (m_measurement.getEffObs() - m_effRangeInt.min())/m_measurement.getEffPdfPdfSigma();
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
  if ( getNObserved()>=0) {
    rval = BeltEstimator::getBeltMin( getNObserved(),
				      getEffPdfDist(),
                                      getEffObs(),
                                      getEffPdfSigma(),
                                      getBkgPdfDist(),
                                      getBkgObs(),
                                      getBkgPdfSigma(),
                                      getIntNorm());
    if (rval<20) rval=20; // becomes less reliable for small n
  }
  if (m_verbose>1) {
    std::cout << "Using max N(belt) = " << rval << std::endl;
    std::cout << "t(test var) = " << getSVar() << std::endl;
  }
  return rval;
}

void Pole::initBeltArrays() {
  //
  m_nBeltUsed = getNObserved();
  if (m_nBeltUsed<2) m_nBeltUsed=2;
  m_nBeltMinUsed = m_nBeltUsed;
  m_nBeltMaxUsed = 0;

  //  unsigned int nbs = static_cast<unsigned int>(m_nBelt);
  //  if (m_muProb.size()!=nbs) {
  m_muProb.resize(m_nBeltUsed,0.0);
  m_bestMuProb.resize(m_nBeltUsed,0.0);
  m_bestMu.resize(m_nBeltUsed,0.0);
  m_lhRatio.resize(m_nBeltUsed,0.0);
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
    int ntst = 1+int((sMax-sMin)/m_bestMuStep);
    double dmus=m_bestMuStep;
    if (ntst>m_bestMuNmax) {
      ntst=m_bestMuNmax;
      dmus=(sMax-sMin)/double(ntst-1);
    }
    //// TEMPORARY CODE - REMOVE /////
    //    ntst = 1000;
    //    m_bestMuStep = (sMax-sMin)/double(ntst);
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
	if (m_verbose>3) std::cout << " calcProb(" << n << ", " << mu_test << ") = " << lh_test << std::endl;
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
  for (int n=0; n<m_nBeltUsed; n++) {
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
  for (int n=0; n<m_nBeltUsed; n++) {
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
  nbMin = 0; // m_nBeltMinLast; // NOTE if using nBeltMinLast -> have to be able to approx prob for s in the range 0...nminlast
  nbMax = m_nBeltUsed-1;
  double norm_p_low=0;// m_sumPlowLast; // approx
  //  double norm_p_up=0;
  //
  //
  //  double sumPmin=0.0;
  //  double sumPtot=0.0;
  // find lower and upper
  bool expandedBelt=false;
  while (!upNfound) {
    if (n==m_nBeltUsed) { // expand nBelt
      expandedBelt=true;
      m_nBeltUsed++;
      nbMax++;
      m_muProb.push_back(0.0);
      m_lhRatio.push_back(0.0);
      if (usesFHC2()) {
        m_bestMu.push_back(0.0);
        m_bestMuProb.push_back(0.0);
        findBestMu(n);
      }
    }
    m_muProb[n] =  m_measurement.calcProb(n, s);
    if (m_verbose>0) {
      //      if ((s<5.02) &&(m_muProb[n]>0.00001)) std::cout << "s = " << s << "     : m_muProb[" << n << "] = " << m_muProb[n] << std::endl;
    }
    lhSbest = getLsbest(n);
    m_lhRatio[n]  = m_muProb[n]/lhSbest;
    if (m_verbose>8) {
      std::cout << "LHRATIO: " << "\t" << n << "\t" << s << "\t" << m_muProb[n] << "\t" << lhSbest << std::endl;
    }
    norm_p += m_muProb[n]; // needs to be renormalised
    //
    if ((!lowNfound) && (norm_p>m_minMuProb)) {
      lowNfound=true;
      nbMin = n;
      norm_p_low = norm_p;
    } else {
      if ((nInBelt>1) && lowNfound && (norm_p>1.0-m_minMuProb)) {
	upNfound = true;
	nbMax = n-1;
      }
    }
    if (lowNfound && (!upNfound)) nInBelt++;
    n++;
    //    std::cout << "calcLhRatio: " << n << ", p = " << m_muProb[n] << ", lhSbest = " << lhSbest << std::endl;
  }
  //
  m_nBeltMinLast = nbMin;
  if (nbMin<m_nBeltMinUsed) m_nBeltMinUsed = nbMin;
  if (nbMax>m_nBeltMaxUsed) m_nBeltMaxUsed = nbMax;
  //
  if (m_verbose>2) {
    if (expandedBelt) {
      std::cout << "calcLhRatio() : expanded belt to " << m_nBeltUsed << std::endl;
      std::cout << "                new size = " << m_muProb.size() << std::endl;
    }
  }
  return norm_p;//-norm_p_min;
}


int Pole::calcLimit(double s) {
  int k,i;
  int nBeltMinUsed;
  int nBeltMaxUsed;
  //
  // Reset belt probabilities
  //
  m_prevSumProb = m_sumProb;
  m_sumProb     = 0.0;
  //
  // Get RL(n,s)
  //
  if (m_verbose>2) {
    std::cout << std::endl;
    std::cout << "=== calcLimit( s = " << s << " )" << std::endl;
  }
  m_scanBeltNorm = calcLhRatio(s,nBeltMinUsed,nBeltMaxUsed);
  if (m_verbose>3) {
    std::cout << "--- Normalisation over n for s = " << s << " is " << m_scanBeltNorm << std::endl;
    if ((m_scanBeltNorm>1.5) || (m_scanBeltNorm<0.5)) {
      std::cout << "--- Normalisation off (" << m_scanBeltNorm << ") for s= " << s << std::endl;
      for (int n=0; n<nBeltMaxUsed; n++) {
	std::cout << "--- muProb[" << n << "] = " << m_muProb[n] << std::endl;
      }
    }
  }
  //
  // Get N(obs)
  //
  k = getNObserved();
  //
  // If k is outside the belt obtained above, the likelihood ratio is ~ 0 -> set it to 0
  //
  if ((k>nBeltMaxUsed) || (k<nBeltMinUsed)) {
    if (m_verbose>2) {
      std::cout << "--- Belt used for RL does not contain N(obs) - skip. Belt = [ "
                << nBeltMinUsed << " : " << nBeltMaxUsed << " ]" << std::endl;
    }
    m_lhRatio[k]   = 0.0;
    m_sumProb      = 1.0; // should always be >CL
    m_scanBeltNorm = 1.0;
    return 1; 
  }

  if (k>=m_nBeltUsed) {
    k=m_nBeltUsed; // WARNING::
    std::cout << "--- FATAL :: n_observed is larger than the maximum n used for R(n,s)!!" << std::endl;
    std::cout << "             -> increase nbelt such that it is more than n_obs = " << getNObserved() << std::endl;
    std::cout << "             *** IF THIS MESSAGE IS SEEN, IT'S A BUG!!! ***" << std::endl;
    exit(-1);
  }									\
  if (m_verbose>2) {
    std::cout << "--- Got nBelt range: " << nBeltMinUsed << ":" << nBeltMaxUsed << "( max = " << m_nBeltUsed-1 << " )" << std::endl;
    std::cout << "--- Will now make the likelihood-ratio ordering and calculate a probability" << std::endl;
  }


  // Calculate the probability for all n and the given s.
  // The Feldman-Cousins method dictates that for each n a
  // likelihood ratio (R) is calculated. The n's are ranked according
  // to this ratio. Values of n are included starting with that giving
  // the highest R and continuing with decreasing R until the total probability
  // matches the searched CL.
  // Below, the loop sums the probabilities for a given s and for all n with R>R0.
  // R0 is the likelihood ratio for n_observed.
  // When N(obs)=k lies on a confidencebelt curve (upper or lower), the likelihood ratio will be minimum.
  // If so, m_sumProb> m_cl
  i=nBeltMinUsed;
  bool done=false;
  if (m_verbose>9) std::cout << "\nSearch s: = " <<  s << std::endl;
  while (!done) {
    if(i != k) { 
      if(m_lhRatio[i] > m_lhRatio[k])  {
	m_sumProb  +=  m_muProb[i];
      }
      if (m_verbose>9) {
	std::cout << "RL[" << i << "] = " << m_lhRatio[i]
                  << ", RLmax[" << k << "] = " << m_lhRatio[k]
                  << ", sumP = " << m_sumProb << std::endl;
      }
    }
    i++;
    done = ((i>nBeltMaxUsed) || m_sumProb>m_cl); // CHANGE 11/8
  }
  //
  // Return:
  //  0 : if |diff|<0.001
  // +1 : if diff>0.001
  // -1 : if diff<-0.001
  //
  double diff  = (m_sumProb-m_cl)/(m_sumProb+m_cl);
  int rval;
  if (fabs(diff)<0.001) rval = 0;                // match!
  else                  rval = (diff>0 ? +1:-1); // need still more precision
  if (m_verbose>2) {
    std::cout << "--- Done. Results:" << std::endl;
    std::cout << "---    Normalisation  : " << m_scanBeltNorm << std::endl;
    std::cout << "---    Sum(prob)      : " << m_sumProb << std::endl;
    std::cout << "---    Direction      : " << rval      << std::endl;
    std::cout << "=== calcLimit() DONE!" << std::endl;
  }
  return rval;
}


bool Pole::limitsOK() {
  bool rval=false;
  if (m_lowerLimitFound && m_upperLimitFound) {
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

void Pole::calcPower() {
  double muTest;
  std::vector< std::vector<double> > fullConstruct;
  std::vector< double > probVec;
  std::vector< int > n1vec;
  std::vector< int > n2vec;
  int n1,n2;
  int nhyp = m_hypTest.n();
  double sumP;
  m_nBeltMinLast = 0;
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
      for ( n=n02+1; n<m_nBeltUsed; n++) {
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

void Pole::calcConstruct() {
  double mu_test;
  int i = 0;
  bool done = (i==m_hypTest.n());
  m_nBeltMinLast = 0;
  //
  while (!done) {
    mu_test = m_hypTest.min() + i*m_hypTest.step();
    calcConstruct(mu_test,true);
    i++;
    done = (i==m_hypTest.n()); // must loop over all hypothesis
  }
}

void Pole::calcNMin() { // calculates the minimum N rejecting s = 0.0
  m_nBeltMinLast = 0;
  m_rejs0P = calcBelt(0.0,m_rejs0N1,m_rejs0N2,false);//,-1.0);
}

void Pole::calcBelt() {
  double mu_test;
  int i = 0;
  bool done = (i==m_hypTest.n());
  //
  m_nBeltMinLast = 0;
  int n1,n2;
  while (!done) {
    mu_test = m_hypTest.min() + i*m_hypTest.step();
    calcBelt(mu_test,n1,n2,true);
    i++;
    done = (i==m_hypTest.n()); // must loop over all hypothesis
  }
}

bool Pole::calcLimits() {
  //  static bool heavyDBG=false;
  m_maxNorm = -1.0;
  m_lowerLimitFound = false;
  m_upperLimitFound = false;
  m_lowerLimit = 0;
  m_upperLimit = 0;
  m_lowerLimitNorm = 0;
  m_upperLimitNorm = 0;
  m_nBeltUsed = getNObserved();
  if (m_nBeltUsed<2) m_nBeltUsed=2;
  m_nBeltMinLast=0;
  m_nBeltMinUsed=m_nBeltUsed;
  m_nBeltMaxUsed=0;
  m_sumProb=0;
  m_prevSumProb=0;
  m_scanBeltNorm=0;
  //
  //
  //
  //
  // Obtain belt for s=0 - useful since lower limit is 0 if N(obs)<N2
  //
  calcNMin();
  //
  // If N(obs) is BELOW belt at s=0.0, then we have an empty belt.
  //
  if (getNObserved()<=m_rejs0N1) {
    m_lowerLimitFound = true;
    m_lowerLimitNorm  = 1.0;
    m_lowerLimit      = 0;
    m_upperLimitFound = true;
    m_upperLimitNorm  = 1.0;
    m_upperLimit      = 0;
    m_sumProb         = 0.0;
    return true;
  }
  double mutest;
  int    i;
  double p0;
  bool done = false;
  //
  if (m_verbose>0) std::cout << "   Scanning over all hypothesis" << std::endl;

  //
  // Find hypothesis start - s = (Nobs-bkg)/eff - should be in ~the middle of the belt
  //
  double mustart = getObservedSignal();
  double muhigh  = mustart;
  double mulow   = 0.0;
  double dmu;
  int dir;
  //
  calcLimit(mustart);
  p0 = m_sumProb;

  // check if N(obs) < N2(s=0.0)
  // if so, the lower limit is 0.0
  if (getNObserved()<m_rejs0N2) {
    m_lowerLimitFound = true;
    m_lowerLimitNorm  = 1.0;
    m_lowerLimit      = 0;
  } else {
    // First scan for lower limit
    done = false;
    m_prevSumProb = p0;
    i = 1;
    while (!done) {
      mutest = (muhigh+mulow)/2.0;
      dmu = fabs(muhigh-mulow)/(muhigh+mulow);
      //
      dir = calcLimit(mutest);
      if ((dir==0) || (dmu<0.001)) { // hardcoded precision CHANGE! same as the one deciding for dir==0
        done = true;
        m_lowerLimitFound = true;
        m_lowerLimitNorm  = m_scanBeltNorm;
        m_lowerLimit      = mutest;
      } if (dir==1) {
        mulow = mutest;
      } else {
        muhigh = mutest;
      }
      i++;
    }
  }

  // Scan for upper limit
  done = false;
  m_prevSumProb = p0;
  i = 1;

  // find a rough range guaranteed to cover the upper limit
  dir = -1;
  mulow  = mustart;
  muhigh = 2.0*(mustart+1.0)-m_lowerLimit; // pick an initial upper limit
  // loop
  while (dir==-1) {
    dir = calcLimit(muhigh); // what if out of reach for nbelt???? CHECK
    if (dir==-1) muhigh *= 1.1; // scale up
  }
  // found a rough upper limit, now go down
  while (!done) {
    mutest = (mulow+muhigh)/2.0;
    dmu = fabs(muhigh-mulow)/(muhigh+mulow);
    dir = calcLimit(mutest);
    if ((dir==0) || (dmu<0.001)) { // see comment above for lower limit
      done = true;
      m_upperLimitFound = true;
      m_upperLimitNorm  = m_scanBeltNorm;
      m_upperLimit      = mutest;
    } if (dir==1) {
      muhigh = mutest;
    } else {
      mulow = mutest;
    }
    i++;
  }
  //
  return limitsOK();
}

// bool Pole::calcLimits() {
//   //  static bool heavyDBG=false;
//   m_maxNorm = -1.0;
//   m_lowerLimitFound = false;
//   m_upperLimitFound = false;
//   m_lowerLimit = 0;
//   m_upperLimit = 0;
//   m_lowerLimitNorm = 0;
//   m_upperLimitNorm = 0;
//   m_nBeltMinLast=0;
//   m_nBeltMinUsed=m_nBelt;
//   m_nBeltMaxUsed=0;
//   m_sumProb=0;
//   m_prevSumProb=0;
//   m_scanBeltNorm=0;
//   //
//   // Obtain belt for s=0
//   //
//   calcNMin();
//   //
//   // If N(obs) is outside belt, fail.
//   //
//   if (getNObserved()>=m_nBelt) return false;
//   double mutest;
//   int    i;
//   double p0;
//   bool done = false;
//   //
//   if (m_verbose>0) std::cout << "   Scanning over all hypothesis" << std::endl;

//   //
//   // Find hypothesis start - s = (Nobs-bkg)/eff - should be in ~the middle of the belt
//   //
//   double mustart = getObservedSignal();

//   if (calcLimit(mustart,true)) {
//     p0 = m_sumProb;
//   } else {
//     std::cerr << "FATAL: scan with default s0 failed! BUG!?!" << std::endl;
//     return false;
//   }

//   // First scan for lower limit
//   done = false;
//   m_prevSumProb = p0;
//   i = 1;
//   while (!done) {
//     mutest = mustart - i*m_limitHypStep;
//     if (calcLimit(mutest,true)) {
//       done = m_lowerLimitFound;
//     } else {
//       done = true;
//     }
//     i++;
//   }

//   // Scan for upper limit
//   done = false;
//   m_prevSumProb = p0;
//   i = 1;
//   while (!done) {
//     mutest = mustart + i*m_limitHypStep;
//     if (calcLimit(mutest,false)) {
//       done = m_upperLimitFound;
//     } else {
//       done = true;
//     }
//     i++;
//   }
//   return limitsOK();
// }

bool Pole::calcCoverageLimits() {
  m_maxNorm = -1.0;
  m_lowerLimitFound = false;
  m_upperLimitFound = false;
  m_lowerLimit = 0;
  m_upperLimit = 0;
  m_lowerLimitNorm = 0;
  m_upperLimitNorm = 0;
  m_nBeltMinLast=0;
  m_nBeltMinUsed=m_nBelt;
  m_nBeltMaxUsed=0;
  m_sumProb=0;
  m_prevSumProb=0;
  m_scanBeltNorm=0;

  int dir = calcLimit(m_sTrue);
  return (dir==-1); // s(true) inside belt - p(s)<cl
}

// bool Pole::calcCoverageLimitsOLD() {
//   //
//   // Same as calcLimits() but it will stop the scan as soon it is clear
//   // that the true s (m_sTrueMean) is inside or outside the range.
//   //
//   bool decided = false;
//   m_maxNorm = -1.0;
//   m_lowerLimitFound = false;
//   m_upperLimitFound = false;
//   m_lowerLimit = 0;
//   m_upperLimit = 0;
//   m_lowerLimitNorm = 0;
//   m_upperLimitNorm = 0;
//   m_nBeltMinLast=0;
//   m_nBeltMinUsed=m_nBelt;
//   m_nBeltMaxUsed=0;
//   m_sumProb=0;
//   m_prevSumProb=0;
//   m_scanBeltNorm=0;

//   //
//   // If N(obs) is outside belt, fail.
//   //
//   if (getNObserved()>=m_nBelt) return false;
//   double mutest;
//   int    i;
//   double p0;
//   bool done = false;
//   //
//   if (m_verbose>0) std::cout << "   Scanning over all hypothesis" << std::endl;

//   //
//   // Find hypothesis start - s = (Nobs-bkg)/eff - should be in ~the middle of the belt
//   //
//   double mustart = getObservedSignal();

//   if (calcLimit(mustart,true)) {
//     p0 = m_sumProb;
//   } else {
//     std::cerr << "FATAL: scan with default s0 failed! BUG!?!" << std::endl;
//     return false;
//   }

//   //
//   // Decide where to scan - below or above mustart
//   //
//   bool inLowerHalf = (m_sTrue<=mustart);

//   if (inLowerHalf) {
//     // The truth is in the lower half - scan for lower limit
//     done = false;
//     m_prevSumProb = p0;
//     i = 1;
//     while (!done) {
//       mutest = mustart - i*m_limitHypStep;
//       if (calcLimit(mutest,true)) {
//         done = m_lowerLimitFound;
//         if (mutest<m_sTrue) { // test if s(true) is within the range
//           decided = true;
//           done    = true;
//         }
//       } else {
//         done = true;
//       }
//       i++;
//     }
//   } else {
//     // The truth is in the upper half - scan for upper limit
//     done = false;
//     m_prevSumProb = p0;
//     i = 1;
//     while (!done) {
//       mutest = mustart + i*m_limitHypStep;
//       if (calcLimit(mutest,false)) {
//         done = m_upperLimitFound;
//         if (mutest>m_sTrue) { // test if s(true) is within the range
//           decided = true;
//           done    = true;
//         }
//       } else {
//         done = true;
//       }
//       i++;
//     }
//   }

//   return decided;
// }

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
    calcBelt();
  }
  if (m_verbose>0) std::cout << "Calculating limit" << std::endl;
  if (m_coverage) {
    rval=calcCoverageLimits();
  } else {
    rval=calcLimits();
  }
  // Should not do this - if probability is OK then the belt is also OK...?
  // The max N(Belt) is defined by a cutoff in probability (very small)
  //  if (m_nBeltMaxUsed==m_nBelt) rval=false; // reject limit if the full belt is used
  return rval;
}

void Pole::printLimit(bool doTitle) {
  std::string cmtPre;
  std::string linePre;
  bool simple;
  switch (m_printLimStyle) {
  case 0:
    cmtPre = "#";
    linePre = "  ";
    simple=false;
    break;
  case 1:
    cmtPre = "";
    linePre = " ";
    simple=false;
    break;
  case 2:
    cmtPre = "";
    linePre = " ";
    simple=true;
    break;
  default:
    cmtPre = "";
    linePre = " ";
    simple=true;
  }
  if (doTitle && (!simple)) {
    std::cout << cmtPre << "-------------------------------------------------------------------------------------" << std::endl;
    std::cout << cmtPre << " Max N(belt) set  : " << m_nBelt << std::endl;
    std::cout << cmtPre << " Max N(belt) used : " << m_nBeltMaxUsed << std::endl;
    std::cout << cmtPre << " Min N(belt) used : " << m_nBeltMinUsed << std::endl;
    std::cout << cmtPre << " Prob(belt s=0)   : " << m_rejs0P << std::endl;
    std::cout << cmtPre << " N1(belt s=0)     : " << m_rejs0N1 << std::endl;
    std::cout << cmtPre << " N2(belt s=0)     : " << m_rejs0N2 << std::endl;
    
    std::cout << cmtPre << "-------------------------------------------------------------------------------------" << std::endl;
    std::cout << cmtPre << " N(obs)  eff(mean)     eff(sigma)      bkg(mean)       bkg(sigma)      low     up" << std::endl;
    std::cout << cmtPre << "-------------------------------------------------------------------------------------" << std::endl;
  }
  std::cout << linePre;
  if (!simple) {
    TOOLS::coutFixed(6,getNObserved());   std::cout << "\t  ";
    TOOLS::coutFixed(6,getEffObs());      std::cout << "\t";
    TOOLS::coutFixed(6,getEffPdfSigma()); std::cout << "\t";
    TOOLS::coutFixed(6,getBkgObs());      std::cout << "\t";
    TOOLS::coutFixed(6,getBkgPdfSigma()); std::cout << "\t";
    TOOLS::coutFixed(2,m_scaleLimit*m_lowerLimit);     std::cout << "\t";
    TOOLS::coutFixed(2,m_scaleLimit*m_upperLimit);     std::cout << std::endl;
  }
  if (simple) {
    std::cout << "Limits = [ ";
    TOOLS::coutFixed(2,m_scaleLimit*m_lowerLimit); std::cout << ", ";
    TOOLS::coutFixed(2,m_scaleLimit*m_upperLimit); std::cout << " ]";
    std::cout << std::endl;
  }
}

void Pole::printSetup() {
  std::cout << "\n";
  std::cout << "================ P O L E ==================\n";
  std::cout << " Confidence level   : " << m_cl << std::endl;
  std::cout << " N observed         : " << getNObserved() << std::endl;
  std::cout << "----------------------------------------------\n";
  std::cout << " Coverage friendly  : " << yesNo(m_coverage) << std::endl;
  std::cout << " True signal        : " << m_sTrue << std::endl;
  std::cout << "----------------------------------------------\n";
  std::cout << " Efficiency meas    : " << getEffObs() << std::endl;
  std::cout << " Efficiency sigma   : " << getEffPdfSigma() << std::endl;
  std::cout << " Efficiency dist    : " << PDF::distTypeStr(getEffPdfDist()) << std::endl;
  std::cout << "----------------------------------------------\n";
  std::cout << " Background meas    : " << getBkgObs() << std::endl;
  std::cout << " Background sigma   : " << getBkgPdfSigma() << std::endl;
  std::cout << " Background dist    : " << PDF::distTypeStr(getBkgPdfDist()) << std::endl;
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
  std::cout << " Limit hyp. step    : " << m_limitHypStep << std::endl;
  std::cout << "----------------------------------------------\n";
  std::cout << " *Test hyp. min     : " << m_hypTest.min() << std::endl;
  std::cout << " *Test hyp. max     : " << m_hypTest.max() << std::endl;
  std::cout << " *Test hyp. step    : " << m_hypTest.step() << std::endl;
  std::cout << "----------------------------------------------\n";
  std::cout << " Min. probability   : " << m_minMuProb << std::endl;
  std::cout << "----------------------------------------------\n";
  if (m_method==RL_FHC2) {
    std::cout << " Step mu_best       : " << m_bestMuStep << std::endl;
    std::cout << " Max N, mu_best     : " << m_bestMuNmax << std::endl;
    std::cout << "----------------------------------------------\n";
  }
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
  std::cout << "----------------------------------------------\n";
  std::cout << " Parameters prefixed with a * above are not   \n";
  std::cout << " relevant for limit calculations.\n";
  std::cout << "==============================================" << std::endl;
  
  //
}

void Pole::printFailureMsg() {
  std::cout << "ERROR: limit calculation failed. Possible causes:" << std::endl;
  std::cout << "1. nbelt is too small (set nbelt = " << getNBelt() << ", max used = " << getNBeltMaxUsed() << ")" << std::endl;
  std::cout << "2. precision in integrations (eff,bkg) not sufficient" << std::endl;
  std::cout << "3. limit hypothesis step size too large ( step = " << m_limitHypStep << " )" << std::endl;
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

