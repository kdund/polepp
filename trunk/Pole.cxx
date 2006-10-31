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
#include <algorithm>

#include "Pole.h"


void Pole::initDefault() {
  m_thresholdBS    = 0.001;
  m_thresholdPrec = 0.01;
  m_scaleLimit=1.0;
  m_inputFile="";
  m_inputFileLines=0;
  m_printLimStyle=0; // computer readable
  m_verbose=0;
  m_coverage = false;
  m_method = RL_FHC2;
  m_normMaxDiff = 0.01;
  m_lowerLimitNorm = -1.0;
  m_upperLimitNorm = -1.0;
  m_lowerLimitPrec = -1.0;
  m_upperLimitPrec = -1.0;
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
  setTestHyp(0.0,-1.0,0.01);
  //
  m_validBestMu = false;
  m_nBeltUsed    = 5; // not really needed to be set here - set in initBeltArrays()
  m_nBeltMinUsed = m_nBeltUsed; // idem
  m_nBeltMaxUsed = 0;
}

void Pole::execute() {
  if (m_inputFile.size()>0) {
    exeFromFile();
  } else {
    m_printLimStyle = 2;
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
  double y,z,t,u;
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
      //      sstr >> n >> ed >> em >> es >> bd >> bm >> bs;
      z=0;
      sstr >> n >> y >> z >> t >> u;
      //
      if ((sstr) && (z>0)) {
        nread++;
        if ((m_verbose>0) && (nlines<10)) {
          std::cout << "LINE: ";
          //          std::cout << n << "\t" << ed << "\t" << em << "\t" << es << "\t" << bd << "\t" << bm << "\t" << bs << std::endl;
          std::cout << n << "\t" << y << "\t" << z << "\t" << t << "\t" << u  << std::endl;
        }
        m_measurement.setEffScale(1.0/u);
        m_measurement.setBkgScale(1.0/t);
        setNObserved(static_cast<int>(n+0.5));
        setEffPdf( z, z, PDF::DIST_POIS );
        setEffObs();
        setBkgPdf( y, y, PDF::DIST_POIS );
        setBkgObs();
        //        setScaleLimit(u/t);

//         setNObserved(n);
//         setEffPdf( em, es, PDF::DISTYPE(ed) );
//         setEffObs();
//         setBkgPdf( bm, bs, PDF::DISTYPE(bd) );
//         setBkgObs();
//         setScaleLimit(1.0);
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

// void Pole::setTestHyp(double step) {
//   std::cout << "ERROR: Should not be called! setTestHyp(step)" << std::endl;
//   exit(-1);
//   //
//   // Find hypothesis test range based on the input measurement
//   // MUST be called after the measurement is set.
//   //
//   if (step<=0) {
//     step = m_hypTest.step();
//     if (step<=0) step = 0.01;
//   }
  
//   double low = BeltEstimator::getSigLow(getNObserved(),
// 					getEffPdfDist(), getEffObs(), getEffPdfSigma(),
// 					getBkgPdfDist(), getBkgObs(), getBkgPdfSigma(), m_measurement.getNuisanceIntNorm());
//   double up  = BeltEstimator::getSigUp( getNObserved(),
// 					getEffPdfDist(), getEffObs(), getEffPdfSigma(),
// 					getBkgPdfDist(), getBkgObs(), getBkgPdfSigma(), m_measurement.getNuisanceIntNorm());
//   m_hypTest.setRange(low,up,step);
// }

void Pole::setTestHyp(double low, double high, double step) {
  //
  // Set explicitely the test range.
  // * step<=0  => step = (high-low)/100
  // * high<low => set high = max( observed signal, 1.0 )
  //
  if (high<low) {
    low = 0.0;
    high = 1.0 + getObservedSignal();
  }
  if (step<=0) step = (high-low)/100.0; // default 100 pts
  m_hypTest.setRange(low,high,step);
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

void Pole::initBeltArrays() {
  //
  m_nBeltUsed = getNObserved();
  if (m_nBeltUsed<2) m_nBeltUsed=2;
  m_nBeltMinUsed = m_nBeltUsed;
  m_nBeltMaxUsed = 0;
  //
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
  if(n<getBkgObs()*getBkgScale()) {
    m_bestMu[n] = 0; // best mu is 0
    m_bestMuProb[n] = m_measurement.calcProb(n,0);
  } else {
    //    sMax = double(n)-m_measurement.getBkgObs(); // OLD version
    sMax = (double(n) - getBkgIntMin()*getBkgScale())/(getEffObs()*getEffScale());
    if(sMax<0) {sMax = 0.0;}
    //    sMin = (double(n) - m_bkgRangeInt.max())/m_effRangeInt.max();
    //    sMin = sMax/2.0; //
    sMin = (double(n) - getBkgIntMax()*getBkgScale())/(getEffIntMax()*getEffScale());
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
    if (m_verbose>1) std::cout << "FindBestMu range: "
	  << " I(bkg) " << getBkgIntMax()
	  << " I(eff)" << getEffIntMax()
	  << " N " << n
	  << " Bkg " << m_measurement.getBkgObs()
	  << " ntst " << ntst << " [" << sMin << "," << sMax << "] => ";
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
    for (unsigned int i=0; i<(m_bestMu.size()<10 ? m_bestMu.size():10); i++) {
      std::cout << m_bestMu[i] << "\t" << m_bestMuProb[i] << std::endl;
    }
  }
  m_validBestMu = true;
}

double Pole::calcLhRatio(double s, int & nbMin, int & nbMax) {
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
    //
    lhSbest = getLsbest(n);
    if (lhSbest>0) {
       m_lhRatio[n]  = m_muProb[n]/lhSbest;
    } else {
       m_lhRatio[n]  = 0.0;
    }
    if (m_verbose>8) {
      std::cout << "LHRATIO: " << "\t" << n << "\t"
                << s << "\t"
                << m_lhRatio[n] << "\t"
                << m_muProb[n] << "\t"
                << getSbest(n) << "\t"
                << lhSbest
                << std::endl;
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
    std::cout << "LHRATIOSUM: " << norm_p << std::endl; 
    if (expandedBelt) {
      std::cout << "calcLhRatio() : expanded belt to " << m_nBeltUsed << std::endl;
      std::cout << "                new size = " << m_muProb.size() << std::endl;
    }
  }
  return norm_p;//-norm_p_min;
}


int Pole::calcLimit(double s, double & prec) {
  int k,i;
  int nBeltMinUsed;
  int nBeltMaxUsed;
  //
  // Reset belt probabilities
  //
  m_sumProb     = 0.0;
  //
  // Get RL(n,s)
  //
  m_scanBeltNorm = calcLhRatio(s,nBeltMinUsed,nBeltMaxUsed);
  if (m_verbose>2) {
    std::cout << std::endl;
    std::cout << "=== calcLimit( s = " << s << " )" << std::endl;
    std::cout << "--- Normalisation over n for s  is " << m_scanBeltNorm << std::endl;
  }
  if (m_verbose>3) {
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
    prec           = 1.0;
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
  int nsum=0;
  if (m_verbose>9) std::cout << "\nSearch s: = " <<  s << std::endl;
  int nbmin,nbmax;
  nbmin=100000;
  nbmax=0;
  while (!done) {
    if(i != k) { 
      if(m_lhRatio[i] > m_lhRatio[k])  {
        if (i<nbmin) nbmin=i;
        if (i>nbmax) nbmax=i;
        m_sumProb  +=  m_muProb[i];
        nsum++;
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
  if (nsum==0) nbmin=0;
  //
  // Return:
  //  0 : if |diff|<0.001
  // +1 : if diff>0.001
  // -1 : if diff<-0.001
  //
  double diff = (m_sumProb-m_cl)/m_cl; //
  //  double diffOpt = 1.0 - ((1.0-m_sumProb)/(1.0-m_cl));
  int rval;
  if (fabs(diff)<m_thresholdPrec) rval = 0;                // match!
  else                             rval = (diff>0 ? +1:-1); // need still more precision
  prec = (nsum>0 ? fabs(diff) : -1.0 );
  if (m_verbose>3) {
    std::cout << "s = " << s
              << " ; nsum = " << nsum
              << " ; nbmin = " << nbmin
              << " ; nbmax = " << nbmax
              << " ; Sum(prob) = " << m_sumProb
              << " ; Dir = " << rval << " prec = " << prec << std::endl;
  }
  if (m_verbose>2) {
    std::cout << "--- Done. Results:" << std::endl;
    std::cout << "---    Normalisation  : " << m_scanBeltNorm << std::endl;
    std::cout << "---    Sum(prob)      : " << m_sumProb << std::endl;
    std::cout << "---    Diff           : " << diff      << std::endl;
    std::cout << "---    Direction      : " << rval      << std::endl;
    std::cout << "=== calcLimit() DONE!" << std::endl;
  }
  return rval;
}


bool Pole::limitsOK() {
  bool rval=false;
  if (m_lowerLimitFound && m_upperLimitFound) {
    rval = ( (normOK(m_lowerLimitNorm) && normOK(m_upperLimitNorm)) &&
             (m_lowerLimitPrec<m_thresholdPrec) &&
             (m_upperLimitPrec<m_thresholdPrec) );
  }
  return rval;
}

//*********************************************************************//
//*********************************************************************//
//*********************************************************************//

void Pole::calcConstruct(double s, bool title) {
  int i;
  int nb1,nb2;
  //
  double norm_p = calcLhRatio(s,nb1,nb2);
  if (title) {
    std::cout << "#========================================================================" << std::endl;
    std::cout << "#          Signal       N       R(L)            P(s)            Norm     " << std::endl;
    std::cout << "#========================================================================" << std::endl;
  }
  for (i=nb1; i<nb2; i++) {
    std::cout << "CONSTRUCT: ";
    TOOLS::coutFixed(4,s);
    std::cout << "\t";
    TOOLS::coutFixed(6,i);
    std::cout << "\t";
    TOOLS::coutFixed(6,m_lhRatio[i]);
    std::cout << "\t";
    TOOLS::coutFixed(6,m_muProb[i]);
    std::cout << "\t";
    TOOLS::coutFixed(6,norm_p);
    std::cout << std::endl;
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

double Pole::calcBelt(double s, int & n1, int & n2, bool verb, bool title) { //, double muMinProb) {
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
    if (title) {
      std::cout << "#==================================================" << std::endl;
      std::cout << "#         Signal       N1      N2      P(n1,n2|s)  " << std::endl;
      std::cout << "#==================================================" << std::endl;
    }
    std::cout << "CONFBELT: ";
    TOOLS::coutFixed(4,s);
    std::cout << "\t";
    TOOLS::coutFixed(6,n1);
    std::cout << "\t";
    TOOLS::coutFixed(6,n2);
    std::cout << "\t";
    TOOLS::coutFixed(6,sumProb);
    std::cout << std::endl;
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
    sumP = calcBelt(muTest,n1,n2,false,false);
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
  int i1 = int(getTrueSignal()/m_hypTest.step());
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
  bool title = true;
  m_nBeltMinLast = 0;
  //
  while (!done) {
    mu_test = m_hypTest.min() + i*m_hypTest.step();
    calcConstruct(mu_test,title);
    i++;
    title = false;
    done = (i==m_hypTest.n()); // must loop over all hypothesis
  }
}

void Pole::calcNMin() { // calculates the minimum N rejecting s = 0.0
  m_nBeltMinLast = 0;
  m_rejs0P = calcBelt(0.0,m_rejs0N1,m_rejs0N2,false,false);//,-1.0);
}

void Pole::calcBelt() {
  double mu_test;
  int i = 0;
  bool done = (i==m_hypTest.n());
  bool title = true;
  //
  m_nBeltMinLast = 0;
  int n1,n2;
  while (!done) {
    mu_test = m_hypTest.min() + i*m_hypTest.step();
    calcBelt(mu_test,n1,n2,true,title);
    title = false;
    i++;
    done = (i==m_hypTest.n()); // must loop over all hypothesis
  }
}

bool Pole::scanLowerLimit( double mustart, double p0 ) {
  //
  //------------------------------
  // Start search for lower limit
  //------------------------------
  // This routine will scan for a lower limit starting at mustart.
  //  < p0 = m_sumProb after a call to calcLimit(mustart). >
  //
  // Special case:
  //  check if N(obs) < N2(s=0.0)
  //  if so, the lower limit is 0.0
  //
  if (getNObserved()<m_rejs0N2) {
    m_lowerLimitFound = true;
    m_lowerLimitNorm  = 1.0;
    m_lowerLimit      = 0;
    if (m_verbose>1) {
      std::cout << "*** Lower limit scan: limit concluded from the fact that N(obs) = " << getNObserved()
                << " and Belt(s=0) = [ " << m_rejs0N1 << " : " << m_rejs0N2 << " ]" << std::endl;
    }
  } else {
    //
    // Scan for lower limit
    //
    bool done = false;
    int  dir = -1;
    //
    // Special case:
    //  N(obs) == N2(s=0)
    //  if (dir!=1) then the lower limit is at s=0
    //  if (dir==1) m_sumProb>cl => lower limit at some s>0
    //
    if (getNObserved()==m_rejs0N2) {
      if (m_verbose>1)
        std::cout << "*** Lower limit scan: special case -> N(obs) = N2(s=0), test for limit at s=0" << std::endl;
      dir = calcLimit(0);
      if (m_verbose>1)
        std::cout << "*** Lower limit scan: calcLimit at sL = 0 ==> P = " << m_sumProb
                  << " and dir = " << dir
                  << std::endl;
      if (dir<1) {
        done = true;
        m_lowerLimitFound = true;
        m_lowerLimitNorm  = m_scanBeltNorm;
        m_lowerLimit = 0;
        if (m_verbose>1)
          std::cout << "*** Lower limit scan: limit concluded to be 0"
                    << std::endl;
      }
    }
    //
    // Lower limit binary search
    //
    int    cn = 0;                // reset loop counter
    double mutestPrev = mustart;  // start point
    double muhigh     = mustart;  // binary search, high point
    double mulow      = 0.0;      // idem, low point
    double prevPrec   = 100000.0; // a dummy large value
    double precMin    = 100000.0;
    double prec       = 0;
    double dmu        = 0;
    double mutest;
    //
    if ((m_verbose>1) && (!done))
      std::cout << "*** Lower limit scan: binary search, start at mutest = " << mustart << std::endl;
    while (!done) {
      mutest = (muhigh+mulow)/2.0;
      dmu = 2.0*fabs(mutestPrev-mutest)/(mutestPrev+mutest);
      //
      if (m_verbose>1)
        std::cout << "*** Lower limit scan: test at mutest  = " << mutest
                  << " mulow,muhigh = " << mulow << " , " << muhigh
                  << std::endl;
      //
      dir = calcLimit(mutest,prec);
      if (!(prec<0)) {
        if (prec<precMin) {
          precMin           = prec;
          m_lowerLimitNorm  = m_scanBeltNorm;
          m_lowerLimit      = mutest;
          m_lowerLimitPrec  = prec;
        }
      }
      if (m_verbose>1)
        std::cout << "*** Lower limit scan: got dir = " << dir
                  << ", P(s) = " << m_sumProb
                  << " and prec = " << prec << std::endl;
      if ((dir==0) || (dmu<m_thresholdBS)) {
        if (m_verbose>1) {
          if (dir==0) {
            std::cout << "*** Lower limit scan: scan stopped since sufficient precision is reached!"
                      << std::endl;
            std::cout << "*** Lower limit scan: precision = " << prec << " < " << m_thresholdPrec
                      << std::endl;
          } else if (dmu<m_thresholdBS) {
            std::cout << "*** Lower limit scan: scan stopped since change in mutest is small enough!"
                      << std::endl;
            std::cout << "*** Lower limit scan: delta mu = " << dmu << " < " << m_thresholdBS
                      << std::endl;
          }
        }
        done = true;
        m_lowerLimitFound = true;
      } else if (dir==1) {
        mulow = mutest;
      } else {
        muhigh = mutest;
      }
      mutestPrev = mutest;
      prevPrec = prec;
      cn++;
      if (cn>1000) {
        done=true;
        std::cout << "*** WARNING calcLimit() failed when scanning for lower limit! Infinite or long loop encountered - limit at 1000" << std::endl;
      }
    }
    if ((m_lowerLimitFound) && (m_verbose>1)) {
      std::cout << "*** Lower limit scan: limit from full scan; dir = " << dir << " dmu = " << dmu <<  " and sL = " << m_lowerLimit << std::endl;
    }
  }
  return m_lowerLimitFound;
}


bool Pole::scanUpperLimit( double mustart, double p0 ) {
  //
  if (!m_lowerLimitFound) return false;
  //
  // Scan for upper limit
  //
  const double muscale=1.1;
  bool   done     = false;
  int    dir      = -1;
  double mulow    = mustart;
  double muhigh   = 2.0*(mustart+1.0)-m_lowerLimit; // pick an initial upper limit
  double usescale = muscale;
  double prevmu   = muhigh;
  double prec     = 0;
  double prevPrec = 100000.0;
  double precMin  = 100000.0;
  int    cn       = 0;
  double mutest;
  double muRoughLim;
  double dmu;
  //
  // Start at muhigh.
  // calcLimit returns:
  //---------------------------------------------
  // dir | prec | action
  //---------------------------------------------
  //  -  | <0   | sumProb==0, muhigh at s_best
  // +1  | >0   | sumProb>cl
  // -1  | <0   | sumProb<cl
  //  0  |  0   | sumProb = cl (within precision)
  //---------------------------------------------
  if (m_verbose>1)
    std::cout << "*** Upper limit scan: find a rough limit, starting at " << muhigh << std::endl;
  mutest = muhigh;
  muRoughLim = mutest;
  while (!done) { //dir==-1) {
    dir = calcLimit(mutest,prec);
    if (m_verbose>1)
      std::cout << "*** Upper limit scan: test at " << mutest
                << " gave dir = " << dir
                << " , precision = " << prec
                << " and prob = " << m_sumProb
                << std::endl;
    if (prec<0) {
      usescale = muscale;
    } else {
      if (cn>0) { // check direction
        if (prec>prevPrec) { // if precision is degraded => we're going in the wrong direction
          usescale = 1.0/usescale;
        }
      }
      if (prec<precMin) {
        precMin = prec;
        muRoughLim = mutest;
      }
      done = false;
      // check the precision
      if (cn>1) {
        done = (prec>prevPrec);
      }
      // step up counter of ok mu
      cn++;
      prevPrec = prec;
      // check for infinite loop
      if ((!done) && (cn>100)) {
        done = true;
        std::cout << "*** WARNING infinite (?) loop when scanning for the rough upper limit!" << std::endl;
      }
    }
    // continue if not done
    if (!done) {
      prevmu = mutest;
      mutest *= usescale;
    }
  }
  muhigh = muRoughLim;
  //
  // found a rough upper limit, now go down
  // binary search start range is [ mulow = mustart, muhigh ]
  //
  if (m_verbose>1) {
    std::cout << "*** Upper limit scan: found a rough limit at " << muhigh << std::endl;
    std::cout << "*** Upper limit scan: binary search" << std::endl;
  }
  cn = 0;
  done = false;
  prevmu = mustart;  // start point
  precMin = 100000.0;
  while (!done) {
    mutest = (mulow+muhigh)/2.0;
    if (m_verbose>1)
      std::cout << "*** Upper limit scan: test at mutest  = " << mutest
                << " mulow,muhigh = " << mulow << " , " << muhigh
                << std::endl;
    //      dmu = fabs(muhigh-mulow)/(muhigh+mulow);
    dmu = 2.0*fabs(prevmu-mutest)/(prevmu+mutest);
    dir = calcLimit(mutest,prec);
    if (!(prec<0)) {
      if (prec<precMin) {
        precMin           = prec;
        m_upperLimitNorm  = m_scanBeltNorm;
        m_upperLimit      = mutest;
        m_upperLimitPrec  = prec;
      }
    }
    if (m_verbose>1)
      std::cout << "*** Upper limit scan: got dir = " << dir
                << " , P(s) = " << m_sumProb
                << " and prec = " << prec
                << std::endl;
    if ((dir==0) || (dmu<m_thresholdBS)) { // see comment above for lower limit
      if (m_verbose>1) {
        if (dir==0) {
          std::cout << "*** Upper limit scan: scan stopped since sufficient precision is reached!"
                    << std::endl;
          std::cout << "*** Upper limit scan: precision = " << prec << " < " << m_thresholdPrec
                    << std::endl;
        } else if (dmu<m_thresholdBS) {
          std::cout << "*** Upper limit scan: scan stopped since change in mutest is small enough!"
                    << std::endl;
          std::cout << "*** Upper limit scan: delta mu = " << dmu << " < " << m_thresholdBS
                    << std::endl;
        }
      }
      done = true;
      m_upperLimitFound = true;
    } else if (dir==1) {
      muhigh = mutest;
    } else {
      mulow = mutest;
    }
    prevmu   = mutest;
    prevPrec = prec;
    cn++;
    if (cn>1000) { // TODO: change this to a more smart check
      done=true;
      std::cout << "calcLimit() failed when scanning for UL! Infinite or long loop encountered - limit at 1000" << std::endl;
    }
  }
  if ((m_upperLimitFound) && (m_verbose>1)) {
    std::cout << "*** Upper limit scan: limit obtained from full scan; dir = " << dir << " dmu = " << dmu <<  " and sU = " << m_upperLimit << std::endl;
  }
  return m_upperLimitFound;
}

bool Pole::calcLimit() {
  if (m_verbose>1) std::cout << "*** Calculating limits" << std::endl;
  //
  // Reset
  //
  resetCalcLimit();
  //
  // Obtain belt for s=0 - useful since lower limit is 0 if N(obs)<N2
  //
  calcNMin();
  //
  // If N(obs) is BELOW belt at s=0.0, then we have an empty belt.
  //
  if (getNObserved()<m_rejs0N1) {
    m_lowerLimitFound = true;
    m_lowerLimitNorm  = 1.0;
    m_lowerLimit      = 0;
    m_lowerLimitPrec  = 0   ;
    m_upperLimitFound = true;
    m_upperLimitNorm  = 1.0;
    m_upperLimitPrec  = 0   ;
    m_upperLimit      = 0;
    m_sumProb         = 0.0;
    m_coversTruth     = false; // SHOULD THIS BE TRUE IF s=0??? (getTrueSignal()==0);
    std::cout << "WARNING: Empty limit created!" << std::endl;
    return true;
  }
  //
  // Find hypothesis start - s = (Nobs-bkg)/eff - should be in ~the middle of the belt
  //
  //
  if (m_verbose>1) std::cout << "*** Scanning over all hypothesis" << std::endl;
  //
  double mustart = getObservedSignal();
  //
  if (m_verbose>1)
    std::cout << "*** Calculating belt for start s = " << mustart << std::endl;
  //
  // Calculate belt for mustart
  //
  calcLimit(mustart);
  //
  if (m_verbose>1) {
    std::cout << "*** Obtained probability p = " << m_sumProb << std::endl;
    std::cout << "*** N2(s=0)                = " << m_rejs0N2 << std::endl;
  }

  scanLowerLimit(mustart,m_sumProb);

  scanUpperLimit(mustart,m_sumProb);

  //
  bool rval=limitsOK();
  if (rval) {
    m_coversTruth = ((getTrueSignal()>=m_lowerLimit) && (getTrueSignal()<=m_upperLimit));
    //    calcLimit(getTrueSignal());
//     std::cout << "COV: ";
//     TOOLS::coutFixed("s = ",4,getTrueSignal()); std::cout << "   ";
//     TOOLS::coutFixed("p = ",4,m_sumProb); std::cout << "   ";
//     TOOLS::coutFixed("l = ",4,m_lowerLimit); std::cout << "   ";
//     TOOLS::coutFixed("u = ",4,m_upperLimit); std::cout << "   ";
//     std::cout << "c = " << TOOLS::yesNo(m_coversTruth) << std::endl;
  }
  return rval;
}

////////////////////////////// OLD 
// bool Pole::calcLimit() {
//   if (m_verbose>1) std::cout << "*** Calculating limits" << std::endl;
//   resetCalcLimit();
//   //
//   // Obtain belt for s=0 - useful since lower limit is 0 if N(obs)<N2
//   //
//   calcNMin();
//   //
//   // If N(obs) is BELOW belt at s=0.0, then we have an empty belt.
//   //
//   if (getNObserved()<m_rejs0N1) {
//     m_lowerLimitFound = true;
//     m_lowerLimitNorm  = 1.0;
//     m_lowerLimit      = 0;
//     m_upperLimitFound = true;
//     m_upperLimitNorm  = 1.0;
//     m_upperLimit      = 0;
//     m_sumProb         = 0.0;
//     m_coversTruth     = false; // SHOULD THIS BE TRUE IF s=0??? (getTrueSignal()==0);
//     std::cout << "WARNING: Empty limit created!" << std::endl;
//     return true;
//   }
//   //
//   if (m_verbose>1) std::cout << "*** Scanning over all hypothesis" << std::endl;

//   //
//   // Find hypothesis start - s = (Nobs-bkg)/eff - should be in ~the middle of the belt
//   //
//   double mutest;
//   double mutestPrev;
//   int    i;
//   double p0;
//   bool done = false;

//   double mustart = getObservedSignal();
//   double muhigh  = mustart;
//   double mulow   = 0.0;
//   double dmu = 0.0;
//   int dir=100; // invalid number - check for bugs
//   //
//   if (m_verbose>1) {
//     std::cout << "*** Calculating belt for start s = " << mustart << std::endl;
//   }
//   calcLimit(mustart);
//   //
//   p0 = m_sumProb;
//   if (m_verbose>1) {
//     std::cout << "*** Obtained probability p = " << m_sumProb << std::endl;
//     std::cout << "*** N2(s=0)                = " << m_rejs0N2 << std::endl;
//   }
//   //
//   //------------------------------
//   // Start search for lower limit
//   //------------------------------
//   //
//   // Special case:
//   //  check if N(obs) < N2(s=0.0)
//   //  if so, the lower limit is 0.0
//   //
//   if (getNObserved()<m_rejs0N2) {
//     m_lowerLimitFound = true;
//     m_lowerLimitNorm  = 1.0;
//     m_lowerLimit      = 0;
//     if (m_verbose>1) {
//       std::cout << "*** Lower limit concluded from the fact that N(obs) = " << getNObserved()
//                 << " and Belt(s=0) = [ " << m_rejs0N1 << " : " << m_rejs0N2 << " ]" << std::endl;
//     }
//   } else {
//     //
//     // Scan for lower limit
//     //
//     done = false;
//     m_prevSumProb = p0;
//     i = 1;
//     //
//     // Special case:
//     //  N(obs) == N2(s=0)
//     //
//     if (getNObserved()==m_rejs0N2) {
//       if (m_verbose>1)
//         std::cout << "*** Special case: N(obs) = N2(s=0), test for limit at s=0" << std::endl;
//       dir = calcLimit(0);
//       if (m_verbose>1)
//         std::cout << "*** calcLimit at sL = 0 ==> P = " << m_sumProb
//                   << " and dir = " << dir
//                   << std::endl;
//       if (dir<1) {
//         done = true;
//         m_lowerLimitFound = true;
//         m_lowerLimitNorm  = m_scanBeltNorm;
//         m_lowerLimit = 0;
//         if (m_verbose>1)
//           std::cout << "*** Lower limit concluded to be sL = 0"
//                     << std::endl;
//       }
//     }
//     //
//     // Lower limit binary search
//     //
//     mutestPrev = mustart;
//     if ((m_verbose>1) && (!done))
//       std::cout << "*** Lower limit: binary search, start at mutest = " << mustart << std::endl;
//     while (!done) {
//       mutest = (muhigh+mulow)/2.0;
//       dmu = 2.0*fabs(mutestPrev-mutest)/(mutestPrev+mutest);
//       //
//       if (m_verbose>1)
//         std::cout << "*** Test at mutest  = " << mutest << std::endl;
//       //
//       dir = calcLimit(mutest);
//       if (m_verbose>1)
//         std::cout << "*** Got dir = " << dir << " and P(s) = " << m_sumProb << std::endl;
//       if ((dir==0) || (dmu<m_thresholdBS)) {
//         done = true;
//         m_lowerLimitFound = true;
//         m_lowerLimitNorm  = m_scanBeltNorm;
//         m_lowerLimit      = mutest;
//       } else if (dir==1) {
//         mulow = mutest;
//       } else {
//         muhigh = mutest;
//       }
//       mutestPrev = mutest;
//       i++;
//       if (i>1000) {
//         done=true;
//         std::cout << "*** WARNING calcLimit() failed  when scanning for LL! Infinite or long loop encountered - limit at 1000" << std::endl;
//       }
//     }
//     if ((m_lowerLimitFound) && (m_verbose>1)) {
//       std::cout << "*** Lower limit obtained from full scan; dir = " << dir << " dmu = " << dmu <<  " and sL = " << m_lowerLimit << std::endl;
//     }
//   }

//   if (m_lowerLimitFound) {
//     // Scan for upper limit
//     done = false;
//     m_prevSumProb = p0;
    
//     // find a rough range guaranteed to cover the upper limit
//     dir = -1;
//     mulow  = mustart;
//     muhigh = 2.0*(mustart+1.0)-m_lowerLimit; // pick an initial upper limit
//     // loop
//     if (m_verbose>1)
//       std::cout << "*** Find a rough upper limit, starting at " << muhigh << std::endl;
//     int cn=0;
//     double prec,prevPrec;
//     const double muscale=1.1;
//     double usescale = muscale;
//     double prevmu = muhigh;
//     prevPrec=0;
//     //
//     // Start at muhigh.
//     // calcLimit returns:
//     //---------------------------------------------
//     // dir | prec | action
//     //---------------------------------------------
//     //  -  | <0   | sumProb==0, muhigh at s_best
//     // +1  | >0   | sumProb>cl
//     // -1  | <0   | sumProb<cl
//     //  0  |  0   | sumProb = cl (within precision)
//     //---------------------------------------------
//     while (!done) { //dir==-1) {
//       dir = calcLimit(muhigh,prec);
//       if (prec<0) {
//         usescale = muscale;
//       } else {
//         if (cn>0) { // check direction
//           if (prec>prevPrec) { // if precision is degraded => we're going in the wrong direction
//             usescale = 1.0/usescale;
//           }
//         }
//         done = false;
//         // check the precision
//         if (cn>1) {
//           done = (prec>prevPrec);
//           if (done) muhigh = prevmu;
//         }
//         // step up counter of ok mu
//         cn++;
//         prevPrec = prec;
//         // check for infinite loop
//         if ((!done) && (cn>100)) {
//           done = true;
//           std::cout << "*** WARNING infinite (?) loop when scanning for the rough upper limit!" << std::endl;
//         }
//       }
//       // continue if not done
//       if (!done) {
//         prevmu = muhigh;
//         muhigh *= usescale;
//       }
//     }
//     //
//     // found a rough upper limit, now go down
//     //
//     if (m_verbose>1) {
//       std::cout << "*** Found a rough upper limit at " << muhigh << std::endl;
//       std::cout << "*** Upper limit: binary search" << std::endl;
//     }
//     i = 0;
//     done = false;
//     while (!done) {
//       mutest = (mulow+muhigh)/2.0;
//       if (m_verbose>1)
//         std::cout << "*** Test at mutest  = " << mutest << std::endl;
//       //      dmu = fabs(muhigh-mulow)/(muhigh+mulow);
//       dmu = 2.0*fabs(mutestPrev-mutest)/(mutestPrev+mutest);
//       dir = calcLimit(mutest,prec);
//       if (m_verbose>1)
//         std::cout << "*** Got dir = " << dir << " and P(s) = " << m_sumProb << std::endl;
//       if ((dir==0) || (dmu<m_thresholdBS)) { // see comment above for lower limit
//         done = true;
//         m_upperLimitFound = true;
//         m_upperLimitNorm  = m_scanBeltNorm;
//         m_upperLimit      = mutest;
//       } else if (dir==1) {
//         muhigh = mutest;
//       } else {
//         mulow = mutest;
//       }
//       mutestPrev = mutest;
//       prevPrec
//       i++;
//       if (i>1000) { // TODO: change this to a more smart check
//         done=true;
//         std::cout << "calcLimit() failed when scanning for UL! Infinite or long loop encountered - limit at 1000" << std::endl;
//       }
//     }
//     if ((m_upperLimitFound) && (m_verbose>1)) {
//       std::cout << "Upper limit obtained from full scan; dir = " << dir << " dmu = " << dmu <<  " and sU = " << m_upperLimit << std::endl;
//     }
//   }
//   //
//   bool rval=limitsOK();
//   if (rval) {
//     m_coversTruth = ((getTrueSignal()>=m_lowerLimit) && (getTrueSignal()<=m_upperLimit));
//     //    calcLimit(getTrueSignal());
// //     std::cout << "COV: ";
// //     TOOLS::coutFixed("s = ",4,getTrueSignal()); std::cout << "   ";
// //     TOOLS::coutFixed("p = ",4,m_sumProb); std::cout << "   ";
// //     TOOLS::coutFixed("l = ",4,m_lowerLimit); std::cout << "   ";
// //     TOOLS::coutFixed("u = ",4,m_upperLimit); std::cout << "   ";
// //     std::cout << "c = " << TOOLS::yesNo(m_coversTruth) << std::endl;
//   }
//   return limitsOK();
// }

bool Pole::calcCoverageLimit() {
  resetCalcLimit();
  calcNMin();
  if (getNObserved()<m_rejs0N1) {
    m_coversTruth = false;
  } else {
    calcLimit(getTrueSignal());
    m_coversTruth = (m_sumProb<m_cl); // s(true) inside belt - p(s)<cl
  }
  return true;
}

void Pole::resetCalcLimit() {
  m_maxNorm = -1.0;
  m_lowerLimitFound = false;
  m_upperLimitFound = false;
  m_lowerLimit = 0;
  m_upperLimit = 0;
  m_lowerLimitNorm = 0;
  m_upperLimitNorm = 0;
  m_lowerLimitPrec = -1;
  m_upperLimitPrec = -1;
  m_nBeltUsed = getNObserved();
  if (m_nBeltUsed<2) m_nBeltUsed=2;
  m_nBeltMinLast=0;
  m_nBeltMinUsed=m_nBeltUsed;
  m_nBeltMaxUsed=0;
  m_sumProb=0;
  m_scanBeltNorm=0;
  m_coversTruth=false;
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
    rval=calcCoverageLimit();
  } else {
    rval=calcLimit();
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
    cmtPre = "LIM:#";
    linePre = "LIM: ";
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
    std::cout << cmtPre << " Max N(belt) set  : " << m_nBeltUsed << std::endl;
    std::cout << cmtPre << " Max N(belt) used : " << m_nBeltMaxUsed << std::endl;
    std::cout << cmtPre << " Min N(belt) used : " << m_nBeltMinUsed << std::endl;
    std::cout << cmtPre << " Prob(belt s=0)   : " << m_rejs0P << std::endl;
    std::cout << cmtPre << " N1(belt s=0)     : " << m_rejs0N1 << std::endl;
    std::cout << cmtPre << " N2(belt s=0)     : " << m_rejs0N2 << std::endl;
    
    std::cout << cmtPre << "-----------------------------------------------------------------------------------------"
              << std::endl;
    std::cout << cmtPre << " N(obs)  eff(mean)     eff(sigma)      bkg(mean)       bkg(sigma)      low         up"
              << std::endl;
    std::cout << cmtPre << "-----------------------------------------------------------------------------------------"
              << std::endl;
  }
  std::cout << linePre;
  if (!simple) {
    TOOLS::coutFixed(6,getNObserved());   std::cout << "\t  ";
    TOOLS::coutFixed(6,getEffObs());      std::cout << "\t";
    TOOLS::coutFixed(6,getEffPdfSigma()); std::cout << "\t";
    TOOLS::coutFixed(6,getBkgObs());      std::cout << "\t";
    TOOLS::coutFixed(6,getBkgPdfSigma()); std::cout << "\t";
    TOOLS::coutFixed(6,m_scaleLimit*m_lowerLimit);     std::cout << "\t";
    TOOLS::coutFixed(6,m_scaleLimit*m_upperLimit);     std::cout << std::endl;
  }
  if (simple) {
    std::cout << "Limits = [ ";
    TOOLS::coutFixed(6,m_scaleLimit*m_lowerLimit); std::cout << ", ";
    TOOLS::coutFixed(6,m_scaleLimit*m_upperLimit); std::cout << " ]";
    std::cout << std::endl;
  }
}

void Pole::printSetup() {
  std::cout << "\n";
  std::cout << "================ P O L E ==================\n";
  std::cout << " 1.0 - conf. level  : " << 1.0-m_cl << std::endl;
  std::cout << " N observed         : " << getNObserved() << std::endl;
  std::cout << "----------------------------------------------\n";
  std::cout << " Coverage friendly  : " << TOOLS::yesNo(m_coverage) << std::endl;
  std::cout << " True signal        : " << getTrueSignal() << std::endl;
  std::cout << "----------------------------------------------\n";
  std::cout << " Efficiency meas    : " << getEffObs() << std::endl;
  std::cout << " Efficiency sigma   : " << getEffPdfSigma() << std::endl;
  std::cout << " Efficiency dist    : " << PDF::distTypeStr(getEffPdfDist()) << std::endl;
  std::cout << " Efficiency scale   : " << getEffScale() << std::endl;
  std::cout << "----------------------------------------------\n";
  std::cout << " Background meas    : " << getBkgObs() << std::endl;
  std::cout << " Background sigma   : " << getBkgPdfSigma() << std::endl;
  std::cout << " Background dist    : " << PDF::distTypeStr(getBkgPdfDist()) << std::endl;
  std::cout << " Background scale   : " << getBkgScale() << std::endl;
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
  std::cout << " Binary search thr. : " << m_thresholdBS << std::endl;
  std::cout << " 1-CL threshold     : " << m_thresholdPrec << std::endl;
  std::cout << " Min prob in belt   : " << m_minMuProb << std::endl;
  std::cout << "----------------------------------------------\n";
  std::cout << " *Test hyp. min     : " << m_hypTest.min() << std::endl;
  std::cout << " *Test hyp. max     : " << m_hypTest.max() << std::endl;
  std::cout << " *Test hyp. step    : " << m_hypTest.step() << std::endl;
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
  std::cout << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << "- ERROR: limit calculation failed. Possible causes:" << std::endl;
  std::cout << "- 1. binary search threshold too large" << std::endl;
  std::cout << "- 2. precision threshold too large" << std::endl;
  std::cout << "- 3. precision in integrations (eff,bkg) not sufficient" << std::endl;
  std::cout << "- 4. poisson table insufficient precision -> try without tabulated poisson" << std::endl;
  std::cout << "- 5. minimum probability in belt calculation too large; min = " << m_minMuProb << std::endl;
  std::cout << "-" << std::endl;
  std::cout << "- Results:" << std::endl;
  std::cout << "-    probability (should be ~ CL)  = " << getSumProb() << std::endl;
  std::cout << "-    lower lim norm (should be ~1) = " << getLowerLimitNorm() << std::endl;
  std::cout << "-    upper lim norm (ditto)        = " << getUpperLimitNorm() << std::endl;
  std::cout << "-    lower lim                     = " << getLowerLimit() << std::endl;
  std::cout << "-    upper lim                     = " << getUpperLimit() << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
}

