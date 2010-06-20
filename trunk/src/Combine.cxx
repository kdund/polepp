#include <iomanip>
#include <algorithm>
#include <iterator>
#include "Pole.h"
#include "Combine.h"


Combine::Combine() {
  m_foundLower=false;
  m_foundUpper=false;
  m_sVecMax = -1.0;
  m_sMaxUsed = -1.0;
  m_normProb = 1.0;
}

Combine::~Combine() {
}

bool Combine::isOK(const Pole *point) {
  //
  // Checks if point is OK:
  // * pointer not null
  // * object is properly initialized ( not done yet )
  // * nBelt>0 - no automatic assignment of nBelt allowed
  // * nBelt should be equal for all points
  // * CL should be the same, of course
  //
  bool rval=false;
  //
  // First make sure it is an ok object
  //
  if (point==0)            return rval;
  if (findPoint(point)>-1) return rval;
  //
  // Object is OK - make no further tests if it's the first object added
  //
  if (m_poleList.size()==0) return true;

  //
  // Compare this with the ref. (m_poleRef)
  //

  rval=true;
  return rval;
}

void Combine::add(const Pole *point) {
  if (isOK(point)) {
    if (m_poleList.size()==0) m_poleRef=point;
    m_poleList.push_back(point);
    m_nVecObs.push_back(point->getNObserved());
  }
}

void Combine::add(const std::vector<Pole *> poleList) {
  for (unsigned int i=0; i<poleList.size(); i++) {
    add(poleList[i]);
  }
}
//
// Correlation matrix:
//---------------------
//
//   b11 b12 ... b1N | c11 0    ...   |
//   b21 ...      .  | 0   c22  0 ... |
//    .           .  | .              |
//   bN1 ...     bNN | ...        cNN |
//   ----------------------------------
//   c11 0    ...    | e11 e12 ...e1N |
//   0   c22 0 ...   | e21 ...     .  |
//   .               |  .          .  |
//   0 ...       cNN | eN1 ...    eNN |
//
// where c11 is the correlation between b1 and e1

void Combine::initCorrelations() {
  unsigned int i,n;
  n = 2*m_poleList.size(); // n = 2 * n experiments
  m_corrMat.resize(n);
  m_corrMatFlag.resize(n);
  //
  // Set unit matrix
  //
  for (i=0; i<n; i++) {
    m_corrMat[i].resize(n,0.0);
    m_corrMat[i][i] = 1.0;
    m_corrMatFlag[i].resize(n,0);
    m_corrMatFlag[i][i] = 2;
  }
}
//
// Check is a pole is correlated with any other in a list
//
bool Combine::isCorrelated(const Pole *pole, std::vector<const Pole *> poleList, bool doEff) {
  if (pole==0) return false;
  int psize = int(poleList.size());
  int i=0;
  bool done=(i>=psize);
  //
  bool rval = false;
  int pindRef;
  int pind = findPoint(pole);
  int mofs = (doEff ? m_poleList.size():0);
  //
  if (pind>-1) {
    while (!done) {
      pindRef = findPoint(poleList[i]);
      if (pindRef>-1) {
	if (pind != pindRef) {
	  rval = ( ! isNotCorrelated(m_corrMat[pind+mofs][pindRef+mofs]) );
	}
      }
      i++;
      done = (rval || (i==psize));
    }
  }
  return rval;
}
//
// TODO: Add check that p1 and p2 are of proper distributions
// TODO: Hard-coded for 2 experiments
//
bool Combine::correlation(const Pole *p1, const Pole *p2, const double corr, bool isBkg) {
  int n1, n2, nexp;
  //
  bool rval = false;
  int nb1 = findPoint(p1);
  int nb2 = findPoint(p2);
  if ((nb1>-1) && (nb2>-1)) {
    nexp = m_poleList.size();
    int ne1 = nb1+nexp;
    int ne2 = nb2+nexp;
    if (isBkg) {
      n1 = nb1;
      n2 = nb2;
    } else {
      n1 = ne1;
      n2 = ne2;
    }
    int cflag = corrFlag(corr);
    m_corrMat[n1][n2] = corr;
    m_corrMat[n2][n1] = corr;
    m_corrMatFlag[n1][n2] = cflag;
    m_corrMatFlag[n2][n1] = cflag;
    //
    // Bkg-Eff correlation in experiment
    //
    double beCorr;
    beCorr = p1->getEffBkgCorr();
    cflag = corrFlag(beCorr);
    m_corrMat[ne1][nb1] = beCorr;
    m_corrMat[nb1][ne1] = beCorr;
    m_corrMatFlag[ne1][nb1] = cflag;
    m_corrMatFlag[nb1][ne1] = cflag;
    //
    // ditto for 2nd exp
    //
    beCorr = p2->getEffBkgCorr();
    cflag = corrFlag(beCorr);
    m_corrMat[ne2][nb2] = beCorr;
    m_corrMat[nb2][ne2] = beCorr;
    m_corrMatFlag[ne2][nb2] = cflag;
    m_corrMatFlag[nb2][ne2] = cflag;

    rval = true;
  }
  return rval;
}
//
// Creates the multidim. integral, taking into account the correlations between experiments
// TODO: hard-coded 2d
//
void Combine::makeCorrInt() {
  int nexp = m_poleList.size(); // n experiments
  int neint,nbint;
  const Range<double> *range;
  //
  std::vector<int> nBkgRange;
  std::vector<int> nEffRange;
  //
  // Loop over experiments and create ranges.
  // These are used to make vectors of all eff and bkg (indecis)
  // for the weight integral.
  //
  const Pole *polem1,*pole, *pole1, *pole2;

  std::vector<const Pole*> prevPoles;

  std::cout << "Looping over exps" << std::endl;
  polem1 = 0;
  for (int p=0; p<nexp; p++) {
    std::cout << " #" << p << std::endl;
    pole = m_poleList[p];
    if (!isBkgCorr(pole,prevPoles)) {
      range = pole->getBkgRangeInt();
      nbint = range->n();
      nBkgRange.push_back(nbint-1);
    }
    //
    if (!isEffCorr(pole,prevPoles)) {
      range = pole->getEffRangeInt();
      neint = range->n();
      std::cout << "Got eff range: " << neint << std::endl;
      nEffRange.push_back(neint-1);
    }
    //
    prevPoles.push_back(pole);
  }
  //
  int nbkg = makeIndVector(nBkgRange,m_bkgInt);
  std::cout << "Made bkg int vector = " << nbkg << std::endl;
  int neff = makeIndVector(nEffRange,m_effInt);
  std::cout << "Made eff int vector = " << neff << std::endl;
  std::cout << "Calculate weight" << std::endl;
  m_weightInt.resize(neff);
  std::cout << "Weight size = " << m_weightInt.size() << std::endl;
  // NOTE: if fully correlated this neff is change to the eff range of the first experiment
  //
  // Now, calculate the weight for each combination of bkg and eff.
  //
  ////  double de,eff;
  double pd;
  //
  double effPdf, bkgPdf;
  double de1,de2;
  double db1,db2;
  double bkg1,bkg2;
  double eff1,eff2;
  double bkgMeas1,bkgMeas2;
  double effMeas1,effMeas2;
  double effSigma1,effSigma2;
  double bkgSigma1,bkgSigma2;
  double edetc=0,esdetc=0,eseff1=0,eseff2=0,eveffc=0; // for 2d gauss
  double bdetc=0,bsdetc=0,bseff1=0,bseff2=0,bveffc=0; // for 2d gauss
  const Range<double> *erng1,*erng2;
  const Range<double> *brng1,*brng2;
  double ecorr,bcorr;
  bool effIsCorr, bkgIsCorr;
  double norm;
  //
  norm = 0.0;

  pole1 = m_poleList[0];
  pole2 = m_poleList[1];
  //
  // eff specifics
  //
  erng1 = m_poleList[0]->getEffRangeInt();
  erng2 = m_poleList[1]->getEffRangeInt();
  //
  ecorr = m_corrMat[2][3];
  std::cout << "Eff correlation = " << ecorr << std::endl;
  effIsCorr = isFullyCorrelated(ecorr);
  //
  if (effIsCorr) {
    neff = pole1->getEffRangeInt()->n();
    //    m_weightInt.resize(neff);
  } else {
    edetc = PDF::gGauss2D.getDetC(pole1->getEffSigma(),pole2->getEffSigma(),ecorr);
    esdetc = sqrt(edetc);
    eveffc = PDF::gGauss2D.getVeffCorr(edetc,pole1->getEffSigma(),pole2->getEffSigma(),ecorr);
    eseff1 = PDF::gGauss2D.getVeff(edetc,pole1->getEffSigma());
    eseff2 = PDF::gGauss2D.getVeff(edetc,pole2->getEffSigma());
    std::cout << edetc << " : " << edetc << " : " << eveffc << " : " << eseff1 << " : " << eseff2 << std::endl;
  }
  //
  effMeas1 = pole1->getEffMeas();
  effMeas2 = pole2->getEffMeas();

  effSigma1 = pole1->getEffSigma();
  effSigma2 = pole2->getEffSigma();
  //
  de1 = erng1->step();
  de2 = erng2->step();
  //
  // Bkg specifics
  //
  brng1 = m_poleList[0]->getBkgRangeInt();
  brng2 = m_poleList[1]->getBkgRangeInt();
  //
  bcorr = m_corrMat[0][1];
  std::cout << "Bkg correlation = " << bcorr << std::endl;
  bkgIsCorr = isFullyCorrelated(bcorr);

  if (bkgIsCorr) {
    nbkg = pole1->getBkgRangeInt()->n();
  } else {
    bdetc = PDF::gGauss2D.getDetC(pole1->getBkgSigma(),pole2->getBkgSigma(),bcorr);
    bsdetc = sqrt(bdetc);
    bveffc = PDF::gGauss2D.getVeffCorr(bdetc,pole1->getBkgSigma(),pole2->getBkgSigma(),bcorr);
    bseff1 = PDF::gGauss2D.getVeff(bdetc,pole1->getBkgSigma());
    bseff2 = PDF::gGauss2D.getVeff(bdetc,pole2->getBkgSigma());
    std::cout << bdetc << " : " << bdetc << " : " << bveffc << " : " << bseff1 << " : " << bseff2 << std::endl;
  }
    //
  bkgMeas1 = pole1->getBkgMeas();
  bkgMeas2 = pole2->getBkgMeas();

  bkgSigma1 = pole1->getBkgSigma();
  bkgSigma2 = pole2->getBkgSigma();
    //
  db1 = brng1->step();
  db2 = brng2->step();

  for (int i=0; i<neff; i++) {
    eff1 = erng1->getVal(m_effInt[i][0]); // eff val of b
    eff2 = erng2->getVal(m_effInt[i][1]); // ditto for exp 2
    if (effIsCorr) {
      effPdf = PDF::gGauss.getVal(eff1,effMeas1,effSigma1)*de1;
    } else {
      effPdf = PDF::gGauss2D.getVal2D(eff1,effMeas1,eff2,effMeas2,esdetc,eseff1,eseff2,eveffc)*de1*de2;
      effPdf = effPdf/((eff1+de1)*(eff2+de2)); //PRIOR
    }

    //    std::cout << "e1 = " << eff1 << "  e2 = " << eff2 << std::endl;

    ////////////////////// UNCORRELATED VERSION
    //
    // Calculate PDF: f1(e1)*f2(e2)*...fN(eN)*de1*de2*...deN
    // NOTE: Only uncorrelated gaussian is implemented at the moment
    //
//     for (int p=0; p<nexp; p++) {
//       pole = m_poleList[p];
//       //
//       eff = pole->getEffRangeInt()->getVal(m_effInt[i][p]); // get eff in point i and experiment p
//       de  = pole->getEffRangeInt()->step();
//       effPdf *= PDF::gGauss.getVal(eff,pole->getEffMeas(),pole->getEffSigma())*de; // multiply the efficiencies - no corr, gauss
//     }
    //
    // Background: gaussian, allowing correlation between experiments
    // NOTE: hardcoded at the moment for combining 2 experiments.
    //
    //    std::cout << " sdetc = " << sdetc << " veffc = " << veffc << " seff1 = " << seff1 << " seff2 = " << seff2
    //	      << " bmeas1 = " << bkgMeas1 << " bmeas2 = " << bkgMeas2 << std::endl;
    //
    for (int j=0; j<nbkg; j++) {
      bkg1 = brng1->getVal(m_bkgInt[j][0]); // bkg val of b
      bkg2 = brng2->getVal(m_bkgInt[j][1]); // ditto for exp 2
      if (bkgIsCorr) {
	bkgPdf = PDF::gGauss.getVal(bkg1,bkgMeas1,bkgSigma1)*db1; // fully correlated
      } else {
	bkgPdf  = PDF::gGauss2D.getVal2D(bkg1,bkgMeas1,bkg2,bkgMeas2,bsdetc,bseff1,bseff2,bveffc)*db1*db2;
	bkgPdf = bkgPdf/((bkg1+db1)*(bkg2+db2)); //PRIOR
      }
      //      std::cout << "b1 = " << bkg1 << "  b2 = " << bkg2 << std::endl;
      pd = bkgPdf*effPdf;
      m_weightInt[i].push_back(pd);
      norm += pd;
      //      std::cout << m_bkgInt[j][0] << ":" << m_bkgInt[j][1] << "   p(b)  = " << bkgPdf << " p(e)  = " << effPdf << " w = " << pd << std::endl;
    }
  }
  for (int i=0; i<neff; i++) {
    for (int j=0; j<nbkg; j++) {
      m_weightInt[i][j] = m_weightInt[i][j]/norm;
    }
  }
  std::cout << "Weight size = " << m_weightInt.size() << std::endl;
  std::cout << "Normalisation of weights: " << norm << std::endl;
}

// Same as previous but now:
// n = (e1+e2)*s + b
// e1 = uncorr eff
// e2 = corr eff
// Hack approach:
// * m_poleList[0,1] -> pole objects for the two experiments with the UNCORRELATED effs
// * m_poleCorr      -> a pole object with the CORRELATED effs, bkg and N can be arbitrary
//                      just needed in order to get the integration
//
void Combine::makeCorrIntCDF() {
  int nexp = m_poleList.size(); // n experiments
  int neint,nbint;
  const Range<double> *range;
  //
  std::vector<int> nBkgRange;
  std::vector<int> nEffRange;
  //
  // Loop over experiments and create ranges.
  // These are used to make vectors of all eff and bkg (indecis)
  // for the weight integral.
  //
  const Pole *polem1,*pole, *pole1, *pole2, *poleC;

  std::vector<const Pole*> prevPoles;

  std::cout << "Looping over exps" << std::endl;
  polem1 = 0;
  for (int p=0; p<nexp; p++) {
    std::cout << " #" << p << std::endl;
    pole = m_poleList[p];
    if (!isBkgCorr(pole,prevPoles)) {
      range = pole->getBkgRangeInt();
      nbint = range->n();
      nBkgRange.push_back(nbint-1);
    }
    //
    if (!isEffCorr(pole,prevPoles)) {
      range = pole->getEffRangeInt();
      neint = range->n();
      nEffRange.push_back(neint-1);
      std::cout << "Got eff range: " << neint << std::endl;
    }
    //
    prevPoles.push_back(pole);
  }
  //
  int nbkg = makeIndVector(nBkgRange,m_bkgInt);
  std::cout << "Made bkg int vector = " << nbkg << std::endl;
  int neff = makeIndVector(nEffRange,m_effInt);
  std::cout << "Made eff int vector = " << neff << std::endl;
  std::cout << "Calculate weight" << std::endl;
  //  m_weightInt.resize(neff);
  //  std::cout << "Weight size = " << m_weightInt.size() << std::endl;
  // NOTE: if fully correlated this neff is change to the eff range of the first experiment
  //
  // Now, calculate the weight for each combination of bkg and eff.
  //
  int neffC;
  ////  double de,eff;
  double pd;
  //
  double effPdf, bkgPdf, effPdfC;
  double de1,de2, deC;
  double db1,db2;
  double bkg1,bkg2;
  double eff1,eff2, effC;
  double bkgMeas1,bkgMeas2; //,bkgMeasC;
  double effMeas1,effMeas2,effMeasC;
  double effSigma1,effSigma2,effSigmaC;
  double bkgSigma1,bkgSigma2;
  const Range<double> *erng1,*erng2,*erngC;
  const Range<double> *brng1,*brng2;
  //  double ecorr,bcorr;
  //  bool effIsCorr, bkgIsCorr;
  double norm;
  //
  pole1 = m_poleList[0];
  pole2 = m_poleList[1];
  poleC = m_poleCorr;
  //
  // eff specifics
  //
  erng1 = m_poleList[0]->getEffRangeInt();
  erng2 = m_poleList[1]->getEffRangeInt();
  erngC = m_poleCorr->getEffRangeInt();
  //
  neffC = erngC->n();
  if (neffC!=erng1->n()) std::cout << "WARNING::Correlated and uncorrelated eff range not equal! The end is near..." << std::endl;

  effMeas1 = pole1->getEffMeas();
  effMeas2 = pole2->getEffMeas();
  effMeasC = poleC->getEffMeas();

  effSigma1 = pole1->getEffSigma();
  effSigma2 = pole2->getEffSigma();
  effSigmaC = poleC->getEffSigma();
  //
  de1 = erng1->step();
  de2 = erng2->step();
  deC = erngC->step();
  //
  // Bkg specifics
  //
  brng1 = m_poleList[0]->getBkgRangeInt();
  brng2 = m_poleList[1]->getBkgRangeInt();
  //
  nbkg = brng1->n();

  bkgMeas1 = pole1->getBkgMeas();
  bkgMeas2 = pole2->getBkgMeas();

  bkgSigma1 = pole1->getBkgSigma();
  bkgSigma2 = pole2->getBkgSigma();
    //
  db1 = brng1->step();
  db2 = brng2->step();

  m_effIntC.resize(neffC,0);
  m_wbkg.resize(nbkg,1.0);
  m_weff.resize(neff,1.0);
  m_weffC.resize(neffC,1.0);

  //
  // TMP CDF solution for combining 2 results
  norm = 0.0;
  for (int j=0; j<nbkg; j++) {
    bkg1 = brng1->getVal(m_bkgInt[j][0]); // bkg val of b
    bkg2 = brng2->getVal(m_bkgInt[j][1]); // ditto for exp 2
    bkgPdf = PDF::gGauss.getVal(bkg1,bkgMeas1,bkgSigma1)*
      PDF::gGauss.getVal(bkg2,bkgMeas2,bkgSigma2)*db1*db2; // uncorrelated
    m_wbkg[j] = bkgPdf;
    for (int i=0; i<neff; i++) {
      eff1 = erng1->getVal(m_effInt[i][0]); // eff val of b
      eff2 = erng2->getVal(m_effInt[i][1]); // ditto for exp 2
      effPdf = PDF::gGauss.getVal(eff1,effMeas1,effSigma1)*
	PDF::gGauss.getVal(eff2,effMeas2,effSigma2)*de1*de2;
      //      std::cout << "b1 = " << bkg1 << "  b2 = " << bkg2 << std::endl;
      m_weff[i]=effPdf;
      for (int c=0; c<neffC; c++) {
	effC = erngC->getVal(c);
	m_effIntC[c] = effC;
	effPdfC = PDF::gGauss.getVal(effC,effMeasC,effSigmaC)*deC;
	m_weffC[c]=effPdfC;
	//
	pd = bkgPdf*effPdf*effPdfC;
	//	m_weightInt[i].push_back(pd);
	norm += pd;
      }
    }
  }
  //
  // "Normalise" one weight array since it's the full product which is truly normalized
  //
//   for (int j=0; j<nbkg; j++) {
//     m_wbkg[j] = m_wbkg[j]/norm;
//   }

  std::cout << "Weight size = " << m_weightInt.size() << std::endl;
  std::cout << "Normalisation of weights: " << norm << std::endl;
}

const double Combine::calcProb(std::vector<int> nvec, double s) const {  
  double prob;
  double sumProb = 0.0;
  double g;
  //
  // loop over all eff and bkg vectors
  // loop over all experiments
  // p = Sum{e,b,p} w[e][b]*Po(n[p]|e*s+b)
  //
  const Range<double> *effr,*bkgr;
  for(unsigned int i=0; i<m_effInt.size(); i++) {
    for(unsigned int j=0; j<m_bkgInt.size(); j++) {
      prob = 1.0;
      if (m_weightInt[i][j]>1e-16) {
	for (unsigned int p=0; p<m_poleList.size(); p++) {
	  effr = m_poleList[p]->getEffRangeInt();
	  bkgr = m_poleList[p]->getBkgRangeInt();
	  g = effr->getVal(m_effInt[i][p])*s + bkgr->getVal(m_bkgInt[j][p]);
	  prob *= PDF::gPoisson.getVal(nvec[p],g);
	}
	prob *= m_weightInt[i][j];
	sumProb += prob;
      }
    }
  }
  //  if (m_verbose>9) std::cout << "calcProb() : " << n[0] << ", " << s << " => p = " << p << std::endl;
  return sumProb;
}

const double Combine::setNormProbCDF() {
  double rval=0.0;
  m_normProb = 1.0;
  for (unsigned int n=0; n<m_nVectors.size(); n++) {
    rval += calcProbCDF(m_nVectors[n],0.0);
  }
  m_normProb = rval;
  return rval;
}

// TEMP SOLUTION FOR COMBINING CDF : Must call setNormProbCDF() before
const double Combine::calcProbCDF(std::vector<int> nvec, double s) const {  
  double prob;
  double sumProb = 0.0;
  double g;
  //
  // loop over all eff and bkg vectors
  // loop over all experiments
  // p = Sum{e,b,p} w[e][b]*Po(n[p]|e*s+b)
  //
  const Range<double> *effr,*bkgr;
  double w;
  for(unsigned int j=0; j<m_bkgInt.size(); j++) {
    for(unsigned int i=0; i<m_effInt.size(); i++) {
      for(unsigned int c=0; c<m_effIntC.size(); c++) {
	prob = 1.0;
	w = m_wbkg[j]*m_weff[i]*m_weffC[c];
	//	std::cout << "w,b,e,c: " << w << " --- " << m_wbkg[j] << ", " << m_weff[i] << ", " << m_weffC[c] << std::endl;
	if (w>1e-16) {
	  for (unsigned int p=0; p<m_poleList.size(); p++) {
	    effr = m_poleList[p]->getEffRangeInt();
	    bkgr = m_poleList[p]->getBkgRangeInt();
	    g = (m_effIntC[c] * effr->getVal(m_effInt[i][p]))*s + bkgr->getVal(m_bkgInt[j][p]);
	    prob *= PDF::gPoisson.getVal(nvec[p],g);
	  }
	  prob *= w;
	  sumProb += prob;
	}
      }
    }
  }
  //  if (m_verbose>9) std::cout << "calcProb() : " << n[0] << ", " << s << " => p = " << p << std::endl;
  return sumProb/m_normProb;
}

int Combine::makeIndVector(const int ndim, std::vector< std::vector<int> > & nvec) {
  //
  int nmeas = m_poleList.size();
  nvec.clear();
  //
  std::vector< int > jj;
  jj.resize(nmeas,0);
  nvec.push_back(jj);  
  while (Combination::next_vector(jj,ndim)) {
    nvec.push_back(jj);
  }
  //
  // sort them
  //
  std::sort(nvec.begin(),nvec.end());

  return nvec.size();
}

int Combine::makeIndVector(const std::vector<int> & ndim, std::vector< std::vector<int> > & nvec) {
  //
  int nmeas = m_poleList.size();
  nvec.clear();
  //
  std::vector< int > jj;
  jj.resize(nmeas,0);
  nvec.push_back(jj);
  while (Combination::next_vector(jj,ndim)) {
    nvec.push_back(jj);
  }
  //
  // sort them
  //
  std::sort(nvec.begin(),nvec.end());

  return nvec.size();
}

void Combine::tabulateLikelihoodCorr() {
  //
  // * Fill m_nVectors
  // * Calculate a proper range of s
  // * Calculate L(n,s) and save them
  //
  std::cout << "Tabulating likelihood function (CORR)" << std::endl;
  std::cout << "- Making nVector...." << std::endl;
  int nmeas = m_poleList.size();
  const Pole *pole;
  //
  std::vector<int> nrange;
  int nBeltMax=-1;
  int nbelt;
  int nvSize=1;
  for (int i=0; i<nmeas; i++) {
    pole = m_poleList[i];
    nbelt = pole->getNBelt();
    nrange.push_back(nbelt);
    if (nbelt>nBeltMax) nBeltMax = nbelt;
    nvSize *=(nbelt+1);
    std::cout << "  " << i+1 << ". N(belt) = " << nbelt << std::endl;
  }
  //
  makeIndVector(nrange,m_nVectors);
  //
  //
  // -------------  nVectors done!
  //
  //
  // Obtain a proper range of s (m_sVector)
  // Use the maximum n and find the estimated s(up) for all points in poleList.
  // OLD: Scale the maximum s with some safety-factor
  // NEW: Use the minimum s - the combined limit should be less than that.
  //
  std::cout << "- Obtaining range of s...." << std::endl;
  double smax=0;
  //  double smin=1000000.0;
  double s=0;
  double h;
  double hmax=0;
  double hmin=1000000.0;
  for (int i=0; i<nmeas; i++) {
    pole = m_poleList[i];
    //
    // Find hypthesis scan range
    //
    h = pole->getHypTest()->max();
    if (h>hmax) hmax=h;
    if (h<hmin) hmin=h;
    //
    // Find max s for the tabulation of L(n,s)
    //
    s = (static_cast<double>(pole->getNBelt()) - pole->getBkgIntMin())/pole->getEffIntMax();
    if (s>smax) smax = s;
  }
  if (m_sVecMax>0.0) smax = m_sVecMax;
  m_sMaxUsed = -1.0;
  //
  // Get hypthesis step size from ref.
  //
  double hstep = m_poleRef->getHypTest()->step();
  if (hstep<=0) hstep = 0.01;
  //  sstep=0.1;
  //
  // Fill sVector and sHypRange with values
  //
  m_sVector.clear();
  m_sHypRange.setRange(0.0,hmin,hstep); // make range
  //
  double smaxStep = 0.1;
  int svSize = static_cast<int>((smax + smaxStep)/smaxStep + 0.5);
  for (int i=0; i<svSize; i++) {
    m_sVector.push_back(static_cast<double>(i)*smaxStep);
  }
  //
  std::cout << "- Hypothesis scan range     : [0.0, " << hmin << "]" << std::endl;
  std::cout << "- s(max) in tabulated L(n,s): [0.0, " << smax << "]" << std::endl;
  //
  //
  // -------- sVector done
  //
  //
  // Loop over all (nVectors, sVector) and calulate L(n,s)
  //
  // Save only those L(n,s)>phLim (=0.0001) in order to reduce size of matrix.
  //
//   std::cout << "Size of nVectors   : " << m_nVectors.size() << std::endl;
//   std::cout << "Size of nVectors[0]: " << m_nVectors[0].size() << std::endl;
//   std::cout << "Size of sVector    : " << m_sVector.size() << std::endl;
//   std::cout << "Size of poleList   : " << m_poleList.size() << std::endl;
//
  std::cout << "--- Calculate L(n,s) ---" << std::endl;
  std::cout << "- Max number of L(n,s) to be calculated = " << m_nVectors.size()*m_sVector.size() << std::endl;
  std::cout << "- Normalizing PDF ...." << std::endl;
  setNormProbCDF();
  //
  m_likeliHood.resize(m_nVectors.size());
  m_sMinInd.resize(m_nVectors.size());
  m_sMaxInd.resize(m_nVectors.size());
  int nok=0;
  int ntot=0;
  double plh,lhmax,lhsmax;
  int ismin,ismax; // indecis of min and max s for the current n
  bool isminOK;
  const double plhLim = 0.0000001; // min Lh
  std::vector<double> tmpLh;  // temporary storage
  tmpLh.resize(m_sVector.size(),0.0);
  double sumLhTot=0.0;
  double sumLhS=0.0;
  std::vector<double> sumLhN;
  sumLhN.resize(m_sVector.size(),0.0);
  //
  for (unsigned int n=0; n<m_nVectors.size(); n++) {
    lhmax = -1.0;
    lhsmax = -1.0;
    if (n%10==0)
      std::cout << "  N vector " << n << " out of " << m_nVectors.size() << std::endl;
    m_sMinInd[n] = -1;
    m_sMaxInd[n] = -1;
    isminOK = false;
    ismin = -1;
    ismax = -1;
    sumLhS=0.0;
    for (unsigned int m=0; m<m_sVector.size(); m++) {
      s = m_sVector[m];
      // NOTE: DO NOT FORGET TO NORMALIZE PDF BEFORE...
      plh = calcProbCDF(m_nVectors[n],s);     // L(nvec,s(m)) // TEMP CDF SOLUTION
      tmpLh[m] = plh; // save it temporarily
      sumLhS += plh;
      sumLhTot += plh;
      sumLhN[m] += plh;
      //
      // lh might fluctuate around min before steadily remaining above threshhold.
      // The same is true when the upper s is reached.
      //
      if ((plh>plhLim) && (!isminOK)) { // first time above threshhold - save s
	ismin = m;
	isminOK = true;
      }
      if (plh>plhLim) {
	ismax = m;
	//	std::cout << "s(m) = " << s << " => " << plh << std::endl;
      }
      if (plh>lhmax) {
	lhmax = plh;
	lhsmax = s;
      }
      ntot++;
    }
    //    std::cout << "Sum Lh (n=" << n << ",s=" << s << ") = " << sumLhS << std::endl;
    for (int m=ismin; m<ismax; m++) {
      m_likeliHood[n].push_back(tmpLh[m]/sumLhS);
      nok++;
    }
    //    std::cout << "N: " << n << "   s range = " << ismin << ":" << ismax << std::endl;
    //    std::cout << "       lhmax   = " << lhmax << " at s = " << lhsmax << std::endl;
    m_sMinInd[n] = ismin;
    m_sMaxInd[n] = ismax;
  }
  for (unsigned int k=0; k<sumLhN.size(); k++) {
    std::cout << "Sum Lh over N : " << k << " , " << sumLhN[k] << std::endl;
  }
  std::cout << "- Actual number of L(n,s) calculated = " << nok << std::endl;
  std::cout << "- Sum of L(n,s) = " << sumLhTot << std::endl;
  std::cout << "- Tabulating DONE!" << std::endl;
  //
  // -------- Likelihood construction done
  //
}


void Combine::tabulateLikelihood() {
  //
  // * Fill m_nVectors
  // * Calculate a proper range of s
  // * Calculate L(n,s) and save them
  //
  std::cout << "Tabulating likelihood function" << std::endl;
  std::cout << "- Making nVector...." << std::endl;
  int nmeas = m_poleList.size();
  const Pole *pole;
  //
  std::vector<int> nrange;
  int nBeltMax=-1;
  int nbelt;
  int nvSize=1;
  for (int i=0; i<nmeas; i++) {
    pole = m_poleList[i];
    nbelt = pole->getNBelt();
    nrange.push_back(nbelt);
    if (nbelt>nBeltMax) nBeltMax = nbelt;
    nvSize *=(nbelt+1);
    std::cout << "  " << i+1 << ". N(belt) = " << nbelt << std::endl;
  }
  //
  makeIndVector(nrange,m_nVectors);
  //
//   unsigned int ss;
//   for (unsigned int n=0; n<m_nVectors.size(); n++) {
//     std::cout << n;
//     ss = m_nVectors[n].size();
//     for (unsigned int i=0; i<ss; i++) {
//       std::cout << "\t" << m_nVectors[n][i];
//     }
//     std::cout << std::endl;
//   }
  //
  // -------------  nVectors done!
  //
  //
  // Obtain a proper range of s (m_sVector)
  // Use the maximum n and find the estimated s(up) for all points in poleList.
  // OLD: Scale the maximum s with some safety-factor
  // NEW: Use the minimum s - the combined limit should be less than that.
  //
  std::cout << "- Obtaining range of s...." << std::endl;
  double smax=0;
  //  double smin=1000000.0;
  double s;
  double h;
  double hmax=0;
  double hmin=1000000.0;
  for (int i=0; i<nmeas; i++) {
    pole = m_poleList[i];
    //
    // Find hypthesis scan range
    //
    h = pole->getHypTest()->max();
    if (h>hmax) hmax=h;
    if (h<hmin) hmin=h;
    //
    // Find max s for the tabulation of L(n,s)
    //
    s = (static_cast<double>(pole->getNBelt()) - pole->getBkgIntMin())/pole->getEffIntMax();
    if (s>smax) smax = s;
  }
  if (m_sVecMax>0.0) smax = m_sVecMax;
  m_sMaxUsed = -1.0;
  //
  // Get hypthesis step size from ref.
  //
  double hstep = m_poleRef->getHypTest()->step();
  if (hstep<=0) hstep = 0.01;
  //  sstep=0.1;
  //
  // Fill sVector and sHypRange with values
  //
  m_sVector.clear();
  m_sHypRange.setRange(0.0,hmin,hstep); // make range
  //
  int svSize = static_cast<int>((smax + hstep)/hstep + 0.5);
  for (int i=0; i<svSize; i++) {
    m_sVector.push_back(static_cast<double>(i)*hstep);
  }
  //
  std::cout << "- Hypothesis scan range     : [0.0, " << hmin << "]" << std::endl;
  std::cout << "- s(max) in tabulated L(n,s): [0.0, " << smax << "]" << std::endl;
  //
  //
  // -------- sVector done
  //
  //
  // Loop over all (nVectors, sVector) and calulate L(n,s)
  //
  // Save only those L(n,s)>phLim (=0.0001) in order to reduce size of matrix.
  //
//   std::cout << "Size of nVectors   : " << m_nVectors.size() << std::endl;
//   std::cout << "Size of nVectors[0]: " << m_nVectors[0].size() << std::endl;
//   std::cout << "Size of sVector    : " << m_sVector.size() << std::endl;
//   std::cout << "Size of poleList   : " << m_poleList.size() << std::endl;
//
  std::cout << "- Calculate L(n,s)...." << std::endl;
  m_likeliHood.resize(m_nVectors.size());
  std::cout << "- Max number of L(n,s) to be calculated = " << m_nVectors.size()*m_sVector.size() << std::endl;
  m_sMinInd.resize(m_nVectors.size());
  m_sMaxInd.resize(m_nVectors.size());
  m_lhNorm.resize(m_nVectors.size());
  int nok=0;
  int ntot=0;
  double plh,lh;
  int ismin,ismax; // indecis of min and max s for the current n
  bool isminOK;
  const double plhLim = 0.0001; // min Lh
  std::vector<double> tmpLh;  // temporary storage
  tmpLh.resize(m_sVector.size(),0.0);
  //
  double sumTot=0;
  for (unsigned int n=0; n<m_nVectors.size(); n++) {
    if (n%10==0) std::cout << "  N vector " << n << " out of " << m_nVectors.size() << std::endl;
    m_sMinInd[n] = -1;
    m_sMaxInd[n] = -1;
    m_lhNorm[n] = 1.0;
    isminOK = false;
    ismin = -1;
    ismax = -1;

    for (unsigned int m=0; m<m_sVector.size(); m++) {
      s = m_sVector[m];
      plh = 1.0;
      for (unsigned int p=0; p<m_poleList.size(); p++) { // loop over all points
	pole = m_poleList[p];
	lh = pole->calcProb(m_nVectors[n][p],s);     // L(n(p),s(m))
	plh*=lh;                     // L(n1,s(m))*L(n2,s(m))*..../
	//	if (n==12) std::cout << "n,s,lh = " << p << ": " << m_nVectors[n][p] << ", " << s << ", " << lh << std::endl;
      }
      tmpLh[m] = plh; // save it temporarily
      //
      // lh might fluctuate around min before steadily remaining above threshhold.
      // The same is true when the upper s is reached.
      //
      if ((plh>plhLim) && (!isminOK)) { // first time above threshhold - save s
	ismin = m;
	isminOK = true;
      }
      if (plh>plhLim) ismax = m;
      ntot++;
    }
    double tmpSum=0;
    for (int m=ismin; m<ismax; m++) {
      m_likeliHood[n].push_back(tmpLh[m]);
      tmpSum += tmpLh[m];
      nok++;
    }
    sumTot += tmpSum;
    m_sMinInd[n] = ismin;
    m_sMaxInd[n] = ismax;
    m_lhNorm[n] = tmpSum; // THIS IS NOT A NORM!
  }
  std::cout << "NORM: " << sumTot << std::endl;
  std::cout << "- Actual number of L(n,s) calculated = " << nok << std::endl;
  std::cout << "- Tabulating DONE!" << std::endl;
  //
  // -------- Likelihood construction done
  //
}


unsigned int Combine::getIndVector(std::vector<int> & nvec, std::vector< std::vector<int> > & allVectors) const {
  std::vector<int> nmax =  allVectors.back();     // vector containg the maximum n(belt) per Pole point
  unsigned int nsize = allVectors.front().size(); // number of points
  unsigned int nFullSize = allVectors.size();
  if (nvec.size() != nsize) {
    std::cerr << "WARNING: Input vector of faulty size (" << nvec.size() << " <> " << nsize << ")\n";
    return 0;
  }
  //
  // Find index in array of vectors matching the given vector
  //
  unsigned int ofs = 0;
  unsigned int npp=0;
  int nvecInd;
  for (unsigned int i=0; i<nsize; i++) {
    nvecInd = nsize-i-1;
    if (i==0) {
      npp=1;
    } else {
      npp *= nmax[nvecInd+1]+1;
    }
    if (nvec[nvecInd]>int(nmax[nvecInd])) {
      std::cout << "WARNING: Too large element in input vector (" << nvec[nvecInd] << " > " << nmax[nvecInd] << ")\n";
      std::cout << "         " << nvecInd << ",  nsize = " << nsize << ",  nmax = " << nmax[i] << std::endl;
      return 0;
    }
    ofs += npp*(nvec[nvecInd]);
  }
  if (ofs>=nFullSize) {
    std::cout << "ERROR: Offset too large - input vector not in range or...bug? Ofs = " << ofs << std::endl;
    return 0;
  }
  return ofs;
}

double Combine::getLikelihood(int nind, double s) {
  //  bool verb = (nind==12);
  if ((nind<0) || (nind>int(m_nVectors.size()))) return 0.0;
  int ismin = m_sMinInd[nind];
  int ismax = m_sMaxInd[nind];
  if (ismin==ismax) return 0.0;
  double smin = m_sVector[ismin];
  double smax = m_sVector[ismax];
  //  if (verb) {
  //  std::cout << "getLH -> s, s range: " << s << " => " << smin << "," << smax << std::endl;
  //  }
  if (s<smin) return 0.0;
  if (s>smax) return 0.0;
  //
  //
  // Establish closest s in m_sVector
  //
  int ind0, ind1;
  double step=0.01;
  if (m_sVector.size()>1)
    step= m_sVector[1]-m_sVector[0];
  int index = int(s/step);
  double sA = m_sVector[index];
  double s0, s1, ds;
  int ssize = ismax-ismin; // nr of saved lhR
  if (sA>s) {
    ind1 = index;
    ind0 = (index>1 ? index-1:0);
    s0 = m_sVector[ind0];
    s1 = sA;
  } else {
    ind1 = (index+1<ismax ? index+1:ismax-1);
    ind0 = index;
    s0 = sA;
    s1 = m_sVector[ind1];
  }
  //  std::cout << "getLh: ind0/1 = " << ind0 << ":" << ind1 << std::endl;
  ds = s1-s0;
  ind0 = ind0 - ismin;
  ind1 = ind1 - ismin;
  if (ind1>=ssize) ind1 = ssize-1;
  if (ind0>=ssize) ind0 = ssize-1;
  if (ind0<0)      ind0 = 0;
  double lh0 = m_likeliHood[nind][ind0];
  double lh1 = m_likeliHood[nind][ind1];
  double rval;
  //
  // linear interpolation
  //
  if (ind0==ind1) {
    rval = lh0;
  } else {
    rval = ((lh1-lh0)/(s1-s0))*(s-s0) + lh0;
  }
 //  if (rval>2.0) {
//     std::cout << "N(ind) = " << nind << std::endl;
//     std::cout << "ind0   = " << ind0 << std::endl;
//     std::cout << "ind1   = " << ind1 << std::endl;
//     std::cout << "s0     = " << s0 << std::endl;
//     std::cout << "s1     = " << s1 << std::endl;
//     std::cout << "s      = " << s << std::endl;
//     std::cout << "lh0    = " << lh0 << std::endl;
//     std::cout << "lh1    = " << lh1 << std::endl;
//     std::cout << "lh     = " << rval << std::endl;
//   }
  return rval;
}

void Combine::setBestMuScan(int nind) {
  //
  // loops over all Pole points and finds the minimum range required
  //
  const int nmeas = m_poleList.size();
  const Pole *pole;
  double smin=100000;
  double smax=0;
  double s;
  int n;
  double emin,emax,bmin,bmax,bmeas,emeas;
  //
  for (int i=0; i<nmeas; i++) {
    pole = m_poleList[i];
    n = m_nVectors[nind][i];

    bmeas = pole->getBkgMeas();
    bmin  = pole->getBkgRangeInt()->min();
    bmax  = pole->getBkgRangeInt()->max();

    emin  = pole->getEffRangeInt()->min();
    emax  = pole->getEffRangeInt()->max();
    emeas = pole->getEffMeas();
    s = (static_cast<double>(n) - bmin)/emeas;
    if (s>smax) smax=s;
    //
    s = (static_cast<double>(n) - bmax)/emax;
    if (s<0.0) s=0.0;
    if (s<smin) smin=s;
  }
  m_bestMuScan.setRange(smin,smax,m_poleRef->getDmus());
  //  std::cout << "Set s_best scan range to: [ " << smin << ":" << smax << " ] with ds = " << m_poleRef->getDmus() << std::endl;
}

void Combine::findBestMu(int ind) {
  //
  // Finds s(best) for the set of N(obs) given by m_nVectors.
  // The index ind is used for chosing the element to fill
  // in the m_bestMu/Prob vectors
  //
  if (m_poleList.size()==0) {
    std::cerr << "WARNING: No measurements added" << std::endl;
    return;
  }
  if (m_poleList.size()!=m_nVectors[0].size()) {
    std::cerr << "WARNING: Size of input vector does not match the #of measurements!" << std::endl;
    std::cerr << "nvec = " << m_nVectors[0].size() << std::endl;
    std::cerr << "pole = " << m_poleList.size() << std::endl;
    return;
  }
  //
  setBestMuScan(ind); // set bestMu scan range
  //
  // Find the s that maximizes the likelihood
  //
  double lhMax=-1000.0;
  double sbest = 0;
  double sscan;
  double lh;
  for (int j=0; j<m_bestMuScan.n(); j++) {
    sscan = m_bestMuScan.getVal(j);
    lh = getLikelihood(ind,sscan);
    if (lh>lhMax) {
      lhMax = lh;
      sbest = sscan;
    }
  }
  //   std::cout << "SBEST: ";
  //   for (unsigned int i=0; i<m_poleList.size(); i++) {
  //     std::cout << m_nVectors[ind][i] << " ";
  //   }
  //   std::cout << sbest << " " << std::endl;
  m_bestMu[ind]     = sbest;
  m_bestMuProb[ind] = lhMax;
  //
  //  if (sbest>m_sMaxUsed) m_sMaxUsed=sbest;
}

void Combine::findBestMu() {
  // loop over ALL vectors of n....
  std::cout << "- Finding all s_best..." << std::endl;

  //
  // Loop over all vectors
  //
  for (unsigned int i=0; i<m_nVectors.size(); i++) {
//     std::cout << i << "  -  ";
//     for (unsigned int m=0; m<m_poleList.size(); m++) {
//       std::cout << (m_nVectors[i][m]) << " ";
//     }
//     std::cout << std::endl;
    findBestMu(i);
    //    std::cout << "s(best) = " << m_bestMu[i] << " , P = " << m_bestMuProb[i] << std::endl;
  }
  std::cout << "- Done!" << std::endl;
}

double Combine::calcLimit(double s) {
  double normProb=0;
  double p;
  double pmax,pmin; //,pave;
  int natmax,natmin;
  pmax=-1000;
  pmin=1.0e28;
  natmin=-1;
  natmax=-1;
  for (unsigned n=0; n<m_nVectors.size(); n++) {
    p = getLikelihood(n,s);
    if (m_bestMuProb[n]>0.0) {
      m_lhRatio[n] = p/m_bestMuProb[n];
    } else {
      m_lhRatio[n] = 0.0; //(n==m_indexNobs ? 0.0:0.001); // this is just to make sure that CHECK!!!!
    }
    //    std::cout << "lhRatio[" << n << "] = " << m_lhRatio[n] << " , " << p << std::endl;
    normProb += p;
    if ((p<pmin)||(natmin<0)) {
      pmin=p;
      natmin=n;
    }
    if ((p>pmax)) {
      pmax=p;
      natmax=n;
    }
  }
//   std::cout << "NORMPROB min at " << natmin << " = " << pmin << " , RL = " << m_lhRatio[natmin] << std::endl;
//   std::cout << "NORMPROB max at " << natmax << " = " << pmax << " , RL = " << m_lhRatio[natmax] << std::endl;
//   std::cout << "NORMPROB at s= " << s << " => " << normProb << std::endl;
//   std::cout << "LHRATIO at nind " << m_indexNobs << " = " << m_lhRatio[m_indexNobs] << std::endl; 
 //
  bool done=false;
  unsigned int n=0;
  double sumProb=0;
  double lh=0;
//   if (normProb>0.9) {
//     std::cout << "Index N(obs) = " << m_indexNobs << std::endl;
//   }
//  std::cout << "\nSearch..." << std::endl;
  while (!done) {
    if (n!=m_indexNobs) {
      if (m_lhRatio[n]>m_lhRatio[m_indexNobs]) {
	lh = getLikelihood(n,s);
	sumProb += lh;
      }
    }
    //    std::cout << " n = " << n << " => lh = " << lh << "  RL: " << m_lhRatio[n] << " > " << m_lhRatio[m_indexNobs] << std::endl;
    n++;
    done = ((n==m_lhRatio.size()) || (sumProb>m_poleRef->getCL()));
  }
//   if (normProb>0.9) {
//     std::cout << "Sumprob = " << sumProb << std::endl;
//   }
  //
  double dCL = sumProb-m_poleRef->getCL();
  //  std::cout << "calcLimits(): " << s << " norm = " << normProb << " sump = " << sumProb << std::endl;
  if (dCL<0.0) {
    if (m_foundLower) {
      m_upperLimit = s;
      m_uppNorm = normProb;
      m_uppProb = sumProb;
    } else {
      m_lowerLimit = s;
      m_lowNorm = normProb;
      m_lowProb = sumProb;
      m_foundLower = true;
      m_foundUpper = false;
    }
  } else {
    if (m_foundLower) {
      m_foundUpper = true;
    }
  }
  if (m_foundLower) {
    //    std::cout << s << " : sumProb = " << sumProb << " : normProb = " << normProb << std::endl;
  }
  return normProb;
}

void Combine::calcLimits() {
  std::cout << "- Finding limits..." << std::endl;
  m_foundLower = false;
  m_foundUpper = false;
  m_lowerLimit = 0;
  m_upperLimit = 0;
  m_lowProb=0;
  m_lowNorm=0;
  m_uppProb=0;
  m_uppNorm=0;
  //
  bool done=false;
  int ihmax =  m_sHypRange.n();
  int i=0;
  double htst;
  double pnorm;
  //
  double pnormMax=0;
  int indUpp=0;
  while (!done) {
    htst = m_sHypRange.getVal(i);
    pnorm = calcLimit(htst);
    if (pnorm>pnormMax) pnormMax = pnorm;
    i++;
    done = ((i==ihmax) ||
	    (m_foundUpper));// ||
    //	    (!(fabs(pnorm-1.0)<1.0e-1)));
    //	    (!m_poleRef->normOK(pnorm)));
    if (done) indUpp=i; 
  }
  if (indUpp==ihmax) { // Got the last tested h as upper limit -> WARNING
    std::cout << "WARNING: Upper limit is equal to upper test scan limit!";
    std::cout << "         Increase the scan range to be sure." << std::endl;
  }
  //  std::cout << "Max norm = " << pnormMax << std::endl;
  bool limitsOK = false;
  if (m_foundLower && m_foundUpper) {
    limitsOK = ((fabs(m_lowNorm-1.0)<1e-1) && (fabs(m_uppNorm-1.0)<1.0e-1)); // a bit ugly...
  }
  if (limitsOK) {
    std::cout << std::endl;
    std::cout << "  Limits: [ " << m_lowerLimit << " , " << m_upperLimit << " ]" << std::endl;
    std::cout << "  P(lower limit) = " << m_lowProb << std::endl;
    std::cout << "  P(upper limit) = " << m_uppProb << std::endl;
    std::cout << "  Norm(lower limit) = " << m_lowNorm << "   => should be ~ 1.0" << std::endl;
    std::cout << "  Norm(upper limit) = " << m_uppNorm << "   => should be ~ 1.0" << std::endl;
    //    std::cout << "  s(max) used    = " << m_sMaxUsed << std::endl;
    std::cout << std::endl;
  } else {
    std::cout << std::endl;
    std::cout << "Limits NOT OK!" << std::endl;
    std::cout << "Limits: [ " << m_lowerLimit << " , " << m_upperLimit << " ]" << std::endl;
    std::cout << "ihmax = " << ihmax << std::endl;
    std::cout << "P(lower limit)    = " << m_lowProb << "   => should be <= CL" << std::endl;
    std::cout << "P(upper limit)    = " << m_uppProb << "   => should be <= CL" << std::endl;
    std::cout << "Norm(lower limit) = " << m_lowNorm << "   => should be ~ 1.0" << std::endl;
    std::cout << "Norm(upper limit) = " << m_uppNorm << "   => should be ~ 1.0" << std::endl;
    std::cout << "Found lower = " << m_foundLower << std::endl;
    std::cout << "Found upper = " << m_foundUpper << std::endl;
    //    std::cout << "s(max) used    = " << m_sMaxUsed << std::endl;
    std::cout << std::endl;
  }
}

void Combine::init() {

  tabulateLikelihood();

  m_indexNobs = getIndVector(m_nVecObs,m_nVectors);

  m_bestMu.resize(m_nVectors.size());
  m_bestMuProb.resize(m_nVectors.size());

  m_lhRatio.resize(m_nVectors.size());
}

void Combine::initCorr() {

  makeCorrIntCDF(); //TEMP

  tabulateLikelihoodCorr();

  m_indexNobs = getIndVector(m_nVecObs,m_nVectors);

  m_bestMu.resize(m_nVectors.size());
  m_bestMuProb.resize(m_nVectors.size());

  m_lhRatio.resize(m_nVectors.size());
}

void Combine::doIt() {
  findBestMu();
  calcLimits();
}
