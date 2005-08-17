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
  if (point==0) return rval;

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
  m_nVectors.clear();
  //
  std::vector< int > jj;
  jj.resize(nmeas,0);
  m_nVectors.push_back(jj);  
  while (Combination::next_vector(jj,nrange)) {
    m_nVectors.push_back(jj);
  }
  //
  // sort them
  //
  std::sort(m_nVectors.begin(),m_nVectors.end());

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
    // s(max) estimate using n(max)
//     s = BeltEstimator::getSigUp(pole->getNObserved(),
// 				pole->getEffDist(),pole->getEffMeas(),pole->getEffSigma(),
// 				pole->getBkgDist(),pole->getBkgMeas(),pole->getBkgSigma(),
// 				pole->getIntNorm());
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
  int nok=0;
  int ntot=0;
  double plh,lh;
  int ismin,ismax; // indecis of min and max s for the current n
  bool isminOK;
  const double plhLim = 0.0001; // min Lh
  std::vector<double> tmpLh;  // temporary storage
  tmpLh.resize(m_sVector.size(),0.0);
  //
  for (unsigned int n=0; n<m_nVectors.size(); n++) {
    if (n%1000==0) std::cout << "  N vector " << n << " out of " << m_nVectors.size() << std::endl;
    m_sMinInd[n] = -1;
    m_sMaxInd[n] = -1;
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
    for (int m=ismin; m<ismax; m++) {
      m_likeliHood[n].push_back(tmpLh[m]);
      nok++;
    }
    m_sMinInd[n] = ismin;
    m_sMaxInd[n] = ismax;
  }
  std::cout << "- Actual number of L(n,s) calculated = " << nok << std::endl;
  std::cout << "- Tabulating DONE!" << std::endl;
  //
  // -------- Likelihood construction done
  //
}

unsigned int Combine::getNvecIndex(std::vector<int> & nvec) const {
  std::vector<int> nmax =  m_nVectors.back();     // vector containg the maximum n(belt) per Pole point
  unsigned int nsize = m_nVectors.front().size(); // number of points
  unsigned int nFullSize = m_nVectors.size();
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
  //    std::cout << "s, s range: " << s << " => " << smin << "," << smax << std::endl;
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
  ds = s1-s0;
  ind0 = ind0 - ismin;
  ind1 = ind1 - ismin;
  if (ind1>=ssize) ind1 = ssize-1;
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
//   std::cout << "N(ind) = " << nind << std::endl;
//   std::cout << "s0     = " << s0 << std::endl;
//   std::cout << "s1     = " << s1 << std::endl;
//   std::cout << "s      = " << s << std::endl;
//   std::cout << "lh0    = " << lh0 << std::endl;
//   std::cout << "lh1    = " << lh1 << std::endl;
//   std::cout << "lh     = " << rval << std::endl;
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
  for (unsigned n=0; n<m_nVectors.size(); n++) {
    p = getLikelihood(n,s);
    if (m_bestMuProb[n]>0.0) {
      m_lhRatio[n] = p/m_bestMuProb[n];
    } else {
      m_lhRatio[n] = 0.0; //(n==m_indexNobs ? 0.0:0.001); // this is just to make sure that CHECK!!!!
    }
    normProb += p;
  }
  //
  bool done=false;
  unsigned int n=0;
  double sumProb=0;
  while (!done) {
    if (n!=m_indexNobs) {
      if (m_lhRatio[n]>m_lhRatio[m_indexNobs]) {
	sumProb += getLikelihood(n,s);
      }
    }
    n++;
    done = ((n==m_lhRatio.size()) || (sumProb>m_poleRef->getCL()));
  }
  //
  double dCL = sumProb-m_poleRef->getCL();
  //  std::cout << "findLimits(): " << s << " norm = " << normProb << " sump = " << sumProb << std::endl;
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
  return normProb;
}

void Combine::findLimits() {
  std::cout << "- Finding limits..." << std::endl;
  m_foundLower = false;
  m_foundUpper = false;
  m_lowerLimit = 0;
  m_upperLimit = 0;
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
	    (m_foundUpper) ||
	    (!m_poleRef->normOK(pnorm)));
    if (done) indUpp=i; 
  }
  if (indUpp==ihmax) { // Got the last tested h as upper limit -> WARNING
    std::cout << "WARNING: Upper limit is equal to upper test scan limit!";
    std::cout << "         Increase the scan range to be sure." << std::endl;
  }
  //  std::cout << "Max norm = " << pnormMax << std::endl;
  bool limitsOK = false;
  if (m_foundLower && m_foundUpper) {
    limitsOK = (m_poleRef->normOK(m_lowNorm) && m_poleRef->normOK(m_uppNorm)); // a bit ugly...
  }
  if (limitsOK) {
    std::cout << std::endl;
    std::cout << "  Limits: [ " << m_lowerLimit << " , " << m_upperLimit << " ]" << std::endl;
    std::cout << "  P(lower limit) = " << m_lowProb << std::endl;
    std::cout << "  P(upper limit) = " << m_uppProb << std::endl;
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

  m_indexNobs = getNvecIndex(m_nVecObs);

  m_bestMu.resize(m_nVectors.size());
  m_bestMuProb.resize(m_nVectors.size());

  m_lhRatio.resize(m_nVectors.size());
}

void Combine::doIt() {
  findBestMu();
  findLimits();
}
