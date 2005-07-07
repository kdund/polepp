#include <iomanip>
#include <algorithm>
#include <iterator>
#include "Pole.h"
#include "Combine.h"
#include "Permutation.h"

Combine::Combine() {
  m_foundLower=false;
  m_foundUpper=false;
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
  if (point==0) return rval;
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

void Combine::tabulateLikelihood() {
  //
  // * Fill m_nVectors
  // * Calculate a proper range of s
  // * Calculate L(n,s) and save them
  //
  std::cout << "Tabulating likelihood function" << std::endl;
  std::cout << "- Making nVector...." << std::endl;
  int nmeas = m_poleList.size();
  std::vector<int> nvecRange;
  //
  int nrange=m_poleRef->getNBelt();
  nvecRange.resize(nrange);
  m_nVectors.clear();
  // all possible n
  for (int i=0; i<nrange; i++) nvecRange[i]=i;
  // permutate
  AnalysisUtils::Permutation<const std::vector<int> > nvecPerm(&nvecRange,nmeas);

  std::vector< int > jj;

  // save permutations

  while (nvecPerm.get(jj)) {
    m_nVectors.push_back(jj);
  }
  //
  // now make vectors where all elements are equal
  //
  jj.resize(nmeas);
  for (int i=0; i<nrange; i++) {
    for (int m=0; m<nmeas; m++) {
      jj[m] = i;
    }
    m_nVectors.push_back(jj);
  }
  //
  // sort them
  //
  std::sort(m_nVectors.begin(),m_nVectors.end());
  //
  // -------------  nVectors done!
  //
  //
  // Obtain a proper range of s (m_sVector)
  // Use the maximum n and find the estimated s(up) for all points in poleList.
  // Scale the maximum s with some safety-factor
  //
  std::cout << "- Obtaining range of s...." << std::endl;
  const Pole *pole;
  double smax=0;
  double s;
  for (int i=0; i<nmeas; i++) {
    pole = m_poleList[i];
    // s(max) estimate using n(max)
    s = BeltEstimator::getSigUp(pole->getNObserved(), pole->getBkgMeas());
    if (s>smax) smax=s;
  }
  smax *= 1.1; // increase the range by 10% - just to be sure (?)
  double sstep = m_poleRef->getHypTest()->step();
  if (sstep<=0) sstep = 0.01;
  //  sstep=0.1;
  //
  m_sVector.clear();
  m_sRange.setRange(0.0,smax,sstep); // make range
  std::cout << "- Range in s: [0.0, " << smax << "]" << std::endl;
  for (int i=0; i<m_sRange.n(); i++) {
    m_sVector.push_back(m_sRange.getVal(i));
  }
  //
  // -------- sVector done
  //
  //
  // Loop over all (nVectors, sVector) and calulate L(n,s)
  //
  double lh;
  std::cout << "Size of nVectors   : " << m_nVectors.size() << std::endl;
  std::cout << "Size of nVectors[0]: " << m_nVectors[0].size() << std::endl;
  std::cout << "Size of sVector    : " << m_sVector.size() << std::endl;
  std::cout << "Size of poleList   : " << m_poleList.size() << std::endl;


  std::cout << "- Calculate L(n,s)...." << std::endl;
  m_likeliHood.resize(m_nVectors.size());
  for (unsigned int n=0; n<m_nVectors.size(); n++) {
    //    std::cout << " n=" << n << std::endl;
    m_likeliHood[n].resize(m_sVector.size(),1.0);
    for (unsigned int m=0; m<m_sVector.size(); m++) {
      s = m_sVector[m];
      for (unsigned int p=0; p<m_poleList.size(); p++) { // loop over all points
	pole = m_poleList[p];
	lh = pole->calcProb(m_nVectors[n][p],s);     // L(n(i),s(j))
	m_likeliHood[n][m] *=lh;                     // L(n1,s(j))*L(n2,s(j))*..../
	//	std::cout << "n,s,lh = " << p << ": " << m_nVectors[n][p] << ", " << s << ", " << lh << std::endl;
      }
    }
  }
  std::cout << "- Tabulating DONE!" << std::endl;
  //
  // -------- Likelihood construction done
  //
}

unsigned int Combine::getNvecIndex(std::vector<int> & nvec) const {
  unsigned int nmax =  m_nVectors.back().front();
  unsigned int nsize = m_nVectors.front().size();
  if (nvec.size() != nsize) {
    std::cerr << "WARNING: Input vector of faulty size (" << nvec.size() << " <> " << nsize << ")\n";
    return 0;
  }
  //
  unsigned int ofs = 0;
  unsigned int npp=0;
  int nvecInd;
  for (unsigned int i=0; i<nsize; i++) {
    nvecInd = nsize-i-1;
    if (i==0) {
      npp=1;
    } else {
      npp *= (nmax+1);
    }
    if (nvec[nvecInd]>int(nmax)) {
      std::cerr << "WARNING: Too large element in input vector (" << nvec[nvecInd] << " > " << nmax << ")\n";
      std::cerr << "         " << nvecInd << ",  nsize = " << nsize << ",  nmax = " << nmax << std::endl;
      return 0;
    }
    ofs += npp*nvec[nvecInd];
  }
  return ofs;
}

double Combine::getLikelihood(int nind, double s) {
  if (s<0.0) return 0.0;
  if (s>m_sVector.back()) return 0.0;
  if ((nind<0) || (nind>int(m_nVectors.size()))) return 0.0;
  //
  // Establish closest s in m_sVector
  //
  int ind0, ind1;
  int index = int(s/m_sRange.step());
  double sA = m_sRange.getVal(index);
  double s0, s1, ds;
  int ssize = m_likeliHood[0].size(); // size in s
  if (sA>s) {
    ind1 = index;
    ind0 = (index>1 ? index-1:0);
    s0 = m_sRange.getVal(ind0);
    s1 = sA;
  } else {
    ind1 = (index+1<ssize ? index+1:ssize-1);
    ind0 = index;
    s0 = sA;
    s1 = m_sRange.getVal(ind1);
  }
  ds = s1-s0;
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
  //
  for (int i=0; i<nmeas; i++) {
    pole = m_poleList[i];
    s = pole->getDmus() + (static_cast<double>(m_nVectors[nind][i]) - pole->getBkgRangeInt()->min())/pole->getEffMeas();
    if (s>smax) smax=s;
    //
    s = (static_cast<double>(m_nVectors[nind][i]) - pole->getBkgRangeInt()->max())/pole->getEffRangeInt()->max();
    if (s<0.0) s=0.0;
    if (s<smin) smin=s;
    
  }
  m_bestMuScan.setRange(smin,smax,m_poleRef->getDmus());
  //  std::cout << "Set s_best scan range to: [ " << smin << ":" << smax << " ] with ds = " << m_poleRef->getDmus() << std::endl;
}

// void Combine::findBestMu(int ind) {
//   //
//   // Finds s(best) for the set of N(obs) given by m_nVectors.
//   // The index ind is used for chosing the element to fill
//   // in the m_bestMu/Prob vectors
//   //
//   if (m_poleList.size()==0) {
//     std::cerr << "WARNING: No measurements added" << std::endl;
//     return;
//   }
//   if (m_poleList.size()!=m_nVectors.size()) {
//     std::cerr << "WARNING: Size of input vector does not match the #of measurements!" << std::endl;
//     std::cerr << "nvec = " << m_nVectors.size() << std::endl;
//     std::cerr << "pole = " << m_poleList.size() << std::endl;
//     return;
//   }
//   //
//   setBestMuScan(); // set bestMu scan range
//   //
//   const int nmeas = m_poleList.size();
//   //
//   const Pole *pole;
//   std::vector<double> lhoodTot;
//   double sScan, p;
//   //
//   lhoodTot.resize(m_bestMuScan.n(),1.0); // init all elements to 1.0
//   //
//   // loop over all points. Note nBelt should be equal in all cases
//   //
//   // Make: Ltot = L(n1,s)*L(n2,s)*... for each s
//   //
//   for (int i=0; i<nmeas; i++) { // loop over all points
//     pole = m_poleList[i];
//     for (int j=0; j<m_bestMuScan.n(); j++) { // loop over scan range
//       sScan = m_bestMuScan.getVal(j);        // s(j)
//       p = pole->calcProb(m_nVectors[i],sScan);     // L(n(i),s(j))
//       lhoodTot[j] *= p;                      // L(n1,s(j))*L(n2,s(j))*...
//       //      std::cout << " L(" << m_nVectors[i] << ", " << sScan << ") = " << p << "  => tot(s) = " << lhoodTot[j] << std::endl;
//     }
//   }
//   //
//   // Find the s that maximizes the likelihood
//   //
//   double lhMax=-1000.0;
//   double sbest = 0;
//   for (int j=0; j<m_bestMuScan.n(); j++) {
//     if (lhoodTot[j]>lhMax) {
//       lhMax = lhoodTot[j];
//       sbest = m_bestMuScan.getVal(j);
//     }
//   }
// //   std::cout << "SBEST: " << sbest << " ";
// //   for (int i=0; i<nmeas; i++) {
// //     std::cout << m_nVectors[i] << " ";
// //   }
// //   std::cout << std::endl;
//   m_bestMu[ind]     = sbest;
//   m_bestMuProb[ind] = lhMax;
// }

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
}

// void Combine::findBestMu() {
//   // loop over ALL vectors of n....
//   std::cout << "Finding all s_best..." << std::endl;
//   int i;
//   // **** REMOVE - now in tabulateLikelihood()
//   int nmeas = m_poleList.size();
//   std::vector<int> nvecRange;
//   //  nvecRange.resize(m_poleRef->getNBelt());
//   //
//   int nrange=m_poleRef->getNBelt();
//   nvecRange.resize(nrange);
//   m_nVectors.clear();
//   //
//   for (i=0; i<nrange; i++) nvecRange[i]=i;
  
//   AnalysisUtils::Permutation<const std::vector<int> > nvecPerm(&nvecRange,nmeas);

//   std::vector< int > jj;

//   while (nvecPerm.get(jj)) {
//     m_nVectors.push_back(jj);
//   }
//   //
//   // now make vectors where all elements are equal
//   //
//   jj.resize(nmeas);
//   for (i=0; i<nrange; i++) {
//     for (int m=0; m<nmeas; m++) {
//       jj[m] = i;
//     }
//     m_nVectors.push_back(jj);
//   }
//   //
//   std::sort(m_nVectors.begin(),m_nVectors.end());
//   // **** END OF REMOVE
//   //
//   // Loop over all vectors
//   //
//   for (i=0; i<int(m_nVectors.size()); i++) {
// //     std::cout << ind << "  -  ";
// //     for (int m=0; m<nmeas; m++) {
// //       std::cout << (m_nVectors[i][m]) << " ";
// //     }
// //     std::cout << std::endl;
//     findBestMu(i);
//   }
// }

void Combine::findBestMu() {
  // loop over ALL vectors of n....
  std::cout << "Finding all s_best..." << std::endl;

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
  }
  std::cout << "Done!" << std::endl;
}

double Combine::calcLimit(double s) {
  double normProb=0;
  double p;
  for (unsigned n=0; n<m_nVectors.size(); n++) {
    p = getLikelihood(n,s);
    m_lhRatio[n] = p/m_bestMuProb[n];
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
  if (sumProb<m_poleRef->getCL()) {
    if (m_foundLower) {
      m_upperLimit = s;
      m_uppNorm = normProb;
    } else {
      m_lowerLimit = s;
      m_lowNorm = normProb;
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
  std::cout << "Finding limits..." << std::endl;
  m_foundLower = false;
  m_foundUpper = false;
  m_lowerLimit = 0;
  m_upperLimit = 0;
  //
  bool done=false;
  int ismax =  int(m_sVector.size());
  int i=0;
  double stst;
  double pnorm;
  //
  while (!done) {
    stst = m_sVector[i];
    pnorm = calcLimit(stst);
    done = ((i==ismax) ||
	    (m_foundUpper) ||
	    (!m_poleRef->normOK(pnorm)));
    i++;
  }
  bool limitsOK = false;
  if (m_foundLower && m_foundUpper) {
    limitsOK = (m_poleRef->normOK(m_lowNorm) && m_poleRef->normOK(m_uppNorm)); // a bit ugly...
  }
  if (limitsOK) {
    std::cout << "Limits: [ " << m_lowerLimit << " , " << m_upperLimit << " ]" << std::endl;
  } else {
    std::cout << "Limits not OK!" << std::endl;
    std::cout << "Limits: [ " << m_lowerLimit << " , " << m_upperLimit << " ]" << std::endl;
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
