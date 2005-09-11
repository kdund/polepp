#ifndef COMBINE_H
#define COMBINE_H
//
//
//
#include "Combination.h"

class Pole;

class Combine {
 public:
  Combine();
  virtual ~Combine();
  //
  void add(const Pole *point);
  void add(const std::vector<Pole *> poleList);

  bool correlation(const Pole *p1, const Pole *p2, const double corr, bool isBkg);
  bool corrEff(const Pole *p1, const Pole *p2, const double corr) { return correlation(p1,p2,corr,false); }
  bool corrBkg(const Pole *p1, const Pole *p2, const double corr) { return correlation(p1,p2,corr,true); }
  //
  void setSmax(double s) {m_sVecMax = s;}
  ///////////////////////////////////////////
  void initCorrelations(); //TODO: init of correlation matrix - name confusing???
  void makeCorrInt(); // TODO: keep this one but needs to be modified for the general case.
  void initCorr(); // TODO: temp. solution for correlated experiments - to be removed
  ///////////////////////////////////////////
  void init(); // called AFTER all points are added
  void doIt();
  //
  void printResult();
  //
  int makeIndVector(const int ndim,                std::vector< std::vector<int> > & nvec);
  int makeIndVector(const std::vector<int> & ndim, std::vector< std::vector<int> > & nvec);

  double getSmax(double s) {return m_sVecMax;}

  const int findPoint(const Pole *pole) {
    unsigned int i=0;
    bool foundIt = false;
    while ((i<m_poleList.size()) && (!foundIt)) {
      foundIt = (m_poleList[i]==pole);
      i++;
    }
    int ind = -1;
    if (foundIt) ind = i-1;
    return ind;
  }
  //
 private:
  const Pole *m_poleRef; // first Pole object added
  std::vector<const Pole *> m_poleList;
  std::vector<int>          m_nVecObs;
  //
  std::vector< std::vector<double> >       m_corrMat;   // [2*expInd][2*expInd]
  std::vector< std::vector<int> >          m_corrMatFlag;
  //
  std::vector< std::vector<int> >          m_bkgInt;    // [bkgInd][expInd]
  std::vector< std::vector<int> >          m_effInt;
  std::vector< std::vector<double> >       m_weightInt; // [effInd][bkgInd]

  //
  // Vectors for the s(best) scan in n-space
  //
  Range m_bestMuScan;               // scan range in s
  std::vector<double> m_bestMu;     // s(best)
  std::vector<double> m_bestMuProb; // L(n,s(best))
  //
  // Tabulated Likelihood
  //
  Range                              m_sHypRange;
  std::vector<double>                m_sVector;
  std::vector<int>                   m_sMinInd;
  std::vector<int>                   m_sMaxInd;
  std::vector< std::vector<int> >    m_nVectors; // n-vectors used
  std::vector< std::vector<double> > m_likeliHood; // L(n[i],s[j]) = [i][j]
  std::vector<double>                m_lhNorm; // Norm of L for a fixed n
  double m_sVecMax; // max s in sVec (L(n,s)). If < 0 => automatic
  double m_sMaxUsed;// max used sVec == max s(max) encountered

  std::vector<double> m_lhRatio;
  double m_lowerLimit;
  double m_upperLimit;
  bool   m_foundLower;
  bool   m_foundUpper;
  double m_lowNorm;
  double m_uppNorm;
  double m_lowProb;
  double m_uppProb;
  //
  unsigned int m_indexNobs; // index in nVectors of Nobs vector
  //
//   // returns index of a symmetric matrix
//   // vector = [a00,a01,a02...a0(N-1),a11,a12,...,a(N-1)(N-1)]
//   // return -1 if failure
//   //
//   const int getMatrixIndex(int n, int m, int N) const {
//     if (!((n<N) && (m<N))) return -1;
//     int nn,mm;
//     if (n<m) {
//       nn = n;
//       mm = m;
//     } else {
//       nn = m;
//       mm = n;
//     }
//     return (nn*(N-1)+mm-(nn*(nn-1)/2));
//   }
//   //
//   const double getMatrixVal(std::vector<double> & mat, int n, int m, int N) const {
//     int ind = getMatrixIndex(n,m,N);
//     int ms = mat.size();
//     double rval=0.0;
//     if ((ind>-1) && (ind<ms)) {
//       rval = mat[ind];
//     }
//     return rval;
//   }
//   void setMatrixVal(double val, std::vector<double> & mat, int n, int m, int N) const {
//     int ind = getMatrixIndex(n,m,N);
//     int ms = mat.size();
//     if ((ind>-1) && (ind<ms)) {
//       mat[ind]=val;
//     }
//   }
  //
  // TODO: This is a quick solution to demonstrate the combination of 2 experiments with correlated
  // backgrounds.
  //
  void tabulateLikelihoodCorr();
  //
  void tabulateLikelihood();
  double getLikelihood(int nind, double s);
  unsigned int getIndVector(std::vector<int> & nvec, std::vector< std::vector<int> > & allVectors) const;

  bool isOK(const Pole *point);

  bool isFullyCorrelated(double corr, double eps=1e-16) const { return (((fabs(fabs(corr)-1.0)) < eps)); }
  bool isNotCorrelated(double corr, double eps=1e-16)   const { return (fabs(corr) < eps); }

  bool isCorrelated(const Pole *pole, std::vector<const Pole *> poleList, bool doEff);
  bool isEffCorr(const Pole *pole, std::vector<const Pole *> poleList) { return isCorrelated(pole,poleList,true); }
  bool isBkgCorr(const Pole *pole, std::vector<const Pole *> poleList) { return isCorrelated(pole,poleList,false); }

  const double calcProb(std::vector<int> nvec, double s) const;
  
  const bool corrOK(DISTYPE d1, DISTYPE d2) const {
    bool p1ok = ((d1 == DIST_GAUS) || (d1 == DIST_GAUSCORR));
    bool p2ok = ((d2 == DIST_GAUS) || (d2 == DIST_GAUSCORR));
    return (p1ok && p2ok);
  }
  //
  // Returns:
  // sgn(corr)*2 : fully correlated (c=+1.0 or -1.0)
  // sgn(corr)*1 : correlated  (c<>0.0 and abs(c)<1.0)
  //           0 : no correlation
  //
  const int corrFlag(const double corr) const {
    int flag;
    if (isFullyCorrelated(corr)) {
      flag = (corr<0.0 ? -2:2);
    } else {
      if (isNotCorrelated(corr)) {
	flag = 0;
      } else {
	flag = (corr<0.0 ? -1:1);
      }
    }
    return flag;
  }
  //
  void setBestMuScan(int ind);
  void findBestMu(int ind); // input vector must be of same size as the poleList
  void findBestMu();

  double calcLimit(double s);
  void findLimits();

};

#endif
