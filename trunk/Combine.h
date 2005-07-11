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
  //
  void init(); // called AFTER all points are added
  void doIt();
  //
  void printResult();
  //
 private:
  std::vector<const Pole *> m_poleList;
  std::vector<int> m_nVecObs;
  const Pole *m_poleRef; // first Pole object added
  //
  // Vectors for the s(best) scan in n-space
  //
  Range m_bestMuScan;               // scan range in s
  std::vector<double> m_bestMu;     // s(best)
  std::vector<double> m_bestMuProb; // L(n,s(best))
  //
  // Tabulated Likelihood
  //
  Range m_sRange;
  std::vector<double>                m_sVector;
  std::vector< std::vector<int> >    m_nVectors; // n-vectors used
  std::vector< std::vector<double> > m_likeliHood; // L(n[i],s[j]) = [i][j]

  std::vector<double> m_lhRatio;
  double m_lowerLimit;
  double m_upperLimit;
  bool   m_foundLower;
  bool   m_foundUpper;
  double m_lowNorm;
  double m_uppNorm;
  //
  unsigned int m_indexNobs; // index in nVectors of Nobs vector
  //
  void tabulateLikelihood();
  double getLikelihood(int nind, double s);
  unsigned int getNvecIndex(std::vector<int> & nvec) const;

  bool isOK(const Pole *point);

  void setBestMuScan(int ind);
  void findBestMu(int ind); // input vector must be of same size as the poleList
  void findBestMu();

  double calcLimit(double s);
  void findLimits();

};

#endif
