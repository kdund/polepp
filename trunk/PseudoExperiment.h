#ifndef PSEUDOEXPERIMENT_H
#define PSEUDOEXPERIMENT_H

#include "Random.h"
#include "Measurement.h"

class PseudoExperiment: public Measurement {
 public:
  PseudoExperiment() {};
  PseudoExperiment( const Measurement & m ):m_fixed(false) { copy(m); m_sTrue = m.getSignal(); }

  PseudoExperiment( const PseudoExperiment &exp ) { copy(exp); m_sTrue = exp.getTrueSignal(); m_fixed = exp.isFixed();}
  ~PseudoExperiment() {};

  void setMeasurement( const Measurement & m ) { copy(m); m_fixed = false; m_sTrue = m.getSignal(); }
  void setTrueSignal( double s ) { m_sTrue = s; }
  void setRandomGenerator(RND::Random & rndGen) { m_rnd = rndGen; }
  void setFixed(bool flag=true) { m_fixed=flag; } // true: nobs is just (e*s + b), else nobs = Po(e*s+b)
  //
  const RND::Random & getRandomGenerator() const { return m_rnd; }
  RND::Random *       getRandomGenerator() { return &m_rnd; }
  const double   getTrueSignal() const { return m_sTrue; }
  const bool     isFixed() const { return m_fixed; }
  //
  void generateMeasurement(Measurement & m);
  //
  PseudoExperiment & operator=(const PseudoExperiment & exp) {
    if (this!=&exp) {
      Measurement::copy(exp);
      m_sTrue = exp.getTrueSignal();
      m_rnd   = exp.getRandomGenerator();
    }
    return *this;
  }
  //
  bool operator==(const PseudoExperiment & exp) { return ((static_cast<Measurement>(*this)==static_cast<Measurement>(exp)) && (m_sTrue ==exp.getTrueSignal())); }

 private:
  double m_sTrue;
  RND::Random m_rnd;
  bool   m_fixed;
};

#endif
