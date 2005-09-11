#ifndef POWER_H
#define POWER_H

#include "Pole.h"
#include "Random.h"
#include "Measurement.h"
#include "PseudoExperiment.h"

class Power {
 public:
  Power();
  Power(Pole *pole);
  ~Power();

  void setPole(Pole *pole) { m_pole = pole; if (pole) m_experiment.setMeasurement(pole->getMeasurement()); }
  void setTrueSignal( double smin, double smax, double sstep ) { m_sRange.setRange(smin,smax,sstep); }
  void setTrueSignal( const Range & rng ) { m_sRange = rng; }
  void setHypSignal(double s) { m_sHyp = s; }
  void setLoops(int n) { m_nLoops = n; }

  const double   getMinSignal()   const { return m_sRange.min(); }
  const double   getMaxSignal()   const { return m_sRange.max(); }
  const double   getStepSignal()  const { return m_sRange.step(); }
  const double   getHypSignal()   const { return m_sHyp; }
  const double   getPower()       const { return m_power; }

  //
  bool calculate(double s);
  //
  void resetPower();
  void updatePower();
  const double calcPower() {
    m_power = (m_nTotal>0 ? double(m_nOutside)/double(m_nTotal):0.0);
    m_powerUnc = sqrt(m_power*(1.0-m_power)/static_cast<double>(m_nTotal));
    return m_power;
  }

  void doLoop();

  void outputResult();
  //
 private:
  Pole   *m_pole;
  PseudoExperiment m_experiment;
  Measurement m_measurement;
  int     m_nLoops;
  //
  double  m_sHyp;
  Range   m_sRange;

  double  m_probHyp;
  double  m_probTrue;
  int     m_n1hyp; // ranges of n belt
  int     m_n2hyp;
  int     m_n1true;
  int     m_n2true;

  double m_sumP;
  double m_sumPOutside;
  //
  int     m_nOutside;
  int     m_nTotal;
  double  m_power;
  double  m_powerUnc;


};

#endif
