#ifndef BELTESTIMATOR_H
#define BELTESTIMATOR_H

/*!
 * @class BeltEstimator
 * This code is mainly useful for estimating the number of n in the confidence belt construction (nBelt).
 * The constants used are obtained from studies of the dependencies of the belt on signal, background
 * and efficiency.
 * For the purpose of estimating nBelt, only signal and background are important.
 * Variation in efficiency does not affect the estimate to a large degree.
 */
class BeltEstimator {
public:
  static double fLow(double sig, double bkg) { return m_sl*sig+m_bl*bkg+m_cl; }
  static double fUp(double sig, double bkg)  { return m_su*sig+m_bu*bkg+m_cu; }
  static int    nLow(double sig, double bkg) { return int(0.5+fLow(sig,bkg)); }
  static int    nUp(double sig, double bkg)  { return int(0.5+fUp(sig,bkg)); }
  //
  static double getSigUp(int nobs, double bkg)   { return (double(nobs) - m_cl - m_bl*bkg)/m_sl; }
  static double getSigLow(int nobs, double bkg)  { return (double(nobs) - m_cu - m_bu*bkg)/m_su; }
  static int    getBeltMin(int nobs, double bkg) { return int(0.5+fUp(getSigUp(nobs,bkg),bkg)); }
  //
private:
  static const double m_sl= 0.4;
  static const double m_su= 1.9;
  static const double m_bl= 0.2;
  static const double m_bu= 1.6;
  static const double m_cl=-0.3;
  static const double m_cu= 1.3;
};

#endif
