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
  // signal from n=e*s+b
  static double getS(int n, double e, double b) { return (double(n)-b)/e; }
  // uncertainty^2 on s through error propagation (n is poissonian)
  static double getSigma2(int n, double e, double se, double b, double sb) {
    double s = getS(n,e,b);
    return (double(n)+sb*sb+s*s*se*se)/(e*e);
  }
  // ditto sqrt()
  static double getSigma(int n, double e, double se, double b, double sb) {
    return sqrt(getSigma2(n,e,se,b,sb));
  }
  // test variable t = s+sigma(s)
  // this seems to be a good variable to use when estimating the boundaries of the limits
  static double getT(int n, double e, double se, double b, double sb) {
    return getS(n,e,b)+getSigma(n,e,se,b,sb);
  }
  //
  // Lower/upper limit boundaries
  //
  // Lower: seems ok for most cases - might overestimate when se/e and s+sigma(s) are large.
  static double getSigLow(int n, double e, double se, double b, double sb) {
    double r=m_al1*getT(n,e,se,b,sb) + m_bl1;
    return (r<0.0 ? 0.0:r);
  }
  // Upper: a tighter boundary is used for when se/e is low (<0.2)
  // No separate classification is needed for cases with large sb/b.
  static double getSigUp(int n, double e, double se, double b, double sb) {
    if (se/e > 0.2) return m_au2*getT(n,e,se,b,sb) + m_bu2;
    return m_au1*getT(n,e,se,b,sb) + m_bu1;
  }
  //
  // Belt estimations
  //
  static double fLow(double sig, double bkg)  { return m_sl*sig+m_bl*bkg+m_cl; }
  static double fUp(double sig, double bkg)  { return m_su*sig+m_bu*bkg+m_cu; }
  static int    nLow(double sig, double bkg) { return int(0.5+fLow(sig,bkg)); }
  static int    nUp(double sig, double bkg)  { return int(0.5+fUp(sig,bkg)); }
  //
  //  static double getSigUp(int nobs, double bkg)   { return (double(nobs) - m_cl - m_bl*bkg)/m_sl; } //OLD
  //  static double getSigLow(int nobs, double bkg)  { return (double(nobs) - m_cu - m_bu*bkg)/m_su; }
  static int    getBeltMin(int n, double e, double se, double b, double sb) { return int(0.5+fUp(getSigUp(n,e,se,b,sb),b)); }
  //
private:
  // consts for estimation of upper limit
  static const double m_au1= 1.2;
  static const double m_bu1= 3.0;
  static const double m_au2= 2.2;
  static const double m_bu2= 3.0;
  // ditto for lower limit
  static const double m_al1=  0.2;
  static const double m_bl1= -1.0;
  // for N(belt)
  static const double m_sl= 0.4;
  static const double m_su= 1.9;
  static const double m_bl= 0.2;
  static const double m_bu= 1.6;
  static const double m_cl=-0.3;
  static const double m_cu= 1.3;

};

#endif
