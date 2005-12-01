#ifndef BELTESTIMATOR_H
#define BELTESTIMATOR_H

#include "Measurement.h"
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
  static double getSC(int n, double e, double b) { // same, but will return 1/e when s==0
    double ss =  (double(n)-b);
    if (ss<=0.0) return 1.0/e;
    return (ss/e);
  }
  // uncertainty^2 on s through error propagation (n is poissonian)
  static double getSigma2(int n, double e, double se, double b, double sb) {
    double s = getS(n,e,b);
    return (double(n)+sb*sb+s*s*se*se)/(e*e);
  }
  // ditto sqrt()
  static double getSigma(int n, double e, double se, double b, double sb) {
    return sqrt(getSigma2(n,e,se,b,sb));
  }
  static double getSigmaC(int n, double e, double se, double b, double sb) {
    double s  = getSC(n,e,b);
    double ss2 = (double(n)+sb*sb+s*s*se*se)/(e*e);
    return sqrt(ss2);
  }
  // test variable t = s+sigma(s)
  // this seems to be a good variable to use when estimating the boundaries of the limits
  // ..well, not very:
  // * fails when the efficiency distribution is cut off at 0 and the norm is <1
  //   -> ie for large uncertainties
  // * large n
  // The second getT() is a better one.
  //
//   static double getT(double s, double ds, double e) {
//     if (s<=0.0) s=1.0/e;
//     return (s+ds); // removed /e
//   }
//   static double getT(double s, double ds, double e, double se, double norm=1.0) {
//     if (s==0.0) s=1.0/e;
//     double nf;
//     // Take care of cases with large se/e
//     // nf = 1 if norm=1 or se=0
//     // nf ~ 1/se if norm<1
//     // When the normalisation gets
//     // for most cases when se/e is small, norm=1 and nf will be unity
//     // for large se/e, nf will be scaled up with 1/se
//     if (se>0.0) {
//       const double sc=2.0;
//       nf = sc*(fabs(1.0-norm)/se) + norm;
//     } else {
//       nf = 1.0;
//     }
//     //
//     // In order to smooth the dependance of n, another factor is introduced:
//     // ncorr = sqrt(n-b)
//     //
//     double ncorr = double(n) - b;
//     if ((n==0) || (ncorr<0.5)) ncorr=1.0;    
//     return nf*(s+ds)/sqrt(ncorr);
//   }
  //
  static double getT(int n, double e, double se, double b, double sb, double norm=1.0) {
    double s  = getSC(n,e,b);
    double ds = getSigmaC(n,e,se,b,sb);
    //
    double nf;
    // Take care of cases with large se/e
    // nf = 1 if norm=1 or se=0
    // nf ~ 1/se if norm<1
    // When the normalisation gets
    // for most cases when se/e is small, norm=1 and nf will be unity
    // for large se/e, nf will be scaled up with 1/se
    if (se>0.0) {
      const double sc=4.0;
      nf = sc*(fabs(1.0-norm)/se) + norm;
    } else {
      nf = 1.0;
    }
    //
    // In order to smooth the dependance of n, another factor is introduced:
    // ncorr = sqrt(n-b)
    //
    double ncorr = double(n) - b;
    if (ncorr<0.5) {
      ncorr=1.0;
    } else {
      ncorr = pow(ncorr,-0.37); // optimum squeeze for UL estimation
    }
    return nf*(s+ds)*ncorr;

  }
  //
  // Lower/upper limit boundaries
  //
  // Lower: seems ok for most cases - might overestimate when se/e and s+sigma(s) are large.
  static double getSigLow(int n, double e, double se, double b, double sb, double norm=1.0) {
    //    double r=m_al1*getT(n,e,se,b,sb,norm) + m_bl1;
    //    return (r<0.0 ? 0.0:r*e);
    return 0.0;
  }
  static double getSigLow(int n, PDF::DISTYPE de, double e, double se, PDF::DISTYPE db, double b, double sb, double norm=1.0) {
    double rval = getSigLow(n,e,se,b,sb,norm);
    return rval;
  }
  //
  // Upper limit estimation
  //
//   static double getSigUp(int n, double e, double se, double b, double sb, double norm=1.0) {
//     return (m_au1*getT(n,e,se,b,sb,norm)) + m_bu1;
//   }
  //
  // New upper limit estimator
  //
  static double getSigUp(int n, PDF::DISTYPE de, double e, double se, PDF::DISTYPE db, double b, double sb, double norm=1.0) {
    double a0,a1,a2;
    double t,ul=0;
    double deltaN = double(n)-b;
    const double feLim=1.5;
    //
    switch (de) {
    case PDF::DIST_LOGN:
//       a0 = 2.0;
//       a1 = 0.5;
//       a2 = exp(-0.7*double(n)+0.7);
//       t = getT(n,e,se,b,sb,norm);
//       ul = a2*t*t+a1*t+a0;
//       break;
    case PDF::DIST_GAUS:
    default:
      double fe = (se>0.0 ? e/se:10.0);
      if (fe>8.0) {
	a0 = 0.0;
	a1 = 3.9;
	a2 = 0.0;
      } else {
	a0 = 2.0;
	a1 = 2.0;
	a2 = 0.7;
      }
      if ((de==PDF::DIST_GAUS) && (fe<feLim) && (deltaN<0.5)) {
	t = getT(n,feLim*se,se,b,sb,norm); // for e/se<1.5, the upper limit is ~ constant for all t (which are large)
      } else {
	t = getT(n,e,se,b,sb,norm);
	if (t>6.0) {
	  a0 = -15.0;
	  a1 = 10.0;
	  a2 = 0.0;
	}
      }
      ul = (a2*t*t+a1*t+a0);
    }
    return ul;
  }
  //
  // Belt estimations <OLD>
  //
//   static double fLow(double sig, double bkg)  { return m_sl*sig+m_bl*bkg+m_cl; }
//   static double fUp(double sig, double bkg)  { return m_su*sig+m_bu*bkg+m_cu; }
//   static int    nLow(double sig, double bkg) { return int(0.5+fLow(sig,bkg)); }
//   static int    nUp(double sig, double bkg)  { return int(0.5+fUp(sig,bkg)); }
//   //
//   static int    getBeltMin(int n, double e, double se, double b, double sb, double norm=1.0) {
//     return int(0.5+fUp(getSigUp(n,e,se,b,sb,norm),b));
//   }
  static int    getBeltMin(int n, PDF::DISTYPE de, double e, double se, PDF::DISTYPE db, double b, double sb, double norm=1.0) {
    // Belt estimations vary with distributions
    // Log-normal limits are more narrow -> belt rises more steeply -> larger nBelt is required
    //
    double napp=20;
    double deCorr = (de>0 ? pow(se,0.75):1.0);
    double dN = double(n)-b;
    bool nZero = (dN<0.5);
    double t = getT(n,e,se,b,sb,norm)*deCorr;
    double a0,a1,a2;
    if (de==PDF::DIST_LOGN) {
      if (nZero) {
	a2 = 0.7;
	a1 = 3.0;
	a0 = 0.25;
      } else {
	a2 = exp(-0.22*dN+1.45);
	a1 = 1.0;
	a0 = exp( 0.22*dN-1.37);
      }
      napp = a0+a1*t+a2*t*t;
    } else {
      double fe = (se>0.0 ? e/se:10.0);
      if (fe>8.0) {
	napp = 46.0*t+20.0;
      } else {
	napp = 28.0*t+30.0;
      }
    }
    return int(napp+0.5);
  }
  //
private:
  // consts for estimation of upper limit
//   static const double m_au1= 1.2;
//   static const double m_bu1= 3.0;
//   static const double m_au2= 2.5;
//   static const double m_bu2= 3.0;
//   static const double m_au3= 3.5;
//   static const double m_bu3= 3.0;
//   // ditto for lower limit
//   static const double m_al1=  0.2;
//   static const double m_bl1= -1.0;
//   // for N(belt)
//   static const double m_sl= 0.4;
//   static const double m_su= 1.9;
//   static const double m_bl= 0.2;
//   static const double m_bu= 1.6;
//   static const double m_cl=-0.3;
//   static const double m_cu= 1.3;

};

#endif
