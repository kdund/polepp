//
#include "PseudoExperiment.h"

void PseudoExperiment::generateMeasurement(Measurement & rval) {
  //
  // A priori, copy this base measurement
  //
  rval = *this;
  //
   // All parameters are forced positive
  //
  bool notOk = true;
  double effMeas=0;
  double bkgMeas=0;
  int nobs=0;
  //
  switch (m_effDist) {
  case DIST_NONE:
    effMeas = m_effMeas;
    break;
  case DIST_GAUSCORR:
    notOk = true;
    while (notOk) {
      double z1 = m_rnd.gaus(0.0,1.0);
      double z2 = m_rnd.gaus(0.0,1.0);
      effMeas = m_effMeas + z1*m_effSigma;
      bkgMeas = m_bkgMeas + m_bkgSigma*(m_beCorr*z1 + m_beCorrInv*z2);
      notOk = ((bkgMeas<0.0) || (effMeas<0.0));
    }
    break;
  case DIST_LOGN:
    effMeas = m_rnd.logNormal(m_effMeas,m_effSigma);
    break;
  case DIST_FLAT:
    effMeas = m_rnd.flat(m_effMeas, m_effSigma);
    break;
  case DIST_GAUS:
    notOk = true;
    while (notOk) {
      effMeas = m_rnd.gaus(m_effMeas,m_effSigma); // measured efficiency
      notOk   = (effMeas<0.0);
    }
    break;
  default: // ERROR STATE
    effMeas = m_effMeas;
    break;
  }
  //
  switch (m_bkgDist) {
  case DIST_NONE:
    bkgMeas = m_bkgMeas;
    break;
  case DIST_LOGN:
    bkgMeas = m_rnd.logNormal(m_bkgMeas,m_bkgSigma);
    break;
  case DIST_FLAT:
    bkgMeas = m_rnd.flat(m_bkgMeas, m_bkgSigma);
    break;
  case DIST_GAUS:
    notOk = true;
    while (notOk) {
      bkgMeas    = m_rnd.gaus(m_bkgMeas,m_bkgSigma); // measured background
      notOk = (bkgMeas<0.0);
    }
  case DIST_GAUSCORR: // already taken care of above
    break;
  default: // ERROR STATE
    bkgMeas = m_bkgMeas;
    break;
  }
  //
  if (m_fixed) {
    nobs = static_cast<int>(m_effMeas*m_sTrue+m_bkgMeas+0.5);
  } else {
    nobs = m_rnd.poisson(m_effMeas*m_sTrue+m_bkgMeas);
  }
  if (nobs<0) nobs=0;
  //
  rval.setNObserved(nobs);
  rval.setEffMeas(effMeas);
  rval.setBkgMeas(bkgMeas);
}
