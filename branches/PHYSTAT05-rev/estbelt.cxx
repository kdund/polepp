#include <iostream>
#include <fstream>
#include <iomanip>
#include "BeltEstimator.h"


//
// This code is mainly useful for estimating the number of n in the confidence belt construction (nBelt).
// The constants used are obtained from studies of the dependencies of the belt on signal, background
// and efficiency.
// For the purpose of estimating nBelt, only signal and background are important.
// Variation in efficiency does not affect the estimate to a large degree.
//

int main(int argc, char *argv[]) {
  if (argc<7) {
    std::cout << "Usage: estbelt <nobserved> <effdist> <eff> <sigma(eff)> <bkgdist> <bkg> <sigma(bkg)>" << std::endl;
    std::cout << std::endl;
    std::cout << "For a given setup, this tool estimates the limits and the minimum nbelt required using a simple model." << std::endl;
    return 0;
  }
  int nobs=1;
  DISTYPE deff = DIST_NONE;
  double eff=1.0;
  double seff=0.0;
  DISTYPE dbkg = DIST_NONE;
  double bkg=0.0;
  double sbkg=0.0;
  //
  nobs=atoi(argv[1]);
  deff =DISTYPE(atoi(argv[2]));
  eff =atof(argv[3]);
  seff=atof(argv[4]);
  dbkg =DISTYPE(atoi(argv[5]));
  bkg =atof(argv[6]);
  sbkg=atof(argv[7]);

  //
  double s  = BeltEstimator::getSC(nobs,eff,bkg);
  double ds = BeltEstimator::getSigmaC(nobs,eff,seff,bkg,sbkg);
  double t  = BeltEstimator::getT(nobs,eff,seff,bkg,sbkg);
  double sl = BeltEstimator::getSigLow(nobs,deff,eff,seff,dbkg,bkg,sbkg);
  double su = BeltEstimator::getSigUp(nobs,deff,eff,seff,dbkg,bkg,sbkg);
  int nbelt = BeltEstimator::getBeltMin(nobs,deff,eff,seff,dbkg,bkg,sbkg);

  std::cout << "N_observed      : " << nobs << std::endl;
  std::cout << "Efficiency      : " << eff  << " +- " << seff << std::endl;
  std::cout << "Background      : " << bkg  << " +- " << sbkg << std::endl;
  std::cout << "Signal          : " << s    << " +- " << ds   << std::endl;
  std::cout << "Test var.       : " << t    << std::endl;
  std::cout << "Estimated limit : [ " << sl << ", " << su << " ]" << std::endl;
  std::cout << "Est. minimum of nBelt : " << nbelt << std::endl;
}
