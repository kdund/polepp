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
  if (argc<5) {
    std::cout << "Usage: estbelt <nobserved> <eff> <sigma(eff)> <bkg> <sigma(bkg)>" << std::endl;
    std::cout << std::endl;
    std::cout << "For a given setup, this tool estimates the limits and the minimum nbelt required using a simple model." << std::endl;
    return 0;
  }
  int nobs=1;
  double eff=1.0;
  double seff=0.0;
  double bkg=0.0;
  double sbkg=0.0;
  //
  nobs=atoi(argv[1]);
  eff =atof(argv[2]);
  seff=atof(argv[3]);
  bkg =atof(argv[4]);
  sbkg=atof(argv[5]);
  //
  double sl = BeltEstimator::getSigLow(nobs,eff,seff,bkg,sbkg);
  double su = BeltEstimator::getSigUp(nobs,eff,seff,bkg,sbkg);
  int nbelt = BeltEstimator::getBeltMin(nobs,eff,seff,bkg,sbkg);
  std::cout << "N_observed      : " << nobs << std::endl;
  std::cout << "Efficiency      : " << eff  << " +- " << seff << std::endl;
  std::cout << "Background      : " << bkg  << " +- " << sbkg << std::endl;
  std::cout << "Estimated limit : [ " << sl << ", " << su << " ]" << std::endl;
  std::cout << "Est. minimum of nBelt : " << nbelt << std::endl;
}
