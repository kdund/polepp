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
  int nobs=1;
  double bkg=0.0;
  if (argc>1) {
    nobs=atoi(argv[1]);
  }
  if (argc>2) {
    bkg=atof(argv[2]);
  }
  double sl = BeltEstimator::getSigLow(nobs,bkg);
  double su = BeltEstimator::getSigUp(nobs,bkg);
  int nbelt = BeltEstimator::getBeltMin(nobs,bkg);
  std::cout << "N_observed      : " << nobs << std::endl;
  std::cout << "Background      : " << bkg    << std::endl;
  std::cout << "Estimated limit : [ " << sl << ", " << su << " ]" << std::endl;
  std::cout << "Est. minimum of nBelt : " << nbelt << std::endl;
}
