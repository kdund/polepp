#include <iostream>
#include "Measurement.h"
void Measurement::dump() const {
  std::cout << "-----------MEASUREMENT------------------------\n";
  std::cout << " N observed         : " << getNObserved() << std::endl;
  std::cout << " Efficiency meas    : " << getEffMeas() << std::endl;
  std::cout << " Efficiency sigma   : " << getEffSigma() << std::endl;
  std::cout << " Efficiency dist    : " << distTypeStr(getEffDist()) << std::endl;
  std::cout << " Background meas    : " << getBkgMeas() << std::endl;
  std::cout << " Background sigma   : " << getBkgSigma() << std::endl;
  std::cout << " Background dist    : " << distTypeStr(getBkgDist()) << std::endl;
  std::cout << " Bkg-Eff correlation: " << getBEcorr() << std::endl;
  std::cout << "----------------------------------------------\n";
}
