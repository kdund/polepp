#define OBSERVABLE_CXX
#include "ObsNew.h"

namespace OBS {
  Base *makeObservable(PDFN::DISTYPE dist) {
    Base *obs=0;
    switch (dist) {
    case PDFN::DIST_UNDEF:
      break;
    case PDFN::DIST_NONE:
      break;
    case PDFN::DIST_POIS:
      obs=new ObservablePois();
      obs->setPdf(&PDFN::gPoisson);
      break;
    case PDFN::DIST_GAUS:
      obs=new ObservableGauss();
      obs->setPdf(&PDFN::gGauss);
      break;
    case PDFN::DIST_FLAT:
      //      obs=new ObservableFlat();
      std::cout << "WARNING: Not yet implemented - ObservableFlat()" << std::endl;
      break;
    case PDFN::DIST_LOGN:
      //      obs=new ObservableLogN();
      std::cout << "WARNING: Not yet implemented - ObservableLogN()" << std::endl;
      break;
    default:
      std::cout << "WARNING: Unknown distribution = " << distTypeStr(dist) << std::endl;
      break;
    }
    if (obs) {
      obs->setRndGen(&RND::gRandom);
      obs->validate();
    }
    return obs;
  }
};
