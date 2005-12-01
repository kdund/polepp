#define PDF_CXX
#include "Pdf.h"

namespace PDF {

  Poisson  gPoisson;
  PoisTab  gPoisTab(&gPoisson);
  Gauss    gGauss;
  GaussTab gGaussTab(&gGauss);

  Gauss2D   gGauss2D;
  LogNormal gLogNormal;
};

