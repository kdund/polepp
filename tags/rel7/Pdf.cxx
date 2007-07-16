#define PDF_CXX
#include "Pdf.h"

namespace PDF {
  Poisson  gPoisson;
  Poisson  gPoissonTab;
  PoisTab  gPoisTab(&gPoissonTab);
  Gauss    gGauss;
  GaussTab gGaussTab(&gGauss);

  Gauss2D   gGauss2D;
  LogNormal gLogNormal;
  Flat      gFlat;
};

