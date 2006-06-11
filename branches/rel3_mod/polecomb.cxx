#include <iostream>
#include <fstream>
#include <iomanip>
#include "Pole.h"
#include "Combine.h"
#include "Measurement.h"

int main(int argc, char *argv[]) {

  if (argc<3) return 0;

  int n0 = atoi(argv[1]);
  int n1 = atoi(argv[2]);

  PDF::gPoisTab.setRangeMean(100000,0,100);
  PDF::gPoisTab.setRangeX(61,0,60);
  PDF::gPoisTab.tabulate();

  MeasPoisEB expA;
  MeasPoisEB expB;
  Combine combine;
  //
  std::vector<Pole *> poleList;
  //
  Pole poleA;
  Pole poleB;
  Pole *pole;

  //  poleA.setVerbose(100);
  //  poleB.setVerbose(100);
  double es = 0.1/sqrt(2.0);
  double bs = 0.3/sqrt(2.0);
  expA.setEffPdf(0.50,es,PDF::DIST_GAUS);
  expA.setBkgPdf(1.50,bs,PDF::DIST_GAUS);
  expA.setNObserved(n0);
  expA.setEffObs();
  expA.setBkgObs();

  expB.setEffPdf(0.50,es,PDF::DIST_GAUS);
  expB.setBkgPdf(1.50,bs,PDF::DIST_GAUS);
  expB.setNObserved(n1);
  expB.setEffObs();
  expB.setBkgObs();
  
  //
  // Set pole measurements
  //
  poleA.setMeasurement(expA);
  poleB.setMeasurement(expB);

  poleList.push_back(&poleA);
  poleList.push_back(&poleB);
  //
  // Loop over all (2) pole objects
  //
  for (unsigned int i=0; i<poleList.size(); i++) {
    pole = poleList[i];
    pole->setPoisson(&PDF::gPoisTab);
    pole->setGauss(&PDF::gGauss);
    pole->setGauss2D(&PDF::gGauss2D);
    pole->setLogNormal(&PDF::gLogNormal);
    pole->setCL(0.90);
    pole->setMinMuProb();
    pole->setDmus(0.01);
    pole->setBelt(20);
    pole->setEffInt(5.0,20); // integrate +-5 sigma around eff and bkg
    pole->setBkgInt(5.0,20); // integrate +-5 sigma around eff and bkg
    pole->setTestHyp(0.0,20.0,0.01); // should be called AFTER the integration construct has been initiated (needs the norm)
    pole->initAnalysis();

    pole->setMethod(RL_FHC2);
    pole->printSetup();
    combine.add(pole);
  }
  //
  // Uncorrelated combination
  //
  combine.setSmax(5.0);
  combine.init();
  combine.doIt();
}
