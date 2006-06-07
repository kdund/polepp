#include <iostream>
#include <fstream>
#include <iomanip>
#include "Pole.h"
#include "Combine.h"
#include "Measurement.h"

int main(int argc, char *argv[]) {

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

  expA.setEffPdf(1.0,0.1,PDF::DIST_GAUS);
  expA.setBkgPdf(0.0,0.0,PDF::DIST_NONE);
  expA.setNObserved(0);
  expA.setEffObs();
  expA.setBkgObs();

  expB.setEffPdf(1.0,0.1,PDF::DIST_GAUS);
  expB.setBkgPdf(0.0,0.0,PDF::DIST_NONE);
  expB.setNObserved(0);
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
    pole->setCL(0.95);
    pole->setMinMuProb();
    pole->setDmus(0.01);
    pole->setBelt(20);
    pole->setEffInt(5.0,20); // integrate +-5 sigma around eff and bkg
    pole->setBkgInt(5.0,20); // integrate +-5 sigma around eff and bkg
    pole->setTestHyp(0.0,35.0,0.01); // should be called AFTER the integration construct has been initiated (needs the norm)
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
