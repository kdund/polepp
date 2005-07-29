#include <iostream>
#include <fstream>
#include <iomanip>
#include "Pole.h"
#include "Combine.h"
#include "Measurement.h"
#include "PseudoExperiment.h"

int main(int argc, char *argv[]) {

  PDF::gPoisson.init(50000,100,50);
  PDF::gGauss.init(50000,10.0);
  //
  Measurement exp97(7,
		    0.2, 0.05, DIST_GAUS,
		    3.5, 0.0,  DIST_NONE);
  Measurement exp99(4,
		    0.6, 0.24, DIST_GAUS,
		    3.0, 0.0,  DIST_NONE);
  Measurement exp00(1,
		    1.0, 0.10, DIST_GAUS,
		    1.6, 0.0,  DIST_NONE);

  Measurement CDFcc(0,
		    0.786,0.127,DIST_GAUS,
		    0.87, 0.148,DIST_GAUS);

  Measurement CDFcf(0,
		    0.485,0.099,DIST_GAUS,
		    0.66, 0.197,DIST_GAUS);

  Measurement D0(4,
		 0.247,0.077,DIST_GAUS,
		 4.3,  0.279,DIST_GAUS);

  Measurement m1(1,
		     1.0, 0.1, DIST_GAUS,
		     0.0, 0.0, DIST_NONE );
  Measurement m2(1,
		     1.0,  0.1, DIST_GAUS,
		     0.0,  0.0, DIST_NONE );
  Measurement m3(1,
		     1.0,  0.1, DIST_GAUS,
		     0.0,  0.0, DIST_NONE );

  PseudoExperiment pmSig(m1);
  PseudoExperiment pmBkg(m2);
  pmSig.getRandomGenerator()->setSeed(65539);
  //

  Combine combine;
  //
  std::vector<Pole *> poleList;
  //
  Pole *pole;
  Pole poleArr[10];
  Measurement m[10];
  m[0] = exp97;
  m[1] = exp99;
  m[2] = exp00;
  //
  for (unsigned int i=0; i<3; i++) {
    //    pmSig.generateMeasurement(m[i]);
    if (i!=7) {
      poleArr[i].setMeasurement(m[i]);
      poleList.push_back(&poleArr[i]);
    }
  }
  
  for (unsigned int i=0; i<poleList.size(); i++) {
    pole = poleList[i];
    pole->setPoisson(&PDF::gPoisson);
    pole->setGauss(&PDF::gGauss);
    pole->setCL(0.9);
    pole->setDmus(0.01);
    pole->setBelt(0);
    pole->setEffInt(5.0,21); // integrate +-5 sigma around eff and bkg
    pole->setBkgInt(5.0,21); // integrate +-5 sigma around eff and bkg
    pole->printSetup();
    pole->initAnalysis();
    combine.add(pole);
  }
  //
  combine.init();
  combine.doIt();
}
