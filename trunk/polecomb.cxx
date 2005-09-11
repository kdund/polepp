#include <iostream>
#include <fstream>
#include <iomanip>
#include "Pole.h"
#include "Combine.h"
#include "Measurement.h"
#include "PseudoExperiment.h"

int main(int argc, char *argv[]) {

  PDF::gPoisson.init(100000,100,200);
  //  PDF::gGauss.init(50000,10.0);
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
  //
  // CDFcc: limit(pole) -> 
  //
  double ep,dep,epref;
  ep = 1.0/9.9e-8;
  epref = ep;
  ep = ep/epref;
  dep = ep*0.131; // uncorr
  //dep = ep*0.119; // corr
  //dep = ep*0.177; // total
  Measurement CDFcc(0,
		    ep,   dep,  DIST_GAUS,
		    0.81, 0.12, DIST_GAUS);
  //
  // CDFcf: limit(pole) -> [0.00,3.55]
  //
  ep = 1.0/1.57e-7/epref;
  dep = ep*0.113; // uncorr
  //dep = ep*0.119; // corr
  //dep = ep*0.164; // total
  Measurement CDFcf(0,
		    ep,   dep,  DIST_GAUS,
		    0.66, 0.13, DIST_GAUS);

  //
  // DO: limit(pole) ->
  //
  ep = 1.0/6.14e-8/epref;
  dep = ep*0.092;// uncorr
  //dep = ep*0.119; // corr
  //dep = ep*0.150; // total
  Measurement D0(4,
		 ep,   dep, DIST_GAUS,
		 4.3,  1.2, DIST_GAUS);

  Measurement m1(4,
		 1.0, 0.0, DIST_NONE,
		 3.0, 0.0, DIST_NONE );

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
  m[0] = CDFcc;
  m[1] = CDFcf;
  m[2] = D0;
//   m[0] = exp97;
//   m[1] = exp99;
//     m[2] = exp00;
  //
  for (unsigned int i=0; i<2; i++) {
    //    pmSig.generateMeasurement(m[i]);
    if (i!=42) {
      poleArr[i].setMeasurement(m[i]);
      poleList.push_back(&poleArr[i]);
    }
  }
  
  for (unsigned int i=0; i<poleList.size(); i++) {
    pole = poleList[i];
    pole->setPoisson(&PDF::gPoisson);
    pole->setGauss(&PDF::gGauss);
    pole->setCL(0.95);
    pole->setDmus(0.01);
    pole->setBelt(30);
    pole->setEffInt(5.0,21); // integrate +-5 sigma around eff and bkg
    pole->setBkgInt(5.0,21); // integrate +-5 sigma around eff and bkg
    pole->setTestHyp(0.0,3.0,0.01); // should be called AFTER the integration construct has been initiated (needs the norm)
    pole->initAnalysis();

    pole->setNLR(false);
    pole->printSetup();
    combine.add(pole);
  }
  //
  //  combine.initCorrelations();
  //  combine.corrEff(poleList[0],poleList[1],1.0);
  combine.setSmax(5.0);
  combine.init();
  //  combine.initCorr();
  combine.doIt();
}
