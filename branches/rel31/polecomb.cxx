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
  double epC,depC;
  ep = 1.0;
  epref = 1.0/(0.786/1785.0); // 0.786 = uncorr eff, 1785 = n events
  dep  = ep*0.131; // uncorr
  Measurement CDFcc(0,
		    ep,   dep,  DIST_GAUS,
		    0.81, 0.12, DIST_GAUS);
  //
  // CDFcf: limit(pole) -> [0.00,3.55]
  //
  ep = (696.0/0.485)/epref;
  dep = ep*0.113; // uncorr
  Measurement CDFcf(0,
		    ep,   dep,  DIST_GAUS,
		    0.66, 0.13, DIST_GAUS);
  //
  // An aditional pole object for the 100% correlated uncertainty
  //
  epC  = 1.0;
  depC = 0.155;
  Measurement CDFcorr(0,
		      epC,   depC,  DIST_GAUS,
		      0,     0,     DIST_NONE);

  //
  // DO: limit(pole) ->
  //
  ep = 1.0/6.14e-8/epref; // WRONG
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
  Pole poleCorr;
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
  poleCorr.setMeasurement(CDFcorr); // CDF TEMP: correlated bit
  //
  for (unsigned int i=0; i<1; i++) {
    //    pmSig.generateMeasurement(m[i]);
    poleArr[i].setMeasurement(m[i]);
    poleList.push_back(&poleArr[i]);
  }
  
  for (unsigned int i=0; i<poleList.size()+1; i++) { // TEMP
    if (i<poleList.size()) { // TEMP
      pole = poleList[i];
    } else {
      pole = &poleCorr;
    }
    pole->setPoisson(&PDF::gPoisson);
    pole->setGauss(&PDF::gGauss);
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
    if (i<poleList.size()) combine.add(pole); // TEMP
  }
  //
    combine.setSmax(5.0);
  if (false) {
    combine.setPoleCorr(&poleCorr);
    combine.initCorrelations();
    //    combine.corrEff(poleList[0],poleList[1],0.0);
    combine.initCorr();
  } else {
    combine.init();
  }
  combine.doIt();
}
