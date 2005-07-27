#include <iostream>
#include <fstream>
#include <iomanip>
#include <tclap/CmdLine.h> // Command line parser
#include "Pole.h"


using namespace TCLAP;
void processArgs(Pole *pole, int argc, char *argv[]) {
  //
  try {
    // Create a CmdLine object.
    // Arg1 = string printed at the end whenever --help is used or an error occurs
    // Arg2 = delimiter character between opt and its value
    // Arg3 = version number given when --version is used
    CmdLine cmd("Try again, friend.", ' ', "0.99");

    ValueArg<int>    nObs(      "","nobs",     "number observed events",false,1,"int");
    ValueArg<double> confLevel( "","cl",       "confidence level",false,0.9,"float");
    ValueArg<double> sTrue(     "","strue",   "s_true, only used if -C is active",false,1.0,"float");
    SwitchArg        coverage(  "C","coverage", "For coverage studies",false);
    SwitchArg        doNLR("N","nlr", "Use NLR",false);
    SwitchArg        useTabulated("T","tab","Use tabulated poisson",false);
    //
    ValueArg<double> effSigma(  "", "esigma","sigma of efficiency",false,0.2,"float");
    ValueArg<double> effMeas(   "", "emeas",  "measured efficiency",false,1.0,"float");
    ValueArg<int>    effDist(   "","effdist",  "Efficiency distribution",false,1,"int");

    //
    ValueArg<double> bkgSigma(  "", "bsigma","sigma of background",false,0.0,"float");
    ValueArg<double> bkgMeas(   "", "bmeas",   "measured background",false,0.0,"float");
    ValueArg<int>    bkgDist(   "", "bkgdist",  "Background distribution",false,0,"int");
    //
    ValueArg<double> beCorr(    "","corr",    "corr(bkg,eff)",false,0.0,"float");
   //
    ValueArg<double> dMus(      "","dmus",    "step size in findBestMu",false,0.002,"float");
    ValueArg<int>    belt(   "","belt", "maximum n for findBestMu" ,false,0,"int");
    //
    ValueArg<double> hypTestMin( "","hmin",   "hypothesis test min" ,false,0.0,"float");
    ValueArg<double> hypTestMax( "","hmax",   "hypothesis test max" ,false,35.0,"float");
    ValueArg<double> hypTestStep("","hstep",  "hypothesis test step" ,false,0.01,"float");
    //
    ValueArg<double> effIntScale( "","escale","eff n sigma in integral", false,5.0,"float");
    ValueArg<double> effIntStep("","estep",   "eff step in integral",    false,-1.0,"float");
    ValueArg<double> bkgIntScale( "","bscale","bkg n sigma in integral", false,5.0,"float");
    ValueArg<double> bkgIntStep("","bstep",   "bkg step in integral",    false,-1.0,"float");

    ValueArg<int>    doVerbose(   "V","verbose", "verbose pole",    false,0,"int");
    //
    cmd.add(doVerbose);
    cmd.add(doNLR);
    cmd.add(useTabulated);
    //
    cmd.add(hypTestMin);
    cmd.add(hypTestMax);
    cmd.add(hypTestStep);

    cmd.add(effIntScale);
    cmd.add(effIntStep);
    cmd.add(bkgIntScale);
    cmd.add(bkgIntStep);

    cmd.add(belt);
    cmd.add(dMus);

    cmd.add(effSigma);
    cmd.add(effMeas);
    cmd.add(effDist);

    cmd.add(bkgSigma);
    cmd.add(bkgMeas);
    cmd.add(bkgDist);

    cmd.add(beCorr);

    cmd.add(sTrue);
    cmd.add(coverage);
    cmd.add(confLevel);
    cmd.add(nObs);

    //
    cmd.parse(argc,argv);
    //
    pole->setPoisson(&PDF::gPoisson);
    pole->setGauss(&PDF::gGauss);
    pole->setNLR(doNLR.getValue());
    pole->setCL(confLevel.getValue());
    pole->setNobserved(nObs.getValue());
    //
    pole->setEffMeas( effMeas.getValue(), effSigma.getValue(), static_cast<DISTYPE>(effDist.getValue()) );
    pole->setBkgMeas( bkgMeas.getValue(), bkgSigma.getValue(), static_cast<DISTYPE>(bkgDist.getValue()) );
    pole->checkEffBkgDists();
    pole->setEffBkgCorr(beCorr.getValue());

    pole->setTrueSignal(sTrue.getValue());
    pole->setCoverage(coverage.getValue());

    pole->setDmus(dMus.getValue());
    pole->setEffInt(effIntScale.getValue(),effIntStep.getValue());
    pole->setBkgInt(bkgIntScale.getValue(),bkgIntStep.getValue());
    //
    pole->setBelt(belt.getValue()); // call after nObserved is set.
    pole->setBeltMax(belt.getValue()*2); // maximum allocated
    pole->setTestHyp(hypTestMin.getValue(), hypTestMax.getValue(), hypTestStep.getValue());
    //
    if (useTabulated.getValue()) {
      pole->initPoisson(50000,100,50);
      pole->initGauss(50000,10.0);
    }
    pole->initIntArrays();
    pole->initBeltArrays();
    //
    pole->setVerbose(doVerbose.getValue());

  }
  catch (ArgException e) {
    std::cout << "ERROR: " << e.error() << " " << e.argId() << std::endl;
  }
}

int main(int argc, char *argv[]) {

  Pole pole;
  //
  processArgs(&pole, argc, argv);
  //
  if (pole.checkParams()) {
    pole.printSetup();
    if (pole.analyseExperiment()) {
      pole.printLimit(true);
    } else {
      std::cout << "ERROR: limit calculation failed - nbelt is probably too small ("
		<< pole.getNBelt() << ") for Nobs = " << pole.getNObserved() << std::endl;
      std::cout << "       probability    = " << pole.getSumProb() << std::endl;
      std::cout << "       lower lim norm = " << pole.getLowerLimitNorm() << std::endl;
      std::cout << "       upper lim norm = " << pole.getUpperLimitNorm() << std::endl;
    }
  }
}
