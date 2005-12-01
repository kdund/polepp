#include <iostream>
#include <fstream>
#include <iomanip>
#include <tclap/CmdLine.h> // Command line parser
#include "Pole.h"


using namespace TCLAP;
void argsPole(Pole *pole, int argc, char *argv[]) {
  //
  try {
    // Create a CmdLine object.
    // Arg1 = string printed at the end whenever --help is used or an error occurs
    // Arg2 = delimiter character between opt and its value
    // Arg3 = version number given when --version is used
    CmdLine cmd("Try again, friend.", ' ', "0.99");

    ValueArg<int>    nObs(      "","nobs",     "number observed events",false,1,"int",cmd);
    ValueArg<double> confLevel( "","cl",       "confidence level",false,0.9,"float",cmd);
    ValueArg<double> sTrue(     "","strue",   "s_true, only used if -C is active",false,1.0,"float",cmd);
    //    SwitchArg        coverage(  "C","coverage", "For coverage studies",false,cmd);
    ValueArg<int>    method(    "m","method",     "method (1 - FHC2 (def), 2 - MBT)",false,1,"int",cmd);
    SwitchArg        useTabulated("T","tab","Use tabulated poisson",false,cmd);
    //
    ValueArg<double> minProb( "","minp",       "minimum probability",false,-1.0,"float",cmd);
    //
    ValueArg<double> effSigma(  "", "effsigma","sigma of efficiency",false,0.2,"float",cmd);
    ValueArg<double> effMeas(   "", "effmeas",  "measured efficiency",false,1.0,"float",cmd);
    ValueArg<int>    effDist(   "", "effdist",  "Efficiency distribution",false,2,"int",cmd);

    //
    ValueArg<double> bkgSigma(  "", "bkgsigma", "sigma of background",false,0.0,"float",cmd);
    ValueArg<double> bkgMeas(   "", "bkgmeas",  "measured background",false,0.0,"float",cmd);
    ValueArg<int>    bkgDist(   "", "bkgdist",  "Background distribution",false,0,"int",cmd);
    //
    ValueArg<double> beCorr(    "", "corr",     "corr(bkg,eff)",false,0.0,"float",cmd);
   //
    ValueArg<double> dMus(      "", "dmus",     "step size in findBestMu",false,0.002,"float",cmd);
    ValueArg<int>    nMus(      "", "nmus",     "maximum number of steps in findBestMu",false,100,"float",cmd);
    ValueArg<int>    belt(      "", "nbelt",    "maximum n for findBestMu" ,false,0,"int",cmd);
    //
    ValueArg<double> hypTestMin( "","hmin",   "hypothesis test min" ,false,0.0,"float",cmd);
    ValueArg<double> hypTestMax( "","hmax",   "hypothesis test max" ,false,35.0,"float",cmd);
    ValueArg<double> hypTestStep("","hstep",  "hypothesis test step" ,false,0.01,"float",cmd);
    //
    ValueArg<double> effIntScale( "","effscale","eff num sigma in integral", false,5.0,"float",cmd);
    ValueArg<int>    effIntN(     "","effn",    "eff: N points in integral", false,21, "int",cmd);
    ValueArg<double> bkgIntScale( "","bscale",  "bkg num sigma in integral", false,5.0,"float",cmd);
    ValueArg<int>    bkgIntN(     "","bkgn",    "bkg: N points in integral", false,21, "int",cmd);

    ValueArg<int>    doVerbose(   "V","verbose", "verbose pole",    false,0,"int",cmd);
    //
    cmd.parse(argc,argv);
    //
    pole->setPoisson(&PDF::gPoisTab);
    pole->setGauss(&PDF::gGauss);
    pole->setGauss2D(&PDF::gGauss2D);
    pole->setLogNormal(&PDF::gLogNormal);

    pole->setMethod(method.getValue());
    pole->setCL(confLevel.getValue());
    pole->setNObserved(nObs.getValue());
    //
    pole->setEffMeas( effMeas.getValue(), effSigma.getValue(), static_cast<PDF::DISTYPE>(effDist.getValue()) );
    pole->setBkgMeas( bkgMeas.getValue(), bkgSigma.getValue(), static_cast<PDF::DISTYPE>(bkgDist.getValue()) );
    pole->checkEffBkgDists();
    pole->setEffBkgCorr(beCorr.getValue());

    pole->setTrueSignal(sTrue.getValue());
    pole->setCoverage(false); //coverage.getValue());

    pole->setDmus(dMus.getValue());
    pole->setNmusMax(nMus.getValue());

    pole->setEffInt(effIntScale.getValue(),effIntN.getValue());
    pole->setBkgInt(bkgIntScale.getValue(),bkgIntN.getValue());
    //
    pole->setBelt(belt.getValue());
    pole->setTestHyp(hypTestMin.getValue(), hypTestMax.getValue(), hypTestStep.getValue());

    pole->setMinMuProb(minProb.getValue());
    //
    if (useTabulated.getValue()) {
      PDF::gPoisTab.setRangeMean(100000,0,100);
      PDF::gPoisTab.setRangeX(61,0,60);
      PDF::gPoisTab.tabulate();
    }
    pole->initAnalysis();
    //    pole->initIntArrays();
    //    pole->initBeltArrays();
    //
    pole->setVerbose(doVerbose.getValue());

  }
  catch (ArgException e) {
    std::cout << "ERROR: " << e.error() << " " << e.argId() << std::endl;
  }
}
