#include <iostream>
#include <fstream>
#include <iomanip>
#include <tclap/CmdLine.h> // Command line parser
#include "Coverage.h"


using namespace TCLAP;
void argsCoverage(Coverage *coverage, Pole *pole, int argc, char *argv[]) {
  time_t timer;
  time(&timer);
  unsigned int rndSeed=static_cast<unsigned int>(timer);
  //
  try {
    // Create a CmdLine object.
    // Arg1 = string printed at the end whenever --help is used or an error occurs
    // Arg2 = delimiter character between opt and its value
    // Arg3 = version number given when --version is used
    CmdLine cmd("Try again, friend.", ' ', "0.99");

    ValueArg<int>    nLoops(    "","nloops",  "number of loops",    false,1,"int",cmd);

    ValueArg<int>    rSeed(     "","rseed",   "rnd seed" ,          false,rndSeed,"int",cmd);
    ValueArg<int>    rSeedOfs(  "","rseedofs","rnd seed offset" ,   false,0,"int",cmd);

    ValueArg<double> confLevel( "","cl",       "confidence level",false,0.9,"float",cmd);

    ValueArg<int>    method(    "m","method",     "method (1 - FHC2 (def), 2 - MBT)",false,1,"int",cmd);
    SwitchArg        noTabulated("K","notab","Do not use tabulated poisson",false,cmd);
    //
    SwitchArg        doStats("C","stats", "Collect statistics",false,cmd);
    SwitchArg        doFixSig("S","fixsig", "Fixed meas. n_observed",false,cmd);

    ValueArg<double> minProb( "","minp",       "minimum probability",false,-1.0,"float",cmd);
    //
    ValueArg<std::string> dump("","dump",    "dump filename",false,"","string",cmd);

    ValueArg<int>    verboseCov(   "V","verbcov", "verbose coverage",false,0,"int",cmd);
    ValueArg<int>    verbosePol(   "P","verbpol", "verbose pole",    false,0,"int",cmd);

    ValueArg<double> sMin(      "","smin",    "min s_true",         false,1.0,"float",cmd);
    ValueArg<double> sMax(      "","smax",    "max s_true",         false,1.0,"float",cmd);
    ValueArg<double> sStep(     "","sstep",   "step s_true",        false,1.0,"float",cmd);

    ValueArg<int>    effDist(   "", "effdist",  "Efficiency distribution",false,2,"int",cmd);
    ValueArg<double> effSigma(  "","esigma",  "sigma of efficiency",false,0.2,"float",cmd);
    ValueArg<double> effMin(    "","emin",    "min eff true",       false,1.0,"float",cmd);
    ValueArg<double> effMax(    "","emax",    "max eff true",       false,1.0,"float",cmd);
    ValueArg<double> effStep(   "","estep",   "step eff true",      false,1.0,"float",cmd);

    ValueArg<int>    bkgDist(   "", "bkgdist",  "Background distribution",false,0,"int",cmd);
    ValueArg<double> bkgSigma(  "","bsigma",  "sigma of background",false,0.2,"float",cmd);
    ValueArg<double> bkgMin(    "","bmin",    "min bkg true",false,0.0,"float",cmd);
    ValueArg<double> bkgMax(    "","bmax",    "max bkg true",false,0.0,"float",cmd);
    ValueArg<double> bkgStep(   "","bstep",   "step bkg true",false,1.0,"float",cmd);
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
    //
    cmd.parse(argc,argv);
    //
    // First set Pole
    //
    pole->setPoisson(&PDF::gPoisTab);
    pole->setGauss(&PDF::gGauss);
    pole->setGauss2D(&PDF::gGauss2D);
    pole->setLogNormal(&PDF::gLogNormal);

    pole->setMethod(method.getValue());
    pole->setVerbose(verbosePol.getValue());
    pole->setNObserved(0);
    pole->setCoverage(true);
    pole->setCL(confLevel.getValue());
    //
    pole->setDmus(dMus.getValue());
    pole->setNmusMax(nMus.getValue());
    //
    pole->setEffMeas( effMin.getValue(), effSigma.getValue(), static_cast<PDF::DISTYPE>(effDist.getValue()) );
    pole->setBkgMeas( bkgMin.getValue(), bkgSigma.getValue(), static_cast<PDF::DISTYPE>(bkgDist.getValue()) );
    pole->checkEffBkgDists();
    pole->setEffBkgCorr(beCorr.getValue());

    pole->setTrueSignal(sMin.getValue());

    pole->setEffInt(effIntScale.getValue(),effIntN.getValue());
    pole->setBkgInt(bkgIntScale.getValue(),bkgIntN.getValue());
    //
    pole->setBelt(belt.getValue());
    pole->setTestHyp(hypTestMin.getValue(), hypTestMax.getValue(), hypTestStep.getValue());

    pole->setMinMuProb(minProb.getValue());
    //
    if (!noTabulated.getValue()) {
      PDF::gPoisTab.setRangeMean(100000,0,100);
      PDF::gPoisTab.setRangeSigma(1,0,0);
      PDF::gPoisTab.setRangeX(61,0,60);
      PDF::gPoisTab.tabulate();
    }
    pole->initAnalysis();

    //
    // Now coverage
    //
    coverage->setDumpBase(dump.getValue().c_str());
    coverage->setFixedSig(doFixSig.getValue());
    coverage->setEffDist( effMin.getValue(), effSigma.getValue(), static_cast<PDF::DISTYPE>(effDist.getValue()));
    coverage->setBkgDist( bkgMin.getValue(), bkgSigma.getValue(), static_cast<PDF::DISTYPE>(bkgDist.getValue()));
    coverage->setEffBkgCorr(beCorr.getValue());
    coverage->checkEffBkgDists(); // will make sure the settings above are OK - it will update pole if changes are made
    //
    coverage->setTestHyp(hypTestMin.getValue(), hypTestMax.getValue(), hypTestStep.getValue());
    //
    coverage->setVerbose(verboseCov.getValue());
    coverage->setPole(pole);
    //
    coverage->collectStats(doStats.getValue());
    coverage->setNloops(nLoops.getValue());
    coverage->setSeed(rSeed.getValue()+rSeedOfs.getValue());
    //
    coverage->setSTrue(sMin.getValue(), sMax.getValue(), sStep.getValue());
    coverage->setEffTrue(effMin.getValue(), effMax.getValue(), effStep.getValue());
    coverage->setBkgTrue(bkgMin.getValue(), bkgMax.getValue(), bkgStep.getValue());
  }
  catch (ArgException e) {
    std::cout << "ERROR: " << e.error() << " " << e.argId() << std::endl;
  }
}
