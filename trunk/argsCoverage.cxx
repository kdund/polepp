#include <iostream>
#include <fstream>
#include <iomanip>
#include <tclap/CmdLine.h> // Command line parser
#include "Coverage.h"


using namespace TCLAP;
void argsCoverage(Coverage *coverage, Pole *pole, int argc, char *argv[]) {
  std::cout << "Calling argsCoverage()" << std::endl;
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

    ValueArg<int>    effDist(   "","effdist", "Efficiency distribution",false,2,"int",cmd);
    ValueArg<double> effSigma(  "","effsigma",  "sigma of efficiency",false,0.2,"float",cmd);
    ValueArg<double> effMin(    "","effmin",    "min eff true",       false,1.0,"float",cmd);
    ValueArg<double> effMax(    "","effmax",    "max eff true",       false,1.0,"float",cmd);
    ValueArg<double> effStep(   "","effstep",   "step eff true",      false,1.0,"float",cmd);

    ValueArg<int>    bkgDist(   "","bkgdist", "Background distribution",false,0,"int",cmd);
    ValueArg<double> bkgSigma(  "","bkgsigma",  "sigma of background",false,0.2,"float",cmd);
    ValueArg<double> bkgMin(    "","bkgmin",    "min bkg true",false,0.0,"float",cmd);
    ValueArg<double> bkgMax(    "","bkgmax",    "max bkg true",false,0.0,"float",cmd);
    ValueArg<double> bkgStep(   "","bkgstep",   "step bkg true",false,1.0,"float",cmd);
    //
    ValueArg<double> beCorr(    "", "corr",   "corr(bkg,eff)",false,0.0,"float",cmd);
   //
    ValueArg<double> dMus(      "", "dmus",     "step size in findBestMu",false,0.002,"float",cmd);
    ValueArg<int>    nMus(      "", "nmus",     "maximum number of steps in findBestMu",false,100,"float",cmd);
    //
    ValueArg<double> threshBS("","threshbs",  "binary search (limit) threshold" ,false,0.001,"float",cmd);
    ValueArg<double> threshAlpha("","threshalpha",  "threshold for accepting a cl in % of (1-cl)" ,false,0.01,"float",cmd);
    //
//     ValueArg<double> hypTestMin( "","hmin",   "hypothesis test min" ,false,0.0,"float",cmd);
//     ValueArg<double> hypTestMax( "","hmax",   "hypothesis test max" ,false,35.0,"float",cmd);
//     ValueArg<double> hypTestStep("","hstep",  "hypothesis test step" ,false,0.01,"float",cmd);
    //
    ValueArg<double> effIntScale( "","effintscale","eff num sigma in integral", false,5.0,"float",cmd);
    ValueArg<int>    effIntN(     "","effintn",    "eff: N points in integral", false,21, "int",cmd);
    ValueArg<double> bkgIntScale( "","bkgintscale","bkg num sigma in integral", false,5.0,"float",cmd);
    ValueArg<int>    bkgIntN(     "","bkgintn",    "bkg: N points in integral", false,21, "int",cmd);
    //
    ValueArg<double> tabPoisMin( "","poismin",  "minimum mean value in table", false,0.0,"float",cmd);
    ValueArg<double> tabPoisMax( "","poismax",  "maximum mean value in table", false,100.0,"float",cmd);
    ValueArg<int>    tabPoisNM(  "","poisnm",   "number of mean value in table", false,100000,"int",cmd);
    ValueArg<int>    tabPoisNX(  "","poisnx",   "maximum value of N in table", false,200,"int",cmd);
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
    pole->setUseCoverage(true);
    pole->setCL(confLevel.getValue());
    //
    pole->setBestMuStep(dMus.getValue());
    pole->setBestMuNmax(nMus.getValue());
    //
    pole->setEffPdf( effMin.getValue(), effSigma.getValue(), static_cast<PDF::DISTYPE>(effDist.getValue()) );
    pole->setEffObs();
    pole->setBkgPdf( bkgMin.getValue(), bkgSigma.getValue(), static_cast<PDF::DISTYPE>(bkgDist.getValue()) );
    pole->setBkgObs();
    pole->checkEffBkgDists();
    pole->setEffPdfBkgCorr(beCorr.getValue());

    pole->setTrueSignal(sMin.getValue());

    pole->getMeasurement().setEffInt(effIntScale.getValue(),effIntN.getValue());
    pole->getMeasurement().setBkgInt(bkgIntScale.getValue(),bkgIntN.getValue());
    //
    pole->setBSThreshold(threshBS.getValue());
    pole->setAlphaThreshold(threshAlpha.getValue());
    //
    pole->setTestHyp(0.0,1.0,0.1);//hypTestMin.getValue(), hypTestMax.getValue(), hypTestStep.getValue());

    pole->setMinMuProb(minProb.getValue());
    //
    if (!noTabulated.getValue()) {
      PDF::gPoisTab.setRangeMean( tabPoisNM.getValue(), tabPoisMin.getValue(), tabPoisMax.getValue() );
      PDF::gPoisTab.setRangeX(tabPoisNX.getValue()+1,0,tabPoisNX.getValue());
      PDF::gPoisTab.tabulate();
    }
    pole->initAnalysis();

    //
    // Now coverage
    //
    coverage->setPole(pole);
    coverage->setDumpBase(dump.getValue().c_str());
    coverage->setVerbose(verboseCov.getValue());
    //
    coverage->collectStats(doStats.getValue());
    coverage->setNloops(nLoops.getValue());
    coverage->setSeed(rSeed.getValue()+rSeedOfs.getValue());
    //
    coverage->setFixedSig(doFixSig.getValue());
    coverage->setSTrue(sMin.getValue(), sMax.getValue(), sStep.getValue());
    coverage->setEffTrue(effMin.getValue(), effMax.getValue(), effStep.getValue());
    coverage->setBkgTrue(bkgMin.getValue(), bkgMax.getValue(), bkgStep.getValue());
  }
  catch (ArgException e) {
    std::cout << "ERROR: " << e.error() << " " << e.argId() << std::endl;
  }
}
