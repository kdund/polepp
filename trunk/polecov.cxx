#include <iostream>
#include <fstream>
#include <iomanip>
#include <tclap/CmdLine.h> // Command line parser
#include "Coverage.h"

Coverage gCoverage;
Pole     gPole;

using namespace TCLAP;
void processArgs(int argc, char *argv[]) {
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

    ValueArg<int>    nLoops(    "","nloops",  "number of loops",false,1,"int");
    ValueArg<double> confLevel( "","cl",      "confidence level",false,0.9,"float");
    ValueArg<int>    rSeed(     "","rseed",   "rnd seed" ,false,rndSeed,"int");
    ValueArg<int>    rSeedOfs(  "","rseedofs","rnd seed offset" ,false,0,"int");
    ValueArg<int>    nBelt(      "","nbelt",  "n Belt" ,false,-1,"int");
    ValueArg<double> sMin(      "","smin",    "min s_true",false,1.0,"float");
    ValueArg<double> sMax(      "","smax",    "max s_true",false,1.0,"float");
    ValueArg<double> sStep(     "","sstep",    "step s_true",false,1.0,"float");
    ValueArg<double> effSigma(  "","esigma","sigma of efficiency",false,0.2,"float");
    ValueArg<double> effMin(    "","emin",    "min eff true",false,1.0,"float");
    ValueArg<double> effMax(    "","emax",    "max eff true",false,1.0,"float");
    ValueArg<double> effStep(   "","estep",   "step eff true",false,1.0,"float");
    ValueArg<double> bkgSigma(  "","bsigma","sigma of background",false,0.2,"float");
    ValueArg<double> bkgMin(    "","bmin",    "min bkg true",false,0.0,"float");
    ValueArg<double> bkgMax(    "","bmax",    "max bkg true",false,0.0,"float");
    ValueArg<double> bkgStep(   "","bstep",   "step bkg true",false,1.0,"float");
    ValueArg<double> beCorr(    "","corr",    "corr(bkg,eff)",false,1000.0,"float");
    ValueArg<double> dMus(      "","dmus",    "step size in findBestMu",false,0.01,"float");
    ValueArg<double> muTestMin( "","mumin",   "test min" ,false,0.0,"float");
    ValueArg<double> muTestMax( "","mumax",   "test max" ,false,35.0,"float");
    ValueArg<double> muTestStep("","mustep",  "test step" ,false,0.02,"float");
    ValueArg<double> effIntMin( "","eintmin",    "eff min in integral",  false,0.0,"float");
    ValueArg<double> effIntMax( "","eintmax",    "eff max in integral",  false,2.0,"float");
    ValueArg<double> effIntStep("","eintstep",   "eff step in integral", false,-1.0,"float");
    ValueArg<double> bkgIntMin( "","bintmin",    "bkg min in integral",  false,0.0,"float");
    ValueArg<double> bkgIntMax( "","bintmax",    "bkg max in integral",  false,2.0,"float");
    ValueArg<double> bkgIntStep("","bintstep",   "bkg step in integral", false,-1.0,"float");
    //
    SwitchArg        doStats("C","stats", "Collect statistics",false);
    SwitchArg        doLogNorm("L","lognorm", "Use lognormal instead of gauss",false);
    SwitchArg        doFixEff("E","fixeff", "Fixed meas. efficiency",false);
    SwitchArg        doFixBkg("B","fixbkg", "Fixed meas. background",false);
    SwitchArg        doFixSig("S","fixsig", "Fixed meas. n_observed",false);
    SwitchArg        doNoDistEff("e","ndeff",  "No efficiency distribution",false);
    SwitchArg        doNoDistBkg("b","ndbkg",  "No background distribution",false);
    //
    SwitchArg        doExamples("","example", "Print examples",false);
    //
    ValueArg<int>    verboseCov(   "V","verbcov", "verbose coverage",false,0,"int");
    ValueArg<int>    verbosePol(   "P","verbpol", "verbose pole",    false,0,"int");
    //
    cmd.add(verboseCov);
    cmd.add(verbosePol);
    //
    cmd.add(effIntMin);
    cmd.add(effIntMax);
    cmd.add(effIntStep);
    cmd.add(bkgIntMin);
    cmd.add(bkgIntMax);
    cmd.add(bkgIntStep);
    //
    cmd.add(doLogNorm);
    cmd.add(effSigma);
    cmd.add(effMin);
    cmd.add(effMax);
    cmd.add(effStep);
    //
    cmd.add(bkgSigma);
    cmd.add(bkgMin);
    cmd.add(bkgMax);
    cmd.add(bkgStep);
    //
    cmd.add(beCorr);
    //
    cmd.add(sMin);
    cmd.add(sMax);
    cmd.add(sStep);
    //
    cmd.add(nLoops);
    cmd.add(rSeedOfs);
    cmd.add(rSeed);
    cmd.add(nBelt);
    cmd.add(dMus);
    cmd.add(muTestMin);
    cmd.add(muTestMax);
    cmd.add(muTestStep);
    cmd.add(doStats);
    cmd.add(doFixEff);
    cmd.add(doFixBkg);
    cmd.add(doFixSig);
    cmd.add(doNoDistEff);
    cmd.add(doNoDistBkg);
    cmd.add(confLevel);
    //    cmd.add(doExamples);
    //
    cmd.parse(argc,argv);
    //
    if (doExamples.getValue()) {
    }
    //
    gPole.useLogNormal(doLogNorm.getValue());
    gPole.setVerbose(verbosePol.getValue());
    gPole.setNobserved(0);
    gPole.setCoverage(true);
    gPole.setCL(confLevel.getValue());
    gPole.setDmus(dMus.getValue());
    gPole.setEffDist( effMin.getValue(), effSigma.getValue(), doNoDistEff.getValue());
    gPole.setBkgDist( bkgMin.getValue(), bkgSigma.getValue(), doNoDistBkg.getValue());
    gPole.setEffInt(effIntMin.getValue(),effIntMax.getValue(),effIntStep.getValue());
    gPole.setBkgInt(bkgIntMin.getValue(),bkgIntMax.getValue(),bkgIntStep.getValue());
    //
    gPole.setBelt(nBelt.getValue());
    gPole.setBeltMax(nBelt.getValue()*2);
    gPole.setTestHyp(muTestMin.getValue(), muTestMax.getValue(), muTestStep.getValue());
    gPole.printSetup();
    gPole.initPoisson(50000,60,200);
    gPole.initGauss(10000,10.0);
    gPole.initIntArrays();
    gPole.initBeltArrays();

    //
    gCoverage.useLogNormal(doLogNorm.getValue());
    gCoverage.setVerbose(verboseCov.getValue());
    gCoverage.setPole(&gPole);
    //
    gCoverage.collectStats(doStats.getValue());
    gCoverage.setNloops(nLoops.getValue());
    gCoverage.setSeed(rSeed.getValue()+rSeedOfs.getValue());
    //
    gCoverage.setSTrue(sMin.getValue(), sMax.getValue(), sStep.getValue());
    gCoverage.setEffTrue(effMin.getValue(), effMax.getValue(), effStep.getValue());
    gCoverage.setBkgTrue(bkgMin.getValue(), bkgMax.getValue(), bkgStep.getValue());
    //
    gCoverage.setFixedEff(doFixEff.getValue());
    gCoverage.setFixedBkg(doFixBkg.getValue());
    gCoverage.setFixedSig(doFixSig.getValue());
    double c = beCorr.getValue(); // Correlation must be set after
    bool uc = true;
    if (c>999.0) uc=false;
    gCoverage.setCorr(c);
    gCoverage.setUseCorr(uc);
    //
    //    gCoverage.initTabs();
    //
    
  }
  catch (ArgException e) {
    cout << "ERROR: " << e.error() << " " << e.argId() << endl;
  }
}

int main(int argc, char *argv[]) {
  
  processArgs(argc,argv);
  //
  if (gPole.checkParams()) {
    gCoverage.printSetup();
    gCoverage.doLoop();
  }
  //
}
