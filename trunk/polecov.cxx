#include <iostream>
#include <fstream>
#include <iomanip>
#include <signal.h>
#include <tclap/CmdLine.h> // Command line parser
#include "Coverage.h"

Coverage gCoverage;
Pole     gPole;
std::string gDumpFile;

void time_stamp(std::string & stamp) {
  time_t epoch;
  time(&epoch);
  struct tm *time;
  char tst[32];
  time = localtime(&epoch); // time_t == long int
  strftime(tst,32,"%d/%m/%Y %H:%M:%S",time);
  stamp = tst;
}

void my_sighandler(int a) {
  std::string timestamp;
  time_stamp(timestamp);
  //
  if (a==SIGUSR1) {
    gCoverage.calcCoverage();
    std::string header("STATUS ( ");
    header += timestamp;
    header += " ) : ";
    gCoverage.outputCoverageResult(1);
  } else {
    std::cout << "WARNING (" << timestamp << " ) Job aborting (signal = " << a
	      << " ). Will output data from unfinnished loop.\n" << std::endl;
    gCoverage.calcCoverage();
    gCoverage.outputCoverageResult(0);
    exit(-1);
  }
}

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

    ValueArg<int>    nLoops(    "","nloops",  "number of loops",    false,1,"int");
    ValueArg<double> confLevel( "","cl",      "confidence level",   false,0.9,"float");
    ValueArg<int>    rSeed(     "","rseed",   "rnd seed" ,          false,rndSeed,"int");
    ValueArg<int>    rSeedOfs(  "","rseedofs","rnd seed offset" ,   false,0,"int");
    ValueArg<int>    nBelt(     "","nbelt",   "n Belt" ,            false,-1,"int");
    ValueArg<double> sMin(      "","smin",    "min s_true",         false,1.0,"float");
    ValueArg<double> sMax(      "","smax",    "max s_true",         false,1.0,"float");
    ValueArg<double> sStep(     "","sstep",   "step s_true",        false,1.0,"float");
    ValueArg<double> effSigma(  "","esigma",  "sigma of efficiency",false,0.2,"float");
    ValueArg<double> effMin(    "","emin",    "min eff true",       false,1.0,"float");
    ValueArg<double> effMax(    "","emax",    "max eff true",       false,1.0,"float");
    ValueArg<double> effStep(   "","estep",   "step eff true",      false,1.0,"float");
    ValueArg<double> bkgSigma(  "","bsigma",  "sigma of background",false,0.2,"float");
    ValueArg<double> bkgMin(    "","bmin",    "min bkg true",false,0.0,"float");
    ValueArg<double> bkgMax(    "","bmax",    "max bkg true",false,0.0,"float");
    ValueArg<double> bkgStep(   "","bstep",   "step bkg true",false,1.0,"float");
    ValueArg<double> beCorr(    "","corr",    "corr(bkg,eff)",false,0.0,"float");
    ValueArg<double> dMus(      "","dmus",    "step size in findBestMu",false,0.01,"float");
    ValueArg<double> muTestMin( "","mumin",   "test min" ,false,0.0,"float");
    ValueArg<double> muTestMax( "","mumax",   "test max" ,false,35.0,"float");
    ValueArg<double> muTestStep("","mustep",  "test step" ,false,0.02,"float");
    //    ValueArg<double> effIntMin( "","eintmin",    "eff min in integral",  false,0.0,"float");
    //    ValueArg<double> effIntMax( "","eintmax",    "eff max in integral",  false,2.0,"float");
    ValueArg<double> effIntScale( "","eintscale",    "n sigmas in integral (eff)",  false,5.0,"float");
    ValueArg<double> effIntStep("","eintstep",   "eff step in integral", false,-1.0,"float");
    //    ValueArg<double> bkgIntMin( "","bintmin",    "bkg min in integral",  false,0.0,"float");
    //    ValueArg<double> bkgIntMax( "","bintmax",    "bkg max in integral",  false,2.0,"float");
    ValueArg<double> bkgIntScale( "","bintscale",    "n sigmas in integral (bkg)",  false,5.0,"float");
    ValueArg<double> bkgIntStep("","bintstep",   "bkg step in integral", false,-1.0,"float");

    ValueArg<int>    effDist("","effdist",  "Efficiency distribution",false,1,"int");
    ValueArg<int>    bkgDist("","bkgdist",  "Background distribution",false,0,"int");
    //
    SwitchArg        doStats("C","stats", "Collect statistics",false);
    SwitchArg        doFixSig("S","fixsig", "Fixed meas. n_observed",false);
    SwitchArg        doNLR("N","nlr", "Use NLR",false);
    SwitchArg        useTabulated("K","notab","Do not use tabulated poisson",true);
    //
    ValueArg<std::string> dump("","dump",    "dump filename",false,"","string");
    //
    SwitchArg        doExamples("","example", "Print examples",false);
    //
    ValueArg<int>    verboseCov(   "V","verbcov", "verbose coverage",false,0,"int");
    ValueArg<int>    verbosePol(   "P","verbpol", "verbose pole",    false,0,"int");
    //
    cmd.add(useTabulated);
    cmd.add(verboseCov);
    cmd.add(verbosePol);
    cmd.add(dump);
    //
    cmd.add(doNLR);
    //
    //    cmd.add(effIntMin);
    //    cmd.add(effIntMax);
    cmd.add(effIntScale);
    cmd.add(effIntStep);
    //    cmd.add(bkgIntMin);
    //    cmd.add(bkgIntMax);
    cmd.add(bkgIntScale);
    cmd.add(bkgIntStep);
    //
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
    cmd.add(doFixSig);
    cmd.add(effDist);
    cmd.add(bkgDist);
    cmd.add(confLevel);
    //    cmd.add(doExamples);
    //
    cmd.parse(argc,argv);
    //
    if (doExamples.getValue()) {
    }
    //
    gPole.setNLR(doNLR.getValue());
    gPole.setVerbose(verbosePol.getValue());
    gPole.setNobserved(0);
    gPole.setCoverage(true);
    gPole.setCL(confLevel.getValue());
    gPole.setDmus(dMus.getValue());
    gPole.setEffMeas( effMin.getValue(), effSigma.getValue(), static_cast<DISTYPE>(effDist.getValue()));
    gPole.setBkgMeas( bkgMin.getValue(), bkgSigma.getValue(), static_cast<DISTYPE>(bkgDist.getValue()));
    gPole.checkEffBkgDists(); // will make sure the settings above are OK - it will update pole if changes are made
    gPole.setEffBkgCorr(beCorr.getValue());
    gPole.setEffInt(effIntScale.getValue(),effIntStep.getValue());
    gPole.setBkgInt(bkgIntScale.getValue(),bkgIntStep.getValue());
    //
    gPole.setBelt(nBelt.getValue());
    gPole.setBeltMax(nBelt.getValue()*2);
    gPole.setTestHyp(muTestMin.getValue(), muTestMax.getValue(), muTestStep.getValue());
    gPole.printSetup();
    if (useTabulated.getValue()) {
      //      gPole.initPoisson(100000,100,50);
      //      gPole.initPoisson(500000,100,50);
      //      gPole.initGauss(1000000,10.0);
      gPole.initPoisson(50000,100,50);
      //      gPole.initPoisson(50000,100,200);
      gPole.initGauss(50000,10.0);
    }
    gPole.initIntArrays();
    gPole.initBeltArrays();

    //
    gCoverage.setDumpBase(dump.getValue().c_str());
    gCoverage.setFixedSig(doFixSig.getValue());
    gCoverage.setEffDist( effMin.getValue(), effSigma.getValue(), static_cast<DISTYPE>(effDist.getValue()));
    gCoverage.setBkgDist( bkgMin.getValue(), bkgSigma.getValue(), static_cast<DISTYPE>(bkgDist.getValue()));
    gCoverage.setEffBkgCorr(beCorr.getValue());
    gCoverage.checkEffBkgDists(); // will make sure the settings above are OK - it will update pole if changes are made
    //
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
    //
    //    gCoverage.initTabs();
    //
    
  }
  catch (ArgException e) {
    std::cout << "ERROR: " << e.error() << " " << e.argId() << std::endl;
  }
}

int main(int argc, char *argv[]) {
  // Trap LSF specific signals
  //                                value  bkill  memlimit  runlimit  cpulimit  filelimit  job_starter
  //                               ====================================================================
  signal(SIGINT, my_sighandler); //   2      2nd      1st       -         -         -        failure
  // SIGKILL not trapable        //   9      3rd      3rd       -         -         -           -
  signal(SIGUSR2,my_sighandler); //  12       -        -     reached      -         -           -
  signal(SIGTERM,my_sighandler); //  15      1st      2nd       -         -         -           -
  signal(SIGXCPU,my_sighandler); //  24       -        -        -      reached      -           -
  signal(SIGXFSZ,my_sighandler); //  25       -        -        -         -      reached        -
  //                               ====================================================================
  // General signals
  signal(SIGSEGV,my_sighandler); // Segmentation fault
  signal(SIGUSR1,my_sighandler);
  signal(SIGIO,  my_sighandler); // Directory access error
  //
  processArgs(argc,argv);
  //
  if (gPole.checkParams()) {
    gCoverage.printSetup();
    gCoverage.doLoop();
  }
}
