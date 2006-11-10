#include <iostream>
#include <fstream>
#include <iomanip>
#include <signal.h>
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

    ValueArg<int>    nLoops(   "","nloops",  "number of loops",    false,1000,"int",cmd);
    ValueArg<int>    rSeed(    "","rseed",   "rnd seed" ,          false,rndSeed,"int",cmd);
    ValueArg<int>    rSeedOfs( "","rseedofs","rnd seed offset" ,   false,0,"int",cmd);
    ValueArg<double> strue(    "","strue",   "s_true",         false,1.6,"float",cmd);
    ValueArg<double> effSigma( "","effsigma",  "sigma of efficiency",false,0.2,"float",cmd);
    ValueArg<double> effMean(  "","effmean",   "eff true mean",       false,1.0,"float",cmd);
    ValueArg<double> effScale( "","effscale",  "effscale factor",false,1.0,"float",cmd);
    ValueArg<int>    effDist(  "","effdist", "Efficiency distribution",false,2,"int",cmd);
    ValueArg<double> bkgSigma( "","bkgsigma",  "sigma of background",false,0.2,"float",cmd);
    ValueArg<double> bkgMean(  "","bkgmean",   "bkg true mean",false,2.0,"float",cmd);
    ValueArg<double> bkgScale(  "","bkgscale",  "bkgscale factor",false,1.0,"float",cmd);
    ValueArg<int>    bkgDist(  "","bkgdist", "Background distribution",false,2,"int",cmd);
    ValueArg<double> beCorr(   "","corr",    "corr(bkg,eff)",false,0.0,"float",cmd);
    ValueArg<std::string> dump("","dump",    "dump filename",false,"","string",cmd);
    //
    cmd.parse(argc,argv);
    //
    // First set Pole
    //
    gPole.setPoisson(&PDF::gPoisTab);
    gPole.setGauss(&PDF::gGauss);
    gPole.setGauss2D(&PDF::gGauss2D);
    gPole.setLogNormal(&PDF::gLogNormal);

    gPole.setEffPdfScale( effScale.getValue() );
    gPole.setEffPdf( effMean.getValue(), effSigma.getValue(), static_cast<PDF::DISTYPE>(effDist.getValue()) );
    gPole.setEffObs();
    gPole.setBkgPdfScale( bkgScale.getValue() );
    gPole.setBkgPdf( bkgMean.getValue(), bkgSigma.getValue(), static_cast<PDF::DISTYPE>(bkgDist.getValue()) );
    gPole.setBkgObs();
    gPole.checkEffBkgDists();
    gPole.setEffBkgPdfCorr(beCorr.getValue());

    gPole.setTrueSignal(strue.getValue());
    //
    gCoverage.setPole(&gPole);
    gCoverage.setDumpBase(dump.getValue().c_str());
    //
    gCoverage.collectStats(true);
    gCoverage.setNloops(nLoops.getValue());
    gCoverage.setSeed(rSeed.getValue()+rSeedOfs.getValue());
    //
    gCoverage.setSTrue(strue.getValue(), strue.getValue(), 0.0);
    gCoverage.setEffTrue(effMean.getValue(), effMean.getValue(), 0.0);
    gCoverage.setBkgTrue(bkgMean.getValue(), bkgMean.getValue(), 0.0);
    //
  }
  catch (ArgException e) {
    std::cout << "ERROR: " << e.error() << " " << e.argId() << std::endl;
  }
}

int main(int argc, char *argv[]) {
  //
  processArgs(argc,argv);
  //
  gCoverage.printSetup();
  gCoverage.doExpTest();
  gCoverage.dumpExperiments(false);
  //
}
