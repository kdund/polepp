#include <iostream>
#include <fstream>
#include <iomanip>
#include <signal.h>
#include <tclap/CmdLine.h> // Command line parser
#include "Coverage.h"

Coverage gCoverage;
std::string gFilename="";

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

    ValueArg<int>    nLoops(   "","nloops",  "number of loops",    false,1000,"int");
    ValueArg<int>    rSeed(    "","rseed",   "rnd seed" ,          false,rndSeed,"int");
    ValueArg<int>    rSeedOfs( "","rseedofs","rnd seed offset" ,   false,0,"int");
    ValueArg<double> strue(    "","strue",   "s_true",         false,1.0,"float");
    ValueArg<double> effSigma( "","esigma",  "sigma of efficiency",false,0.2,"float");
    ValueArg<double> effTrue(  "","etrue",   "eff true",       false,1.0,"float");
    ValueArg<double> bkgSigma( "","bsigma",  "sigma of background",false,0.2,"float");
    ValueArg<double> bkgTrue(  "","btrue",   "bkg true",false,3.0,"float");
    ValueArg<double> beCorr(   "","corr",    "corr(bkg,eff)",false,0.5,"float");
    ValueArg<int>    effDist(  "","effdist", "Efficiency distribution",false,4,"int");
    ValueArg<int>    bkgDist(  "","bkgdist", "Background distribution",false,4,"int");
    ValueArg<std::string> dump("","dump",    "dump filename",false,"","string");
    //
    cmd.add(effSigma);
    cmd.add(effTrue);
    //
    cmd.add(bkgSigma);
    cmd.add(bkgTrue);
    //
    cmd.add(effDist);
    cmd.add(bkgDist);
    //
    cmd.add(beCorr);
    //
    cmd.add(strue);
    //
    cmd.add(nLoops);
    cmd.add(rSeedOfs);
    cmd.add(rSeed);
    cmd.add(dump);
    //
    cmd.parse(argc,argv);
    //
    gFilename = dump.getValue();
    gCoverage.setEffDist( effTrue.getValue(), effSigma.getValue(), static_cast<DISTYPE>(effDist.getValue()));
    gCoverage.setBkgDist( bkgTrue.getValue(), bkgSigma.getValue(), static_cast<DISTYPE>(bkgDist.getValue()));
    gCoverage.setEffBkgCorr(beCorr.getValue());
    gCoverage.checkEffBkgDists(); // will make sure the settings above are OK - it will update pole if changes are made
    //
    gCoverage.collectStats(true);
    gCoverage.setNloops(nLoops.getValue());
    gCoverage.setSeed(rSeed.getValue()+rSeedOfs.getValue());
    //
    gCoverage.setSTrue(strue.getValue(), strue.getValue(), 0.0);
    gCoverage.setEffTrue(effTrue.getValue(), effTrue.getValue(), 0.0);
    gCoverage.setBkgTrue(bkgTrue.getValue(), bkgTrue.getValue(), 0.0);
    //
    
  }
  catch (ArgException e) {
    cout << "ERROR: " << e.error() << " " << e.argId() << endl;
  }
}

int main(int argc, char *argv[]) {
  //
  processArgs(argc,argv);
  //
  gCoverage.printSetup();
  gCoverage.doExpTest();
  gCoverage.dumpExperiments(gFilename);
  //
}
