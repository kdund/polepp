#include <iostream>
#include <fstream>
#include <iomanip>
#include <signal.h>
#include <tclap/CmdLine.h> // Command line parser
#include "Pole.h"
#include "Power.h"

Power gPower;

using namespace TCLAP;

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
    gPower.calcPower();
    std::string header("STATUS ( ");
    header += timestamp;
    header += " ) : ";
    gPower.outputResult(header.c_str());
  } else {
    std::cout << "WARNING (" << timestamp << " ) Job aborting (signal = " << a
	      << " ). Will output data from unfinnished loop.\n" << std::endl;
    gPower.calcPower();
    gPower.outputResult(0); // always output
    exit(-1);
  }
}

void processArgs(Pole *pole, int argc, char *argv[]) {
  //
  try {
    // Create a CmdLine object.
    // Arg1 = string printed at the end whenever --help is used or an error occurs
    // Arg2 = delimiter character between opt and its value
    // Arg3 = version number given when --version is used
    CmdLine cmd("Try again, friend.", ' ', "0.99");

    ValueArg<double> minProb( "","minp",       "minimum probability",false,-1.0,"float");
    ValueArg<int>    nLoops(    "","nloops",  "number of loops",    false,100,"int");

    ValueArg<int>    nObs(      "","nobs",     "number observed events",false,1,"int");
    ValueArg<double> confLevel( "","cl",       "confidence level",false,0.9,"float");
    ValueArg<double> sTrue(     "","strue",   "s_true = s0",false,1.0,"float");
    ValueArg<double> s1min(     "","s1min",   "s1_min",false,1.0,"float");
    ValueArg<double> s1max(     "","s1max",   "s1_max",false,1.0,"float");
    ValueArg<double> s1step(     "","s1step",   "s1_step",false,0.1,"float");


    SwitchArg        useTabulated("T","tab","Use tabulated poisson",false);
    //
    ValueArg<double> effSigma(  "", "esigma","sigma of efficiency",false,0.2,"float");
    ValueArg<double> effMeas(   "", "emeas",  "measured efficiency",false,1.0,"float");
    ValueArg<int>    effDist(   "","effdist",  "Efficiency distribution",false,1,"int");
    ValueArg<int>    method(    "m","method",     "method (1 - FHC2 (def), 2 - MBT)",false,1,"int");
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
    ValueArg<double> hypTestStep("","hstep",  "hypothesis test step" ,false,0.1,"float");
    //
    ValueArg<double> effIntScale( "","escale","eff n sigma in integral", false,5.0,"float");
    ValueArg<int>    effIntN("","en",   "eff: N points in integral",    false,21,"int");
    ValueArg<double> bkgIntScale( "","bscale","bkg n sigma in integral", false,5.0,"float");
    ValueArg<int>    bkgIntN("","bn",   "bkg: N points in integral",    false,21,"int");

    ValueArg<int>    doVerbPole(  "V","vpole",   "verbose pole",     false,0,"int");
    ValueArg<int>    doVerbPow(   "P","vpower",  "verbose power",    false,0,"int");
    //
    cmd.add(doVerbPole);
    cmd.add(doVerbPow);
    cmd.add(minProb);
    cmd.add(method);
    cmd.add(useTabulated);
    //
    cmd.add(hypTestMin);
    cmd.add(hypTestMax);
    cmd.add(hypTestStep);

    cmd.add(effIntScale);
    cmd.add(effIntN);
    cmd.add(bkgIntScale);
    cmd.add(bkgIntN);

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
    cmd.add(s1min);
    cmd.add(s1max);
    cmd.add(s1step);

    cmd.add(confLevel);
    cmd.add(nObs);
    cmd.add(nLoops);

    //
    cmd.parse(argc,argv);
    //
    pole->setPoisson(&PDF::gPoisson);
    pole->setGauss(&PDF::gGauss);
    pole->setMethod(method.getValue());
    pole->setCL(confLevel.getValue());
    pole->setNObserved(nObs.getValue());
    //
    pole->setEffMeas( effMeas.getValue(), effSigma.getValue(), static_cast<DISTYPE>(effDist.getValue()) );
    pole->setBkgMeas( bkgMeas.getValue(), bkgSigma.getValue(), static_cast<DISTYPE>(bkgDist.getValue()) );
    pole->checkEffBkgDists();
    pole->setEffBkgCorr(beCorr.getValue());

    pole->setTrueSignal(sTrue.getValue()); // s0
    
    pole->setCoverage(false);

    pole->setDmus(dMus.getValue());
    pole->setEffInt(effIntScale.getValue(),effIntN.getValue());
    pole->setBkgInt(bkgIntScale.getValue(),bkgIntN.getValue());
    //
    pole->setBelt(belt.getValue());
    pole->setTestHyp(hypTestMin.getValue(), hypTestMax.getValue(), hypTestStep.getValue());
    //
    if (useTabulated.getValue()) {
      PDF::gPoisson.init(100000,200,100);
      PDF::gGauss.init(0,10.0);
    }
    pole->initIntArrays();
    pole->initBeltArrays();
    pole->setMinMuProb(minProb.getValue());
    //
    pole->setVerbose(doVerbPole.getValue());
    //
    gPower.setPole(pole);
    gPower.setHypSignal(sTrue.getValue());
    gPower.setTrueSignal(s1min.getValue(),s1max.getValue(),s1step.getValue());
    gPower.setVerbose(doVerbPow.getValue());
    gPower.setLoops(nLoops.getValue());

  }
  catch (ArgException e) {
    std::cout << "ERROR: " << e.error() << " " << e.argId() << std::endl;
  }
}

int main(int argc, char *argv[]) {

  Pole pole;
  //
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

  processArgs(&pole, argc, argv);
  //
  //  if (pole.checkParams()) {
  pole.printSetup();

  //
  gPower.doLoop();

//   pole.initAnalysis();
//   if (pole.usesFHC2()) {
//     pole.findAllBestMu();
//   }
//   pole.calcPower();
}
