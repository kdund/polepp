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
    //
    ValueArg<double> effSigma(  "", "esigma","sigma of efficiency",false,0.2,"float");
    ValueArg<double> effMeas(   "", "emeas",  "measured efficiency",false,1.0,"float");
    SwitchArg        effNoDist( "E","ndeff", "No efficiency distribution",false);
    //
    ValueArg<double> bkgSigma(  "", "bsigma","sigma of background",false,0.2,"float");
    ValueArg<double> bkgMeas(   "", "bmeas",   "measured background",false,0.2,"float");
    SwitchArg        bkgNoDist( "B","ndbkg",   "No background distribution",false);
    //
    ValueArg<double> dMus(      "","dmus",    "step size in findBestMu",false,0.002,"float");
    ValueArg<int>    belt(   "","belt", "maximum n for findBestMu" ,false,50,"int");
    //
    ValueArg<double> hypTestMin( "","hmin",   "hypothesis test min" ,false,0.0,"float");
    ValueArg<double> hypTestMax( "","hmax",   "hypothesis test max" ,false,35.0,"float");
    ValueArg<double> hypTestStep("","hstep",  "hypothesis test step" ,false,0.01,"float");
    //
    ValueArg<double> effIntMin( "","emin",    "eff min in integral",  false,0.0,"float");
    ValueArg<double> effIntMax( "","emax",    "eff max in integral",  false,0.0,"float");
    ValueArg<double> effIntStep("","estep",   "eff step in integral", false,-1.0,"float");
    ValueArg<double> bkgIntMin( "","bmin",    "bkg min in integral",  false,0.0,"float");
    ValueArg<double> bkgIntMax( "","bmax",    "bkg max in integral",  false,0.0,"float");
    ValueArg<double> bkgIntStep("","bstep",   "bkg step in integral", false,-1.0,"float");

    ValueArg<int>    doVerbose(   "V","verbose", "verbose pole",    false,0,"int");
    //
    cmd.add(doVerbose);
    //
    cmd.add(hypTestMin);
    cmd.add(hypTestMax);
    cmd.add(hypTestStep);

    cmd.add(effIntMin);
    cmd.add(effIntMax);
    cmd.add(effIntStep);
    cmd.add(bkgIntMin);
    cmd.add(bkgIntMax);
    cmd.add(bkgIntStep);

    cmd.add(belt);
    cmd.add(dMus);

    cmd.add(effSigma);
    cmd.add(effMeas);
    cmd.add(effNoDist);

    cmd.add(bkgSigma);
    cmd.add(bkgMeas);
    cmd.add(bkgNoDist);

    cmd.add(sTrue);
    cmd.add(coverage);
    cmd.add(confLevel);
    cmd.add(nObs);

    //
    cmd.parse(argc,argv);
    //
    pole->setCL(confLevel.getValue());
    pole->setNobserved(nObs.getValue());
    //
    pole->setEffDist( effMeas.getValue(), effSigma.getValue(), effNoDist.getValue() );
    pole->setBkgDist( bkgMeas.getValue(), bkgSigma.getValue(), bkgNoDist.getValue() );

    pole->setTrueSignal(sTrue.getValue());
    pole->setCoverage(coverage.getValue());

    pole->setDmus(dMus.getValue());
    pole->setEffInt(effIntMin.getValue(),effIntMax.getValue(),effIntStep.getValue());
    pole->setBkgInt(bkgIntMin.getValue(),bkgIntMax.getValue(),bkgIntStep.getValue());
    //
    pole->setBelt(belt.getValue()); // call after nObserved is set.
    pole->setBeltMax(belt.getValue()*2); // maximum allocated
    pole->setTestHyp(hypTestMin.getValue(), hypTestMax.getValue(), hypTestStep.getValue());
    //
    //    pole->initPoisson(50000,60,200);
    //    pole->initGauss(10000,10.0);
    pole->initIntArrays();
    pole->initBeltArrays();
    //
    pole->setVerbose(doVerbose.getValue());

  }
  catch (ArgException e) {
    cout << "ERROR: " << e.error() << " " << e.argId() << endl;
  }
}

int main(int argc, char *argv[]) {

  Pole pole;
  //
  processArgs(&pole, argc, argv);
  //
  if (pole.checkParams()) {
    pole.printSetup();
    pole.analyseExperiment();
    pole.printLimit();
    pole.setNobserved(4);
    pole.analyseExperiment();
    pole.printLimit();
  }
  //
}
