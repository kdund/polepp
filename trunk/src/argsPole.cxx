#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <tclap/CmdLine.h> // Command line parser
#include "Pole.h"


using namespace TCLAP;
void argsPole(LIMITS::Pole *pole, int argc, char *argv[]) {
  //
  if (pole==0) return;

  pole->initDefault();
  //
  try {
    // Create a CmdLine object.
    // Arg1 = string printed at the end whenever --help is used or an error occurs
    // Arg2 = delimiter character between opt and its value
    // Arg3 = version number given when --version is used
    std::string postMsg;
    std::ostringstream sstr;
    int dist = static_cast<int>(PDF::DIST_CONST);
    sstr << "Available distributions:" << std::endl;
    while (dist !=  static_cast<int>(PDF::DIST_LAST)) {
      sstr << " " << dist << std::setw(4) << " " << PDF::distTypeStr(static_cast<PDF::DISTYPE>(dist)) << std::endl;
      dist++;
    }
    postMsg = sstr.str();
    CmdLine cmd(postMsg.c_str(), ' ', "0.99");

    ValueArg<int>    nObs(      "","nobs",     "number observed events",false,1,"int",cmd);
    ValueArg<double> confLevel( "","cl",       "confidence level",false,0.9,"float",cmd);
    ValueArg<double> sTrue(     "","strue",   "s_true, only used if -C is active",false,1.0,"float",cmd);
    ValueArg<int>    method(    "m","method",     "method (1 - FHC2 (def), 2 - MBT)",false,1,"int",cmd);
    //
    ValueArg<double> minProb( "","minp",       "minimum probability",false,-1.0,"float",cmd);
    //
    ValueArg<double> effSigma(  "", "effsigma","sigma of efficiency",false,0.2,"float",cmd);
    ValueArg<double> effMeas(   "", "effmeas",  "measured efficiency",false,1.0,"float",cmd);
    ValueArg<int>    effDist(   "", "effdist",  "efficiency distribution",false,2,"int",cmd);
    ValueArg<double> effScale(  "", "effscale", "effscale factor",false,1.0,"float",cmd);
    //
    ValueArg<double> bkgSigma(  "", "bkgsigma", "sigma of background",false,0.0,"float",cmd);
    ValueArg<double> bkgMeas(   "", "bkgmeas",  "measured background",false,0.0,"float",cmd);
    ValueArg<int>    bkgDist(   "", "bkgdist",  "background distribution",false,0,"int",cmd);
    ValueArg<double> bkgScale(  "", "bkgscale", "bkgscale factor",false,1.0,"float",cmd);
    //
    ValueArg<double> beCorr(    "", "corr",     "corr(bkg,eff)",false,0.0,"float",cmd);
   //
    ValueArg<double> dMus(      "", "dmus",     "step size in findBestMu",false,0.01,"float",cmd);
    ValueArg<int>    nMus(      "", "nmus",     "maximum number of steps in findBestMu",false,20,"float",cmd);
    //
    ValueArg<double> threshBS(   "","threshbs",     "binary search (limit) threshold" ,false,0.0001,"float",cmd);
    ValueArg<double> threshPrec( "","threshprec",  "threshold for accepting a cl" ,false,0.0001,"float",cmd);
    //
    ValueArg<double> hypTestMin( "","hmin",   "hypothesis test min" ,false,0.0,"float",cmd);
    ValueArg<double> hypTestMax( "","hmax",   "hypothesis test max" ,false,35.0,"float",cmd);
    ValueArg<double> hypTestStep("","hstep",  "hypothesis test step" ,false,0.01,"float",cmd);
    //
    ValueArg<double> effIntNSigma(  "","effintnsigma","eff num sigma in integral", false,5.0,"float",cmd);
    ValueArg<double> bkgIntNSigma(  "","bkgintnsigma","bkg num sigma in integral", false,5.0,"float",cmd);
    ValueArg<int>    gslIntNCalls(  "","gslintncalls","number of calls in GSL integrator", false,100,"int",cmd);
    //
    ValueArg<double> tabPoleSMin(   "","tabpolesmin", "Pole table: minimum signal", false,0.0,"float",cmd);
    ValueArg<double> tabPoleSMax(   "","tabpolesmax", "Pole table: maximum signal", false,10.0,"float",cmd);
    ValueArg<int>    tabPoleSNStep( "","tabpolesn",   "Pole table: number of signals", false,100,"int",cmd);
    ValueArg<int>    tabPoleNMin(   "","tabpolenmin", "Pole table: minimum N(obs)", false, 0,"int",cmd);
    ValueArg<int>    tabPoleNMax(   "","tabpolenmax", "Pole table: maximum N(obs)", false,10,"int",cmd);
    SwitchArg        tabPole(       "I","poletab",    "Pole table: tabulated",false);
    cmd.add(tabPole);

    ValueArg<double> tabPoisMuMin( "","poismumin", "Poisson table: minimum mean value", false,0.0,"float",cmd);
    ValueArg<double> tabPoisMuMax( "","poismumax", "Poisson table: maximum mean value", false,50.0,"float",cmd);
    ValueArg<int>    tabPoisMuN(   "","poismun",   "Poisson table: number of mean value", false,100,"int",cmd);
    ValueArg<int>    tabPoisNmin(  "","poisnmin",  "Poisson table: minimum value of N", false,0,"int",cmd);
    ValueArg<int>    tabPoisNmax(  "","poisnmax",  "Poisson table: maximum value of N", false,35,"int",cmd);
    SwitchArg        tabPois(     "P","poistab",   "Poisson table: tabulated",false);
    cmd.add(tabPois);

    ValueArg<int>    doVerbose(   "V","verbose", "verbose pole",    false,0,"int",cmd);

    ValueArg<std::string> inputFile( "f" ,"infile", "input file with tabulated data: n eff(dist,mean,sigma) bkg(dist,mean,sigma)", false,"","string",cmd);
    ValueArg<int>    fileLines( "l", "infilelines", "max number of lines to be read, read all if < 1", false,0,"int",cmd);
    //
    cmd.parse(argc,argv);
    //
    pole->setVerbose(doVerbose.getValue());
    pole->setPoisson(&PDF::gPoisson);
    pole->setGauss(&PDF::gGauss);
    pole->setGauss2D(&PDF::gGauss2D);
    pole->setLogNormal(&PDF::gLogNormal);

    pole->setInputFile(inputFile.getValue().c_str());
    pole->setInputFileLines(fileLines.getValue());
    pole->setMethod(method.getValue());
    pole->setCL(confLevel.getValue());
    pole->setNObserved(nObs.getValue());
    //
    pole->setEffPdfScale( effScale.getValue() );
    pole->setEffPdf( effMeas.getValue(), effSigma.getValue(), static_cast<PDF::DISTYPE>(effDist.getValue()) );
    pole->setEffObs();

    pole->setBkgPdfScale( bkgScale.getValue() );
    pole->setBkgPdf( bkgMeas.getValue(), bkgSigma.getValue(), static_cast<PDF::DISTYPE>(bkgDist.getValue()) );
    pole->setBkgObs();

    pole->checkEffBkgDists();
    pole->setEffBkgPdfCorr(beCorr.getValue());

    pole->setTrueSignal(sTrue.getValue());
    pole->setUseCoverage(false); //coverage.getValue());

    pole->setBestMuStep(dMus.getValue());
    pole->setBestMuNmax(nMus.getValue());

    pole->setIntEffNSigma(effIntNSigma.getValue());
    pole->setIntBkgNSigma(bkgIntNSigma.getValue());
    pole->setIntGslNCalls(gslIntNCalls.getValue());

    pole->setIntSigRange( tabPoleSMin.getValue(), tabPoleSMax.getValue(), tabPoleSNStep.getValue());
    pole->setIntNobsRange(tabPoleNMin.getValue(), tabPoleNMax.getValue());
    pole->setTabulateIntegral(tabPole.getValue());

    //
    pole->setBSThreshold(threshBS.getValue());
    pole->setPrecThreshold(threshPrec.getValue());
    //
    pole->setHypTestRange(hypTestMin.getValue(), hypTestMax.getValue(), hypTestStep.getValue());

    pole->setMinMuProb(minProb.getValue());
    //
    if (tabPois.getValue()) {
      PDF::gPrintStat = false;
      PDF::gPoisson.initTabulator();
      PDF::gPoisson.setTabMean( tabPoisMuN.getValue(), tabPoisMuMin.getValue(), tabPoisMuMax.getValue() );
      PDF::gPoisson.setTabN(tabPoisNmin.getValue(),tabPoisNmax.getValue());
      PDF::gPoisson.tabulate();
      PDF::gPoisson.clrStat();
    }
  }
  catch (ArgException e) {
    std::cerr << "ERROR: " << e.error() << " " << e.argId() << std::endl;
  }
}
