#include <iostream>
#include <fstream>
#include <iomanip>
#include <signal.h>
#include <tclap/CmdLine.h> // Command line parser
#include "Coverage.h"
#include "Tools.h"

Coverage gCoverage;
Pole     gPole;


void my_sighandler(int a) {
  std::string timestamp;
  TOOLS::makeTimeStamp(timestamp);
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
    gCoverage.outputCoverageResult(0); // always output
    gCoverage.calcStatistics();        // only if collecting stats
    gCoverage.printStatistics();       // ditto
    gCoverage.dumpExperiments();       // only if dump filename is defined
    exit(-1);
  }
}

extern void argsCoverage(Coverage *coverage, Pole *pole, int argc, char *argv[]);

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
  argsCoverage(&gCoverage,&gPole,argc,argv);
  //
  //  if (gPole.checkParams()) {
  gCoverage.printSetup();
  gCoverage.doLoop();
    //  }
}
