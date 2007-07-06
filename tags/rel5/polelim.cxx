#include <iostream>
#include <fstream>
#include <iomanip>
#include "Pole.h"

extern void argsPole(Pole *pole, int argc, char *argv[]);

int main(int argc, char *argv[]) {

  Pole pole;
  //
  //  processArgs(&pole, argc, argv);
  argsPole(&pole, argc, argv);

  pole.execute();

//   pole.initAnalysis();
//   //
//   //  if (pole.checkParams()) {
//   pole.printSetup();
//   if (pole.analyseExperiment()) {
//     pole.printLimit(true);
//   } else {
//     pole.printFailureMsg();
//   }
//   //  }
}
