#include <iostream>
#include <fstream>
#include <iomanip>
#include "Pole.h"

extern void argsPole(Pole *pole, int argc, char *argv[]);

int main(int argc, char *argv[]) {

  Pole pole;
  //
  argsPole(&pole, argc, argv);
  //
  if (pole.checkParams()) {
    pole.printSetup();
    //
    pole.initIntArrays();
    pole.initBeltArrays();
    pole.initIntegral();
    if (!pole.usesMBT()) {
      pole.findAllBestMu(); // loops
    }
    pole.findBelt();
    //
    pole.printLimit(true);
  }
}
