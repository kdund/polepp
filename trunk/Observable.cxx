#include <iostream>
#include <vector>
#include <cmath>
#include "Random.h"
#include "Observable.h"


int main(int argc, char *argv[]) {
  Random rndGen;
  GaussPdf gpdf;
  PoissonPdf ppdf;
  ObservableGauss myGaus(&gpdf,&rndGen,"eff");
  ObservablePoisson myPois(&ppdf,&rndGen,"N");
  double a;
  int b;
  //
  ppdf.setLambda(4.5);
  gpdf.setMean(0.0);
  gpdf.setSigma(1.0);
  std::cout << gpdf(0.0)  << std::endl;
  std::cout << gpdf(1.0)  << std::endl;
  std::cout << gpdf(-1.0) << std::endl;
  std::cout << myGaus(-1.0) << std::endl;
  std::cout << myPois(1) << std::endl;
  for (int i=0; i<50; i++) {
    a = myGaus();
    b = myPois();
    std:: cout << "G: " << a << " - P: " << b << std::endl;
  }
  return 0;
}
