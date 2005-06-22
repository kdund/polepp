#include <iostream>
#include <vector>
#include <cmath>
#include "Random.h"
#include "Observable.h"


int main(int argc, char *argv[]) {
  Random rndGen;
  PdfGauss gpdf;
  PdfPoisson ppdf;
  ObservableGauss myGaus(&gpdf,&rndGen,"eff");
  ObservablePoisson myPois(&ppdf,&rndGen,"N");
  double a,fa,fb;
  int b;
  //
  ppdf.setLambda(4.5);
  std::cout << "Tabulating Poisson..." << std::flush;
  ppdf.setTableParams(10000,10.0,50);
  ppdf.tabulate();
  std::cout << "done" << std::endl;
  gpdf.setMean(0.0);
  gpdf.setSigma(1.0);
  std::cout << "Tabulating Gauss..." <<std::flush;
  gpdf.setTableParams(10000,5.0);
  gpdf.tabulate();
  std::cout << "done" << std::endl;
  std::cout << gpdf(0.0)  << std::endl;
  std::cout << gpdf(1.0)  << std::endl;
  std::cout << gpdf(-1.0) << std::endl;
  std::cout << myGaus(-1.0) << std::endl;
  std::cout << myPois(1) << std::endl;
  for (int i=0; i<50; i++) {
    a = myGaus();
    b = myPois();
    fa = myGaus(a);
    fb = myPois(b);
    std:: cout << "G: " << a << " " << fa << " - P: " << b << " " << fb << std::endl;
  }
  return 0;
}
