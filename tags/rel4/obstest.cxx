#include <iostream>
#include <fstream>
#include <iomanip>
#include "Random.h"
#include "Pdf.h"
#include "Observable.h"

int main(int argc, char *argv[]) {
  int neve=100;
  int what=0; // generate both sig and bkg
  if (argc>1) {
    neve = atoi(argv[1]);
  }
  if (argc>2) {
    what = atoi(argv[2]);
  }
  if (what>1) what=1;
  if (what<-1) what=-1;
  //
  OBS::ObservableGauss *obsAsig = dynamic_cast<OBS::ObservableGauss *>(OBS::makeObservable(PDF::DIST_GAUS));
  OBS::ObservableFlat  *obsBsig = dynamic_cast<OBS::ObservableFlat  *>(OBS::makeObservable(PDF::DIST_FLAT));
  OBS::ObservableFlat  *obsAbkg = dynamic_cast<OBS::ObservableFlat  *>(OBS::makeObservable(PDF::DIST_FLAT));
  OBS::ObservableGauss *obsBbkg = dynamic_cast<OBS::ObservableGauss *>(OBS::makeObservable(PDF::DIST_GAUS));
  //
  // Signal:
  // A - gauss
  // B - flat
  //
  obsAsig->setPdfMean(0);
  obsAsig->setPdfSigma(1.0);
  obsAsig->setName("A(sig)");

  obsBsig->setPdfMean(0);
  obsBsig->setPdfSigma(5.0);
  obsBsig->setName("B(sig)");
  //
  // Background:
  // A - flat
  // B - gauss
  //
  obsAbkg->setPdfMean(0);
  obsAbkg->setPdfSigma(5.0);
  obsAbkg->setName("A(bkg)");

  obsBbkg->setPdfMean(0);
  obsBbkg->setPdfSigma(1.0);
  obsBbkg->setName("B(bkg)");

  std::cout << "event/I:A/F:B/F:class/I" << std::endl;
  int i=0;
  if (what!=-1) {
    for (i=0; i<neve; i++) {
      std::cout << i << "\t"
                << obsAsig->rnd() << "\t"
                << obsBsig->rnd() << "\t"
                << 1
                << std::endl;
    }
  }

  if (what!=1) {
    int j=0;
    for ( j=0; j<neve; j++) {
      std::cout << i << "\t"
                << obsAbkg->rnd() << "\t"
                << obsBbkg->rnd() << "\t"
                << -1
                << std::endl;
      i++;
    }
  }
  return 0;
}
