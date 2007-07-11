#include <iostream>
#include <fstream>
#include <iomanip>
#include "Random.h"
#include "Pdf.h"
#include "Observable.h"

const char gDelim=' ';

double getPositive( OBS::ObservableGauss *obs ) {
   double r=-1.0;
   while (r<0) {
      r=obs->rnd();
   }
   return r;
}

void pixelOccupancy( int neve, double prob ) {
   const int npixTot=1000;
   OBS::ObservablePois *npix    = dynamic_cast<OBS::ObservablePois *>(OBS::makeObservable(PDF::DIST_POIS));
   OBS::ObservableFlat *whatpix = dynamic_cast<OBS::ObservableFlat *>(OBS::makeObservable(PDF::DIST_FLAT));
   double mean = npixTot*prob;
   std::cout << "N(events)       = " << neve << std::endl;
   std::cout << "N(pixels)       = " << npixTot << std::endl;
   std::cout << "N(bad pix) mean = " << mean << std::endl;
   npix->setPdfMean(mean);
   whatpix->setPdfRange(0,npixTot-1);
   int npixObs;
   int pixInd;
   int npsum=0;
   double p0sum=0;
   double p1sum=0;
   double p1psum=0;
   std::vector<int> occupancy(npixTot,0);
   for (int i=0; i<neve; i++) { // loop over events
      npixObs = npix->rnd();    // get N bad pixels
      if (npixObs>npixTot) npixObs=npixTot;
      npsum += npixObs;
      for (int p=0; p<npixObs; p++) { // fill bad pixels
         pixInd = static_cast<int>(whatpix->rnd()+0.5);
         occupancy[pixInd]++;
      }      
      int n0=0;
      int n1=0;
      int n1p=0;
      for (int i=0; i<npixTot; i++) {
         if (occupancy[i]==0) {
            n0++;
         } else {
            n1p++;
            if (occupancy[i]==1) n1++;
         }
      }
      p0sum  += double(n0)/double(npixTot);
      p1sum  += double(n1)/double(npixTot);
      p1psum += double(n1p)/double(npixTot);
      occupancy.clear();
      occupancy.resize(npixTot,0);
   }
   
   p0sum /= double(neve);
   p1sum /= double(neve);
   p1psum /= double(neve);

   double praw = exp(log(p0sum)/double(mean-1));
   std::cout << "average number of bad pixels = " << double(npsum)/double(neve) << std::endl;
   std::cout << "p(0)  = " << p0sum << std::endl;
   std::cout << "p(1)  = " << p1sum << std::endl;
   std::cout << "p(>1) = " << p1psum << std::endl;
   std::cout << "praw  = " << 1.0-praw << std::endl;
}

void makeToyMC( int neve, int what ) {
  OBS::ObservableGauss *mass     = dynamic_cast<OBS::ObservableGauss *>(OBS::makeObservable(PDF::DIST_GAUS));
  OBS::ObservableLogN  *massLN   = dynamic_cast<OBS::ObservableLogN  *>(OBS::makeObservable(PDF::DIST_LOGN));
  OBS::ObservableFlat  *theta    = dynamic_cast<OBS::ObservableFlat  *>(OBS::makeObservable(PDF::DIST_FLAT));
  OBS::ObservableGauss *bkgY1    = dynamic_cast<OBS::ObservableGauss *>(OBS::makeObservable(PDF::DIST_GAUS));
  OBS::ObservableGauss *bkgY2    = dynamic_cast<OBS::ObservableGauss *>(OBS::makeObservable(PDF::DIST_GAUS));
  OBS::ObservableGauss *detector = dynamic_cast<OBS::ObservableGauss *>(OBS::makeObservable(PDF::DIST_GAUS));

  
  const double realBkgProb=0.6;
  const double realSigProb=1.0-realBkgProb;
  const double realMissClassProb=0.01;

  double bkgProb;
  double sigProb;
  double missClassProb;

  switch(what) {
  case -1: // only bkg
     bkgProb = 1.0;
     sigProb = 0.0;
     missClassProb = 0.0;
     break;
  case 0:  // 'real' data
     bkgProb = realBkgProb;
     sigProb = realSigProb;
     missClassProb = realMissClassProb;
     break;
  case 1: // only sig
     bkgProb = 0.0;
     sigProb = 1.0;
     missClassProb = 0.0;
     break;
  default:
     bkgProb = realBkgProb;
     sigProb = realSigProb;
     missClassProb = realMissClassProb;
     break;
  }
        
  //
  // true mass - very sharp spike
  //
  mass->setPdfMean(1.5);
  mass->setPdfSigma(1.0);
  mass->setName("M");
  //
  // signal (y1,y2) : m^2 = y1^2+y2^2
  //
  // (y1,y2)=m*(cos(t),sin(t))
  //  t := Flat(0,pi/2), y1,y2>0
  // => sigma = pi/2/sqrt(12)
  const double mean  = M_PI/4;
  const double sigma = mean*2.0/sqrt(12.0);
  theta->setPdfMean(mean);
  theta->setPdfSigma(sigma);
  theta->setName("theta");
  //
  // bkg of (y1,y2) - independent gauss
  //
  bkgY1->setPdfMean(0);
  bkgY1->setPdfSigma(0.5);
  bkgY1->setName("YB1");

  bkgY2->setPdfMean(0);
  bkgY2->setPdfSigma(0.5);
  bkgY2->setName("YB2");

  //
  // detector smearing: gaussian with zero mean
  //
  detector->setPdfMean(0.0);
  detector->setPdfSigma(0.2);
  detector->setName("detector");
  //
  //
  //
  bool isMissClassified;
  bool isSignal;
  double massT, thetaT;
  double y1T, y2T;
  double x1, x2;
  int  classid, missid;

  std::cout << "event/I:mass/F:y1/F:y2/F:x1/F:x2/F:missid/I:classid/I" << std::endl;
  for (int i=0; i<neve; i++) {
     isSignal = ( RND::gRandom.rndm() < sigProb );
     missid = 0;
     massT = 0;
     thetaT = 0;
     if (isSignal) { // signal event
        classid = 1;
        isMissClassified = ( RND::gRandom.rndm() < missClassProb );
        massT  = mass->rnd();
        thetaT = theta->rnd();
        if (isMissClassified) {
           missid = 1;
           if (RND::gRandom.rndm()>0.5) { // one of (y1,y2) is bkg; which one???
              y1T = getPositive(bkgY1);
              y2T = massT*sin(thetaT);
           } else {
              y1T = massT*cos(thetaT);
              y2T = getPositive(bkgY2);
           }
        } else { // both correctly identified
           missid = 0;
           y1T = massT*cos(thetaT);
           y2T = massT*sin(thetaT);
        }
     } else { // bkg event
        classid = -1;
        y1T = getPositive(bkgY1);
        y2T = getPositive(bkgY2);
     }
     //
     // apply detector smearing
     //
     x1 = -1.0;
     while (x1<0) {
        x1 = y1T + detector->rnd();
     }
     x2 = -1.0;
     while (x2<0) {
        x2 = y2T + detector->rnd();
     }

     std::cout << i << gDelim
               << massT << gDelim
               << y1T << gDelim
               << y2T << gDelim
               << x1 << gDelim
               << x2 << gDelim
               << missid << gDelim
               << classid << gDelim
               << std::endl;
  }
}

void makeTOT(int neve) {
  OBS::ObservableGauss *charge = dynamic_cast<OBS::ObservableGauss *>(OBS::makeObservable(PDF::DIST_GAUS));
  charge->setPdfMean(20000.0);
  charge->setPdfSigma(2000.0);
  charge->setName("Q");
  //
  std::cout << "event/I:charge/F:tot/F" << std::endl;
  const double A = 70.2;
  const double B = -2.0752e6;
  const double C = 26000.0;
  const double p1 = -0.68;
  const double p2 = 0.17;

  double q,t;
  int tot;
  for (int i=0; i<neve; i++) {
    q = charge->rnd(); 
    t = A + (B/(q+C));
    tot = static_cast<int>(t+0.5);
    std::cout << i << gDelim
	      << q << gDelim
	      << tot << gDelim
	      << std::endl;
  }
}

void makeSampleGG(int neve) {
  OBS::ObservableGauss *obsA = dynamic_cast<OBS::ObservableGauss *>(OBS::makeObservable(PDF::DIST_GAUS));
  OBS::ObservableGauss *obsB = dynamic_cast<OBS::ObservableGauss *>(OBS::makeObservable(PDF::DIST_GAUS));
  OBS::ObservableFlat  *obsC = dynamic_cast<OBS::ObservableFlat  *>(OBS::makeObservable(PDF::DIST_FLAT));
  //
  // A - gauss
  // B - gauss
  // C - flat
  //
  obsA->setPdfMean(0);
  obsA->setPdfSigma(1.0);
  obsA->setName("A");

  obsB->setPdfMean(0);
  obsB->setPdfSigma(5.0);
  obsB->setName("B");

  obsC->setPdfMean(10.0);
  obsC->setPdfSigma(1.0);
  obsC->setName("C");

  std::cout << "event/I:AB/F:AC/F" << std::endl;
  int i=0;
  for (i=0; i<neve; i++) {
     std::cout << i << gDelim
               << obsA->rnd()*obsB->rnd() << gDelim
               << obsA->rnd()*obsC->rnd() << gDelim
               << std::endl;
  }
}

void makeSample1(int neve, int what) {
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
      std::cout << i << gDelim
                << obsAsig->rnd() << gDelim
                << obsBsig->rnd() << gDelim
                << 1
                << std::endl;
    }
  }

  if (what!=1) {
    int j=0;
    for ( j=0; j<neve; j++) {
      std::cout << i << gDelim
                << obsAbkg->rnd() << gDelim
                << obsBbkg->rnd() << gDelim
                << -1
                << std::endl;
      i++;
    }
  }
}

void makeSample2(int neve, int what) {
  OBS::ObservableGauss *obsAsig = dynamic_cast<OBS::ObservableGauss *>(OBS::makeObservable(PDF::DIST_GAUS));
  OBS::ObservableGauss *obsBsig = dynamic_cast<OBS::ObservableGauss *>(OBS::makeObservable(PDF::DIST_GAUS));
  OBS::ObservableGauss *obsAbkg = dynamic_cast<OBS::ObservableGauss *>(OBS::makeObservable(PDF::DIST_GAUS));
  OBS::ObservableGauss *obsBbkg = dynamic_cast<OBS::ObservableGauss *>(OBS::makeObservable(PDF::DIST_GAUS));
  //
  // Signal:
  // A - gauss
  // B - gauss
  //
  obsAsig->setPdfMean(0.0);
  obsAsig->setPdfSigma(1.0);
  obsAsig->setName("A(sig)");

  obsBsig->setPdfMean(0.0);
  obsBsig->setPdfSigma(1.0);
  obsBsig->setName("B(sig)");
  //
  // Background:
  // A - gauss
  // B - gauss
  //
  obsAbkg->setPdfMean(1.0);
  obsAbkg->setPdfSigma(1.0);
  obsAbkg->setName("A(bkg)");

  obsBbkg->setPdfMean(1.0);
  obsBbkg->setPdfSigma(1.0);
  obsBbkg->setName("B(bkg)");

  //  std::cout << "event/I:A/F:B/F:class/I" << std::endl;
  std::cout << "Ab/F:Bb/F:As/F:Bs/F" << std::endl;
  int i=0;
  int j=0;
  for ( j=0; j<neve; j++) {
     std::cout << obsAbkg->rnd() << gDelim
               << obsBbkg->rnd() << gDelim
               << obsAsig->rnd() << gDelim
               << obsBsig->rnd()
               << std::endl;
     i++;
  }
}

int main(int argc, char *argv[]) {
  int neve=100;
  int what=0; // generate both sig and bkg
  double val=0;
  if (argc>1) {
    neve = atoi(argv[1]);
  }
  if (argc>2) {
    what = atoi(argv[2]);
  }
  if (argc>3) {
    val = atof(argv[3]);
  }
  if (what>1) what=1;
  if (what<-1) what=-1;

  pixelOccupancy(neve,val);
  //  makeToyMC(neve,what);
  //  makeSample1(neve,what);
  //  makeSample2(neve,what);
  //  makeSampleGG(neve);
  //  makeTOT(neve);
  //
  return 0;
}
