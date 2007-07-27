#include <iostream>
#include <fstream>
#include <iomanip>
#include <sys/times.h>

#include "Tools.h"
#include "Pdf.h"
#include "Random.h"

RND::Random gRnd;

using namespace std;

void print_time(struct tms *start, struct tms *stop) {
  cout << "User time: " << stop->tms_utime-start->tms_utime << endl;
}

// double checkPoisTime(bool usetab) {
//   static struct tms start,stop;
//   //
//   double rval=0;
//   double lmb;
//   const double lmbmax = 30.0;
//   const int nlmb=200000;
//   double dlmb = lmbmax/double(nlmb);
//   double p;
//   int n,i;
//   //
//   PDF::Poisson *pdf;
//   if (usetab) {
//     pdf = &PDF::gPoisTab;
//     cout << "Timing using tabulated poisson." << endl;
//   } else {
//     pdf = &PDF::gPoisson;
//     cout << "Timing NOT using tabulated poisson." << endl;
//   }
//   times(&start);
//   for (n=0; n<40; n++) {
//     for (i=0; i<nlmb; i++) {
//       lmb = dlmb*double(i);
//       p = pdf->getVal(n,lmb);
//       rval +=p;
//     }
//   }
//   times(&stop);
//   print_time(&start,&stop);
//   return rval;
// }

void checkPois(PDF::Base * poisTab, PDF::Base * poisNoTab ) {

  //
  double p,ptab,s=0;
  double dp;
  double dpsum;
  double dpsum2;
  double dpmax,dpmaxn;
  double ns;
  double smax,smaxn;
  double pmax,pmaxn;
  int nmax=0;
  int nr;
  double ns0  = 726.0;
  double sns0 = 10.0;
  //
  dpmax=0;
  pmax=0;
  smax=0;
  nr=0;
  int nsig = 1000;
  for (int n=0; n<100; n++) {
    dpsum=dpsum2=0;
    smaxn=0;
    dpmaxn=0;
    pmaxn = 0;
    ns = RND::gRandom.gauss(ns0,sns0);
    for (int i=0; i<nsig; i++) {
      nr = RND::gRandom.poisson(ns);
      p = poisNoTab->getVal(nr,ns);
      ptab = poisTab->getVal(nr,ns);
      dp = fabs(2.0*100.0*(p-ptab)/(p+ptab));
      //
      //      if (p>0.0) std::cout << "DIFF = " << p << std::endl;
        dpsum  +=dp;
        dpsum2 +=dp*dp;
	if (dp>dpmax) {
	  dpmax=dp;
	  smax = s;
	  nmax = nr;
	}
	if (dp>dpmaxn) {
	  dpmaxn = dp;
	  smaxn = s;
	  pmaxn = p;
	}
        //      }
//       cout << "TAB: "; // << i << "   ";
//       TOOLS::coutFixed(4,nr); cout << "   ";
//       TOOLS::coutFixed(12,ns); cout << "   ";
//       TOOLS::coutFixed(12,p); cout << "   ";
//       TOOLS::coutFixed(12,ptab); cout << "   ";
//       TOOLS::coutFixed(12,dp); cout << endl;
    }
//     cout << "<"; TOOLS::coutFixed(3,nr); cout << "> : ";
//     TOOLS::coutFixed(6,dpmaxn);
//     cout << "%    s = ";
//     TOOLS::coutFixed(6,smaxn);
//     cout << "     p = ";
//     TOOLS::coutFixed(6,pmaxn);
//     cout << "     v = ";
//     TOOLS::coutFixed(6,(dpsum2-(dpsum*dpsum/double(nsig)))/(double(nsig)-1.0));
//     cout << endl;
    //    
  }
  cout << "Overall max : " << dpmax << "% " << nmax << " " << smax << std::endl;
}

void timePois(PDF::Base * pois ) {
  //
  static struct tms start,stop;
  double ns0  = 650.0;
  double sns0 = 10.0;
  //
  int nsig = 10000;
  times(&start);
  for (int n=0; n<100; n++) {
    double ns = RND::gRandom.gauss(ns0,sns0);
    for (int i=0; i<nsig; i++) {
      int nr = RND::gRandom.poisson(ns);
      double p = pois->getVal(nr,ns);
      if (p<0) std::cout << "*" << std::endl;
    }
  }
  times(&stop);
  print_time(&start,&stop);
}

// void cmpPois(PDF::Poisson & poisA, PDF::Poisson & poisB) {
//   const double *dataA = poisA.getData();
//   const double *dataB = poisB.getData();
//   int nA = poisA.getNdata();
//   int nB = poisB.getNdata();
//   if (nA!=nB) {
//     cout << "Poisson tables do not have the same amount of data points" << endl;
//     return;
//   }
//   //
//   double pA,pB;
//   double dp;
//   double dpmax=-1.0;
//   int imax=0;
//   for (int i=0; i<nA; i++) {
//     pA = dataA[i];
//     pB = dataB[i];
//     dp = fabs(2.0*100.0*(pA-pB)/(pA+pB));
//     //    cout << i << "   ";
//     //    TOOLS::coutFixed(4,dp); cout << endl;
//     if (dp>dpmax) {
//       dpmax = dp;
//       imax = i;
//     }
//   }
//   cout << "Max diff (%) = " << dpmax << endl;
// }

void timePoisNoTab(PDF::Poisson & poisson) {
  cout << "Timing of poisson with no tabulation, getValue" << endl;
  timePois(&poisson);
}

void timePoisNewTab( PDF::Poisson & poisTab) {
  static struct tms start,stop;

  cout << "Timing of new tabulated poisson, init" << endl;
  times(&start);
  poisTab.initTabulator();
  poisTab.setTabN( 400,1000 );
  poisTab.setTabMean(1202, 400.0,1000.0 );
  poisTab.tabulate();
  times(&stop);
  print_time(&start,&stop);
  cout << "Timing of new tabulated poisson, getValue" << endl;
  timePois(&poisTab);
}

void timePoisOldTab(  PDF::PoisTab & poisTabOld ) {

  static struct tms start,stop;
  //
  cout << "Timing of old tabulated Poisson, init" << endl;
  times(&start);
  poisTabOld.setRangeMean( 1202,400.0,1000.0);
  poisTabOld.setRangeX(601,400);
  poisTabOld.tabulateOld();
  times(&stop);
  print_time(&start,&stop);
  cout << "Timing of old tabulated Poisson, getValue" << endl;
  timePois(&poisTabOld);
}

int main(int argc, char *argv[]) {
  PDF::Poisson poisTab;
  PDF::Poisson poisNoTab;
  PDF::Poisson poisDum;
  PDF::PoisTab poisTabOld(&poisDum);
  PDF::gPrintStat = true;

  timePoisNoTab(poisNoTab);

  timePoisNewTab(poisTab);

  timePoisOldTab(poisTabOld);

  poisTab.clrStat();
  poisNoTab.clrStat();
  checkPois( &poisTab, &poisNoTab );
  std::cout << "\nTABULATED:" << std::endl;
  poisTab.printStat();
  std::cout << "\nDIRECT:" << std::endl;
  poisNoTab.printStat();
  PDF::gPrintStat = false;
}
