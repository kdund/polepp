#include <iostream>
#include <fstream>
#include <iomanip>
#include <sys/times.h>

#include "Pole.h"
#include "Pdf.h"
#include "Random.h"

Random gRnd;

using namespace std;

void print_time(struct tms *start, struct tms *stop) {
  cout << "User time: " << stop->tms_utime-start->tms_utime << endl;
}

PDF::Poisson gPois;
PDF::Poisson gPoisTab;

double checkPoisTime(bool usetab) {
  static struct tms start,stop;
  //
  double rval=0;
  double lmb;
  const double lmbmax = 30.0;
  const int nlmb=200000;
  double dlmb = lmbmax/double(nlmb);
  double p;
  int n,i;
  //
  PDF::Poisson *pdf;
  if (usetab) {
    pdf = &gPoisTab;
    cout << "Timing using tabulated poisson." << endl;
  } else {
    pdf = &gPois;
    cout << "Timing NOT using tabulated poisson." << endl;
  }
  times(&start);
  for (n=0; n<40; n++) {
    for (i=0; i<nlmb; i++) {
      lmb = dlmb*double(i);
      p = pdf->getVal(n,lmb);
      rval +=p;
    }
  }
  times(&stop);
  print_time(&start,&stop);
  return rval;
}

void checkPois() {

  //
  double p,ptab,s;
  double dp;
  double dpsum;
  double dpsum2;
  double dpmax,dpmaxn;
  double ns;  double smax,smaxn;
  double pmax,pmaxn;
  int nmax;
  //
  dpmax=0;
  pmax=0;
  for (int n=0; n<20; n++) {
    dpsum=dpsum2=0;
    ns=0;
    smaxn=0;
    dpmaxn=0;
    pmaxn = 0;
    for (int i=0; i<20000; i++) {
      s = 0.0012450178945671*static_cast<double>(i);
      p = gPois.getVal(n,s);
      ptab = gPoisTab.getVal(n,s);
      dp = fabs(2.0*100.0*(p-ptab)/(p+ptab));
      dpsum  +=dp;
      dpsum2 +=dp*dp;
      if (p>0.0001) {
	if (dp>dpmax) {
	  dpmax=dp;
	  smax = s;
	  nmax = n;
	}
	if (dp>dpmaxn) {
	  dpmaxn = dp;
	  smaxn = s;
	  pmaxn = p;
	}
      }
//       cout << ">>>" << i << "<<< ";
//       coutFixed(4,n); cout << "   ";
//       coutFixed(4,s); cout << "   ";
//       coutFixed(4,p); cout << "   ";
//       coutFixed(4,ptab); cout << "   ";
//       coutFixed(4,dp); cout << "   " << endl;
    }
    cout << "<"; coutFixed(3,n); cout << "> : ";
    coutFixed(4,dpmaxn);
    cout << "%    s = ";
    coutFixed(4,smaxn);
    cout << "     p = ";
    coutFixed(4,pmaxn);    cout << endl;
    //    
  }
  cout << "Overall max : " << dpmax << "% " << nmax << " " << smax << std::endl;
}

int main(int argc, char *argv[]) {
  //
  //  gPois.init(1000000,60,200.0);
  gPoisTab.init(100000,60,30.0);
  //
  //  checkPois();
  cout << " p = " << checkPoisTime(false) << endl;;
  cout << " p = " << checkPoisTime(true) << endl;;
  //
}
