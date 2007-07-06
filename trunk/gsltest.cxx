#include <iostream>
#include <cmath>
#include "Integrator.h"
#include "Tools.h"

/* Computation of the integral,
     
I = int (dx dy dz)/(2pi)^3  1/(1-cos(x)cos(y)cos(z))
     
over (-pi,-pi,-pi) to (+pi, +pi, +pi).  The exact answer
is Gamma(1/4)^4/(4 pi^3).  This example is taken from
C.Itzykson, J.M.Drouffe, "Statistical Field Theory -
Volume 1", Section 1.1, p21, which cites the original
paper M.L.Glasser, I.J.Zucker, Proc.Natl.Acad.Sci.USA 74
1800 (1977) */
     
/* For simplicity we compute the integral over the region 
   (0,0,0) -> (pi,pi,pi) and multiply by 8 */
     

inline double gauss( double x, double mu, double sigma ) {
   double s = (x-mu)/sigma;
   return (1.0L/sqrt(2.0*M_PIl)*sigma)*exp(-0.5L*s*s);
}

inline double poisson(int n, double s) {
   double prob;
   double nlnl,lnn,lnf;
   prob = 0.0;
   nlnl = double(n)*log(s);  // n*ln(s)
   lnn  = lgamma(n+1);       // ln(fac(n))
   lnf  = nlnl - lnn - s;
   if (isinf(lnf) || isnan(lnf)) {
      prob=(n==0 ? 1.0:0.0);
   } else {
      prob=exp(lnf);
   }
   if (isnan(prob)) {
      std::cout << "NaN in rawPoisson: " << n << ", " << s << ", " << prob << std::endl;
   }
   return prob;
}

double g(double *k, size_t dim, void *params)
{
   double A = 1.0 / (M_PI * M_PI * M_PI);
   return A / (1.0 - std::cos(k[0]) * std::cos(k[1]) * std::cos(k[2]));
}

double gSignal;
int    gNobs;
double gEffObs;
double gEffErrObs;
double gBkgObs;
double gBkgErrObs;

double gParams[6];
void  *gParamsPtr=&gParams[0];

double poleFun(double *k, size_t dim, void *params)
{
   // k[0] = eff
   // k[1] = bkg
   // params[0] = nobs
   // params[1] = signal truth
   // params[2] = eff obs
   // params[3] = eff obs sigma
   // params[4] = bkg obs
   // params[5] = bkg obs sigma

   //
   double *parptr = static_cast<double *>(params);
   double fe = gauss(k[0], parptr[2] ,parptr[3]);
   double fb = gauss(k[1], parptr[4] ,parptr[5]);
   double lambda = k[0]*parptr[1]  + k[1];
   double fn = poisson(static_cast<int>(parptr[0]+0.5),lambda);
   return fn*fe*fb;
}


double poleFunFast(double *k, size_t dim, void *params)
{
   // k[0] = eff
   // k[1] = bkg
   double fe = 1.0;//gauss(k[0],gEffObs, gEffErrObs);
   double fb = 1.0;//gauss(k[1],gBkgObs, gBkgErrObs);
   double lambda = k[0]*gSignal + k[1];
   double fn = 1.0+lambda;
   return fn*fe*fb;
}

void display_results (char *title, double result, double error, double chisq)
{
   double exact = 1.3932039296856768591842462603255;
   std::cout << "===" << title << "===" << std::endl;
   std::cout << "result  = " << result << std::endl;
   std::cout << "sigma   = " << error << std::endl;
   std::cout << "chisq/N = " << chisq << std::endl;
   std::cout << "diff    = " << std::fabs(result-exact)/error << std::endl;
}

void simpleTst()
{
   std::vector<double> xl(3,0);
   std::vector<double> xu(3,M_PI);
   TOOLS::Timer timer;
   unsigned int ncalls = 1000000;
//    IntegratorVegas vegasInt;
//    vegasInt.setFunction(&g);
//    vegasInt.setFunctionDim(3);
//    vegasInt.setFunctionParams(0);
//    vegasInt.setNcalls(ncalls);
//    vegasInt.setIntRanges(xl,xu);
//    vegasInt.initialize();
//    timer.start();
//    vegasInt.go();
//    timer.stop();
//    display_results ("vegas integration", vegasInt.result(), vegasInt.error(), vegasInt.chisq());
//    timer.printUsedClock(0);

   IntegratorPlain plainInt;
   plainInt.setFunction(&g);
   plainInt.setFunctionDim(3);
   plainInt.setFunctionParams(0);
   plainInt.setNcalls(ncalls);
   plainInt.setIntRanges(xl,xu);
   plainInt.initialize();
   plainInt.tabulate();
   timer.start();
   plainInt.go();
   timer.stop();
   display_results ("plain integration", plainInt.result(), plainInt.error(), plainInt.chisq());
   timer.printUsedClock(0);

   IntegratorMiser miserInt;
   miserInt.setFunction(&g);
   miserInt.setFunctionDim(3);
   miserInt.setFunctionParams(0);
   miserInt.setNcalls(ncalls);
   miserInt.setIntRanges(xl,xu);
   miserInt.initialize();
   timer.start();
   miserInt.go();
   timer.stop();
   display_results ("miser integration", miserInt.result(), miserInt.error(), miserInt.chisq());
   timer.printUsedClock(0);
}

void poleTst()
{
   gSignal = 1.0;
   gNobs   = 1;
   gEffObs = 1.0;
   gEffErrObs = 0.5;
   gBkgObs = 2.0;
   gBkgErrObs = 0.1;
   //
   gParams[0] = 1.0; // nobs
   gParams[1] = 1.0; // signal truth
   gParams[2] = 1.0; // eff
   gParams[3] = 0.5; // eff err
   gParams[4] = 2.0; // bkg
   gParams[5] = 0.1; // bkg err
   //
   const size_t dim=2;
   std::vector<double> xl(dim);
   std::vector<double> xu(dim);
   TOOLS::Timer timer;
   unsigned int ncalls = 1000000;

   xl[0] = gParams[2]-5.0*gParams[3];
   if (xl[0]<0) xl[0]=0;
   xl[1] = gParams[4]-5.0*gParams[5];
   if (xl[1]<0) xl[1]=0;
   xu[0] = gParams[2]+5.0*gParams[3];
   xu[1] = gParams[4]+5.0*gParams[5];

//    IntegratorVegas vegasInt;
//    vegasInt.setFunction(&poleFun);
//    vegasInt.setFunctionDim(dim);
//    vegasInt.setFunctionParams(&gParams[0]);
//    vegasInt.setNcalls(ncalls);
//    vegasInt.setIntRanges(xl,xu);
//    vegasInt.initialize();
//    timer.start();
//    vegasInt.go();
//    timer.stop();
//    display_results ("vegas integration", vegasInt.result(), vegasInt.error(), vegasInt.chisq());
//    timer.printUsedClock(0);

   IntegratorPlain plainInt;
   plainInt.setFunction(&poleFun);
   plainInt.setFunctionDim(dim);
   plainInt.setFunctionParams(&gParams[0]);
   plainInt.setNcalls(ncalls);
   plainInt.setIntRanges(xl,xu);
   plainInt.initialize();
   plainInt.addTabPar(0,0.0,10.0,2.0);
   plainInt.addTabPar(1,0.0,1.0,0.2);
   //   plainInt.addTabPar(2,0.0,5.0,0.2);
   timer.start();
   plainInt.tabulate();
   timer.stop();
   timer.printUsedClock(0);
   timer.start();
   plainInt.go();
   timer.stop();
   display_results ("plain integration", plainInt.result(), plainInt.error(), plainInt.chisq());
   timer.printUsedClock(0);

   IntegratorMiser miserInt;
   miserInt.setFunction(&poleFun);
   miserInt.setFunctionDim(dim);
   miserInt.setFunctionParams(&gParams[0]);
   miserInt.setNcalls(ncalls);
   miserInt.setIntRanges(xl,xu);
   miserInt.initialize();
   timer.start();
   miserInt.go();
   timer.stop();
   display_results ("miser integration", miserInt.result(), miserInt.error(), miserInt.chisq());
   timer.printUsedClock(0);

   IntegratorMiser miserIntFast;
   miserIntFast.setFunction(&poleFunFast);
   miserIntFast.setFunctionDim(dim);
   miserIntFast.setFunctionParams(&gParams[0]);
   miserIntFast.setNcalls(ncalls);
   miserIntFast.setIntRanges(xl,xu);
   miserIntFast.initialize();
   timer.start();
   miserIntFast.go();
   timer.stop();
   display_results ("miser integration, fast", miserIntFast.result(), miserIntFast.error(), miserIntFast.chisq());
   timer.printUsedClock(0);
}     
int main (int argc, char *argv[])
{
   poleTst();
   return 0;
}
