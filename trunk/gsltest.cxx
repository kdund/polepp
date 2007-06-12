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
     

double g(double *k, size_t dim, void *params)
{
   double A = 1.0 / (M_PI * M_PI * M_PI);
   return A / (1.0 - std::cos(k[0]) * std::cos(k[1]) * std::cos(k[2]));
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
     
int main (int argc, char *argv[])
{
   std::vector<double> xl(3,0);
   std::vector<double> xu(3,M_PI);
   TOOLS::Timer timer;
   unsigned int ncalls = 1000000;
   IntegratorVegas vegasInt;
   vegasInt.setFunction(&g);
   vegasInt.setFunctionDim(3);
   vegasInt.setFunctionParams(0);
   vegasInt.setNcalls(ncalls);
   vegasInt.setIntRanges(xl,xu);
   vegasInt.initialize();
   timer.start();
   vegasInt.go();
   timer.stop();
   display_results ("vegas integration", vegasInt.result(), vegasInt.error(), vegasInt.chisq());
   timer.printUsedClock(0);

   IntegratorPlain plainInt;
   plainInt.setFunction(&g);
   plainInt.setFunctionDim(3);
   plainInt.setFunctionParams(0);
   plainInt.setNcalls(ncalls);
   plainInt.setIntRanges(xl,xu);
   plainInt.initialize();
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

   return 0;
}
