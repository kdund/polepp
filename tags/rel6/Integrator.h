#ifndef INTEGRATOR_H
#define INTEGRATOR_H
#include <iostream>
#include <vector>
#include <cmath>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_interp.h>


class Integrator {
 public:
   inline Integrator();
   inline virtual ~Integrator();

   inline void setFunction( double (* f)(double * x, size_t dim, void * params) );
   inline void setFunctionDim( size_t dim );
   inline void setFunctionParams( void * params );
   inline void setIntRanges( std::vector<double> & xl, std::vector<double> & xu );
   inline void setNcalls( unsigned int nc );
   inline void addTabPar( size_t i, double xmin, double xmax, double step );
   inline void setTabParamsVal( std::vector<int> & tabind );
   //
   inline void clrTable();

   inline double getValue( std::vector<double> & tabParams );
   //
   inline virtual void initialize();
   inline virtual void initTable();
   inline virtual void tabulate();
   inline virtual void go( void *params );
   inline virtual double chisq()=0;
   inline double result();
   inline double error();

 protected:
   gsl_rng            *m_gslRange;    // GSL range
   unsigned int        m_ncalls;      // number of iterations
   gsl_monte_function  m_gslMonteFun; // structure for GSL MC integration
   std::vector<double> m_intXL;       // lower integration range for each integrand
   std::vector<double> m_intXU;       // idem, upper

   std::vector<size_t> m_tabIndex;    // index of tabulated var
   std::vector<double> m_tabMin;      // minimum of tabulated value
   std::vector<double> m_tabMax;      // maximum
   std::vector<double> m_tabStep;     // step size
   std::vector<int>    m_tabNsteps;   // number of steps
   std::vector<double> m_tabXvalues;  // tabulated X for each parameter
   std::vector<int>    m_tabXvalInd;  // first index in above vector for param [i]
   std::vector<double> m_tabValues;   // the actual table
   bool                m_tabulated;   // true if tabulate() is called successfully


   double m_result;                   // result of integration
   double m_error;                    // error of idem
};

class IntegratorVegas : public Integrator {
public:
   inline IntegratorVegas();
   inline virtual ~IntegratorVegas();
   inline virtual void go( void *params=0 );
   inline virtual void initialize();
   inline virtual double chisq();
private:
   gsl_monte_vegas_state *m_gslVegasState;
};

class IntegratorPlain : public Integrator {
public:
   inline IntegratorPlain();
   inline virtual ~IntegratorPlain();
   inline virtual void go( void *params=0 );
   inline virtual void initialize();
   inline virtual double chisq();
private:
   gsl_monte_plain_state *m_gslPlainState;
};

class IntegratorMiser : public Integrator {
public:
   inline IntegratorMiser();
   inline virtual ~IntegratorMiser();
   inline virtual void go( void *params=0 );
   inline virtual void initialize();
   inline virtual double chisq();
private:
   gsl_monte_miser_state *m_gslMiserState;
};

#include "Integrator.icc"

#endif
