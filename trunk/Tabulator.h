#ifndef TABULATOR_H
#define TABULATOR_H

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <string>
#include "Combination.h"
#include "ITabulator.h"

/*! @class Tabulator

  @brief A class to tabulate a given function

  Assume you have awrapper class around a function:

  class MyFun {
  public:
    MyFun() {}
    ~MyFun() {}
    double getVal( const double *x ) {
      return x[0]+3.0*x[1]+sin(x[2]);
    }
  }

  To tabulate this function you need to do:

  template<> double Tabulator<MyFun>::calcValue() const {
    // set somehow m_x <= m_parameters
    // and then call getVal() of the class
    return m_function->getVal(&m_parameters);
  }

  The main code will be:

  int main(int argc, char *argv[]) {
    MyFun myfun;              // instance of your function to be tabulated
    Tabulator<MyFun> myTab(); // the tabulator for this class
    myTab.addTabParStep(1.0,200.0,0.1); // first parameter
    myTab.addTabParStep(1.4,20.0,0.01); // second parameter
      ...
    myTab.setFunction(&myfun);          // pointer to actual instance of class
    myTab.tabulate();                   // now tabulate!
  }
 */
template<class T>
class Tabulator : public ITabulator {
public:
   //! main constructor
   inline Tabulator(const char *name, const char *desc=0);
   //! empty constructor
   inline Tabulator();
   //! destructor
   inline virtual ~Tabulator();
   //! set table name
   inline void setName( const char *name );
   //! set table description
   inline void setDescription( const char *desc );
   //! set ptr to function of class T
   inline void setFunction(const T *fun);
   //@}
   /*! @name Tabulating */
   //@{
   //! set number of tabulated parameters
   inline void setNTabPar( size_t npars );
   //! add tabulated parameter def, using step size
   inline void addTabParStep( const char *name, int index, double xmin, double xmax, double step, int parInd=-1 );
   //! add tabulated parameter def, using N(steps)
   inline void addTabParNsteps( const char *name, int index, double xmin, double xmax, size_t nsteps, int parInd=-1 );
   //! clear table
   inline void clrTable();
   //! print table
   inline void printTable() const;
   //! do the tabulation
   inline void tabulate();
   //! get the tabulated value with the given parameter vector
   inline double getValue( const std::vector<double> & valvec );

protected:
   inline void setTabPar( const char *name, int index, double min, double max, double step, size_t nsteps, int parInd=-1 );
   //! initialize table
   inline void initTable();
   //! sets the parameter values given
   inline void setParameters( const std::vector<double> & valvec );
   //! sets the parameter values given (index)
   inline bool setParameters( const std::vector<size_t> & indvec );
   //! sets the parameter values given, using cache
   inline void setParameters( const std::vector<size_t> & indvec, const std::vector<size_t> & indvecLast );
   //! calculate the index in the table for a given vector of values
   inline int calcTabIndex( const std::vector<double> & valvec );
   //! to be called by tabulate() - to be implemented for each specific class
   inline double calcValue() const;
   //! to be called by calcValue() - interpolator - may be either a default function or defined per class
   inline double interpolate( size_t ind ) const;
   
   const T  *m_function;    /**< pointer to function class */
};


#include "Tabulator.icc"

#endif
