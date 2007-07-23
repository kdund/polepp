#ifndef ITABULATOR_H
#define ITABULATOR_H

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <string>

/*! @class ITabulator

  @brief Interface class to the Tabulator class

*/
class ITabulator {
public:
   //! main constructor
   inline ITabulator(const char *name, const char *desc=0) {
      if (name) m_name        = name;
      if (desc) m_description = desc;
   }
   //! empty constructor
   inline ITabulator() {}
   //! destructor
   inline virtual ~ITabulator() {}

   //! set table name
   virtual void setName( const char *name ) = 0;
   //! set table description
   virtual void setDescription( const char *desc ) = 0;

   //@}
   /*! @name Tabulating */
   //@{
   virtual void setTabNPar( size_t npars ) = 0;
   //! add tabulated parameter def, using step size
   virtual void addTabParStep( const char *name, int index, double xmin, double xmax, double step, int parInd=-1 ) = 0;
   //! add tabulated parameter def, using N(steps)
   virtual void addTabParNsteps( const char *name, int index, double xmin, double xmax, size_t nsteps, int parInd=-1 ) = 0;
   //! clear table
   virtual void clrTable() = 0;
   //! print table info
   virtual void printTable() const = 0;
   //! set tabulator verbosity
   virtual void setVerbose(bool v) = 0;
   //! do the tabulation
   virtual void tabulate() = 0;
   //! get the tabulated value with the given parameter vector
   virtual double getValue( const std::vector<double> & tabParams ) = 0;

   //! accessors
   virtual const char *getName()        const = 0;
   virtual const char *getDescription() const = 0;
   virtual size_t getTabNPar()    const = 0;
   virtual double getTabMin( size_t pind ) const = 0;
   virtual double getTabMax( size_t pind ) const = 0;
   virtual double getTabStep( size_t pind ) const = 0;
   virtual size_t getTabNsteps( size_t pind ) const = 0;

   //! check if the table is ok
   virtual bool isTabulated() const = 0;

protected:
   //! set tabulated par
   virtual void setTabPar( const char *name, int index, double min, double max, double step, size_t nsteps, int parInd=-1 ) = 0;
   //! initialize table
   virtual void initTable() = 0;
   //! sets the parameter values given
   virtual void setParameters( const std::vector<double> & valvec ) = 0;
   //! sets the parameter values given (index)
   virtual bool setParameters( const std::vector<size_t> & indvec ) = 0;
   //! sets the parameter values given
   virtual void setParameters( const std::vector<size_t> & indvec, const std::vector<size_t> & indvecLast) = 0;
   //! calculate the index in the table for a given vector of values
   virtual int calcTabIndex( const std::vector<double> & valvec ) = 0;
   //! calculate the parameter index for the given table index
   virtual int calcParIndex( const size_t tabind, const size_t parind ) const = 0;
   //! to be called by tabulate() - to be implemented for each specific class
   virtual double calcValue() = 0;
   //! to be called by getValue() - interpolator - may be either a default function or defined per class
   virtual double interpolate( size_t ind ) const = 0;
   //! first derivative
   virtual double deriv( size_t tabind, size_t parind ) const = 0;
   //! second derivative
   virtual double deriv2( size_t tabind, size_t parind ) const = 0;

   //
   bool                m_verbose;     /**< verbose flag */
   std::vector<int>    m_tabIndex;    /**< vector of parameter indecis - not used internally, just for external book keeping */
   std::vector<std::string> m_tabName;/**< vector of parameter names - not used internally, just for external book keeping */
   std::vector<double> m_tabMin;      /**< minimum of tabulated value */
   std::vector<double> m_tabMax;      /**< maximum */
   std::vector<double> m_tabStep;     /**< step size */
   std::vector<size_t> m_tabNsteps;   /**< number of steps */
   std::vector<size_t> m_tabMaxInd;   /**< number of steps - 1 ; not so nice */
   std::vector<size_t> m_tabPeriod;   /**< parameter period */
   std::vector<size_t> m_tabNTabSteps;/**< parameter: number of steps in table until next value tabNStep = tabPeriod[i-1]  */
   std::vector<double> m_tabValues;   /**< the actual table */
   bool                m_tabulated;   /**< true if tabulate() is called successfully */

   std::string         m_name;        /**< name */
   std::string         m_description; /**< brief description */

   std::vector<double> m_parameters;  /**< vector containing the parameters set in tabulate() and used in calcValue() */
   std::vector<bool>   m_parChanged;  /**< flags which parameters were changed since last vector in tabulate()        */
   std::vector<size_t> m_parIndex;    /**< indecis obtained by calcTabIndex() */

};

#endif
