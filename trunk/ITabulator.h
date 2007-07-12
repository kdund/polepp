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
   //@}
   /*! @name Tabulating */
   //@{
   virtual void setNparams( size_t npars ) {};
   //! add tabulated parameter def, using step size
   virtual void addTabParStep( const char *name, int index, double xmin, double xmax, double step, int parInd=-1 ) {}
   //! add tabulated parameter def, using N(steps)
   virtual void addTabParNsteps( const char *name, int index, double xmin, double xmax, size_t nsteps, int parInd=-1 ) {}
   //! clear table
   virtual void clrTable() {}
   //! do the tabulation
   virtual void tabulate() {}
   //! get the tabulated value with the given parameter vector
   virtual double getValue( const std::vector<double> & tabParams ) { return 0;}

   //! check if the table is ok
   bool isTabulated() { return m_tabulated; }

protected:
   //! set tabulated par
   virtual void setTabPar( const char *name, int index, double min, double max, double step, size_t nsteps, int parInd=-1 ) {}
   //! initialize table
   virtual void initTable() {}
   //! sets the parameter values given
   virtual void setParameters( const std::vector<double> & valvec ) {}
   //! sets the parameter values given (index)
   virtual bool setParameters( const std::vector<size_t> & indvec ) { return 0; }
   //! sets the parameter values given
   virtual void setParameters( const std::vector<size_t> & indvec, const std::vector<size_t> & indvecLast) {}
   //! calculate the index in the table for a given vector of values
   virtual int calcTabIndex( const std::vector<double> & valvec ) const { return 0; }
   //! to be called by tabulate() - to be implemented for each specific class
   virtual double calcValue() const {return 0;}

   //
   std::vector<double> m_parameters;  /**< vector containing the parameters set in tabulate() and used in calcValue() */
   std::vector<double> m_parChanged;  /**< flags which parameters were changed since last vector in tabulate()        */
   std::vector<int>    m_tabIndex;    /**< vector of parameter indecis - not used internally, just for external book keeping */
   std::vector<std::string> m_tabName;/**< vector of parameter names - not used internally, just for external book keeping */
   std::vector<double> m_tabMin;      /**< minimum of tabulated value */
   std::vector<double> m_tabMax;      /**< maximum */
   std::vector<double> m_tabStep;     /**< step size */
   std::vector<size_t> m_tabNsteps;   /**< number of steps */
   std::vector<double> m_tabValues;   /**< the actual table */
   bool                m_tabulated;   /**< true if tabulate() is called successfully */

   std::string         m_name;        /**< name */
   std::string         m_description; /**< brief description */
};

#endif
