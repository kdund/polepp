#ifndef OBSERVABLE_H
#define OBSERVABLE_H

#include <iostream>
#include <vector>
#include <cmath>
#include "Random.h"
#include "Pdf.h"
#include "Range.h"

/*! @namespace OBS
  @ brief contains classes for Observable

  The basic idea is that an Observable should contain:
  - the observed value
  - the assumed pdf
  - the parameters of the assumed pdf
  - the apriori probability
  Note that the given pdf is a pointer to a const, i.e the pdf is NEVER modified by Observable.

  The implementation consists of two base classes and specific implementations:
  - OBS::Base            : virtual base class
  - OBS::BaseType<T>     : templated base class inheriting from OBS::Base; the <T> is either <double> or <int>
  - OBS::ObservableGauss : gaussian pdf inheriting from BaseType<double>
  - OBS::ObservableLogN  : log-normal pdf inheriting from BaseType<double>
  - OBS::ObservablePois  : poisson pdf inheriting from BaseType<int>
  - OBS::ObservableFlat  : flat pdf inheriting from BaseType<double>
  - OBS::ObservableConst : const value 'pdf' inheriting from BaseType<double>

  The user can then use the class to generate random observations using the given pdf characteristics.
  {
     ObservablePois poisObs;
     poisObs.setPdfUseMean(2.5);
     obs = poisObs();            // return a random observation
     obs = poisObs.rnd();        // same
     poisObs.setObservedRnd();   // set the observation to a random value
     poisObs.setObservedValue(3);// set the observation to a fixed value
     poisObs.setObservedValue(); // set the observation to the used mean value of the pdf
     f   = poisObs(4);          // f = Poisson(4|mu=2.5)
     
  An observation can be forced to be locked using Base::lock().
 */
namespace OBS {

  class Base {
  public:
    inline Base();
    inline Base(const char *name, const char *description=0);
    inline Base(const Base & other);
    inline virtual ~Base();
    //
    /////////////////////////////////////////////////////
    // VIRTUAL MEMBERS - some are instances of typed template BaseType<T> with T=int,double
    // TODO: Move these to separate file....
    /////////////////////////////////////////////////////

    virtual void   setObservedRnd(void)=0;
    virtual void   setObservedValue()=0;
    //
    /////////////////////////////////////////////////////
    /////////////////////////////////////////////////////

    inline void   setObservedValue( double v );
    inline void   setObservedValue( int v );

    inline void   getObservedValue( double & v ) const;
    inline void   getObservedValue( int & v ) const;
    inline double getObservedValue() const;

    // setting name, pdf and rdn generator
    inline void           setName(const char *name);
    inline void           setDescription(const char *description);
    inline void           lock();
    inline void           unlock();
    inline void           setPdf(const PDF::Base *pdf);
    inline virtual void   setPdfUseMean(double m);
    inline virtual void   setPdfUseSigma(double m);
    inline void           setRndGen(const RND::Random *rndgen);
    inline virtual void   validate();
    inline virtual double aprioriProb( double x ) const;

    // print out
    inline virtual void dump() const;

    inline const double getPdfVal(double val) const;
    // accessors
    inline const std::string & getName()        const;
    inline const std::string & getDescription() const;
    inline const double        getPdfUseMean()  const;
    inline const double        getPdfUseSigma() const;
    inline const PDF::DISTYPE  getPdfDist()     const;
    inline const PDF::Base    *getPdf()         const;
    inline const RND::Random  *getRndGen()      const;
    
    // status
    inline const bool constant() const;
    inline const bool locked()   const;
    inline const bool valid()    const;

    // variable type
    inline virtual bool isInt()    const;
    inline virtual bool isDouble() const;
    inline virtual bool isFloat()  const;

    virtual Base *clone() const = 0;

  protected:
    inline const double getObsVal() const;
    inline void copy(const Base & other);
    //
    std::string m_name;
    std::string m_description;
    //
    // PDF def + rnd gen
    //
    const PDF::Base   *m_pdf;    //! pointer to pdf function
    double             m_mean;   //! used mean value
    double             m_sigma;  //! used sigma 
    bool               m_locked; //! flag: if true, the observed value is locked
    bool               m_valid;  //! flag: true if valid (???)
    const RND::Random *m_rndGen; //! pointer to random number generator
    double             m_obsVal; //! copy of BaseType<T>::m_observedValue
  };

  template <typename T>
  class BaseType: public Base {
  public:
    inline BaseType();
    inline BaseType(const char *name,const char *desc=0);

    inline BaseType(const BaseType<T> & other);
    inline virtual ~BaseType();
    //
    inline virtual T rnd() const;
    //
    inline BaseType<T> const & operator=(BaseType<T> const & rh);
    inline BaseType<T> *clone() const;
    inline const double getPdfVal(T val);

    inline T       operator()();
    inline double  operator()(T val);
    //
    inline void setObservedRnd();
    inline void setObservedValue(T val);
    inline void setObservedValue();
    //
    inline virtual const double getObservedValue()        const;
    inline virtual const void   getObservedValue(T & val) const;

    inline virtual void dump() const;
    inline virtual bool isInt()    const;
    inline virtual bool isDouble() const;
    inline virtual bool isFloat()  const;

  protected:

    inline void setObsVal( T val);
    inline void copy(const BaseType<T> & other);
    //
    T m_observedValue; //! the observed value
    //
  };

  class ObservableGauss : public BaseType<double> {
  public:
    inline ObservableGauss();
    inline ObservableGauss(const char *name, const char *desc=0);
    inline ObservableGauss(const ObservableGauss & other);
    inline virtual ~ObservableGauss();
    //
    inline ObservableGauss const & operator=(ObservableGauss const & rh);
    //
    inline double rnd() const;

    inline ObservableGauss *clone() const;
  };

  class ObservableLogN : public BaseType<double> {
  public:
    inline ObservableLogN();
    inline ObservableLogN(const char *name, const char *desc=0);
    inline ObservableLogN(const ObservableLogN & other);
    inline virtual ~ObservableLogN();
    //
    inline ObservableLogN const & operator=(ObservableLogN const & rh);
    inline double rnd() const;
    inline ObservableLogN *clone() const;
  };

  class ObservablePois : public BaseType<int> {
  public:
    inline ObservablePois();
    inline ObservablePois(const char *name, const char *desc=0);
    inline ObservablePois(const ObservablePois & other);
    inline virtual ~ObservablePois();

    inline ObservablePois const & operator=(ObservablePois const & rh);
    inline void setExcludeZero();
    inline void setIncludeZero();
    inline bool getExcludeZeroFlag() const;
    inline double aprioriProb( double x ) const;
    inline void setPdfUseMean(double m);
    inline void setPdfUseSigma(double m);
    inline int rnd() const;

    inline ObservablePois *clone() const;

  protected:
    bool m_excludeZero;
  };

  class ObservableFlat : public BaseType<double> {
  public:
    inline ObservableFlat();
    inline ObservableFlat(const char *name, const char *desc=0);
    inline ObservableFlat(const ObservableFlat & other);

    inline virtual ~ObservableFlat();
    //
    inline ObservableFlat const & operator=(ObservableFlat const & rh);
    //

    inline void setPdfUseMean(double m);
    inline void setPdfUseSigma(double s);
    inline void setPdfRange(double xmin, double xmax);
    inline double rnd() const;

    inline ObservableFlat *clone() const;
  private:
    double m_xmin;
    double m_xmax;
  };


  class ObservableConst : public BaseType<double> {
  public:
    inline ObservableConst();
    inline ObservableConst(const char *name, const char *desc=0);
    inline ObservableConst(const ObservableConst & other);
    inline virtual ~ObservableConst();
    //
    inline ObservableConst const & operator=(ObservableConst const & rh);
    //
    inline void setPdfUseSigma(double s);

    inline double rnd() const;

    inline ObservableConst *clone() const;
  };

  // Correlated variables - UNDER DEVELOPMENT
  // The idea is to have a list of OBS::Base ptrs, one for each variable.
  // Each varable must be a gauss.
  // The correlations are given as coefficients.
  // The coefficients are stored in a vector< vector<double> > - matching the observable list.
  // coeff[i][j] = coef[j][i]
  // All coeff[i][j] are initialized to 0 unless i=j (they are 1).
  // A flag matrix flag[i][j] should contain flags for which combinations of observables
  // that have been filled.
  class Correlated {
  public:
    inline Correlated();
    inline virtual ~Correlated();
    //
    inline bool  add(Base *obs);
    inline Base *find(const char *name);
    inline bool  setCorrelation(const char *v1, const char *v2, double c);
    inline bool  setCorrelation(const Base *v1, const Base *v2, double c);
    //
    inline void rnd() const;

  private:
    std::vector<const Base *>          m_observables;
    std::vector< std::vector<double> > m_correlations;
    std::vector< std::vector<bool> >   m_corrflags;
  };

  ///////////////////////////////////
  inline Base *makeObservable(PDF::DISTYPE dist);
  inline Base *clone( const Base *bptr );
};

#include "Observable.icc"
#endif
