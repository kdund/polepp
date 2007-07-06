#include <iostream>
#include <vector>
#include <cmath>
#include "Random.h"
#include "Pdf.h"
#include "Range.h"

//class Range<int>;
//class Range<double>;
// template class Range<int>;
// template class Range<float>;
// template class Range<double>;
//
//Range<int> gRint;

// NOTE: pdf is given as a pointer and can be manipulated by the object
//
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
    inline void           setPdf(PDF::Base *pdf);
    inline void           setPdfDist(const PDF::DISTYPE dist);
    inline virtual void   setPdfMean(double m);
    inline virtual void   setPdfSigma(double m);
    inline void           setRndGen(RND::Random *rndgen);
    inline virtual void   validate();
    inline virtual double aprioriProb( double x );
    // setting integral related
    inline void setIntNpts(int n);
    inline void setIntScale(double s);
    inline void setIntXRange(double xmin, double xmax, double step, int n=0);

    inline virtual double transIntX(double x);

    inline void calcIntRange(double & low, double & high, double mean, double sigma, bool positive=true);

    inline void initInt();

    inline void initInt(double xmin, double xmax, double step, int n=0);

    //! init with zero range
    inline void initIntConst();

    //! init with default range, given by observed value, sigma and scale
    inline void initIntegral();

    inline void fillInt();

    // print out
    inline virtual void dump() const;

    inline const double getPdfVal(double val);
    // accessors
    inline const std::string & getName()        const;
    inline const std::string & getDescription() const;
    inline const double        getPdfMean()     const;
    inline const double        getPdfSigma()    const;
    inline const PDF::DISTYPE  getPdfDist()     const;
    inline PDF::Base          *getPdf()         const;
    inline RND::Random        *getRndGen()      const;
    
    // status
    inline const bool constant() const;
    inline const bool locked()   const;
    inline const bool valid()    const;

    // integral accessors
    inline const bool                 isIntFilled()       const;
    inline const double               getIntScale()       const;
    inline const std::vector<double> *getIntWeight()      const;
    inline const double               getIntWeight(int i) const;
    inline const double               getIntegral()       const;
    inline const Range<double>       *getIntXRange()      const;
    inline const double               getIntdX()          const;
    inline const std::vector<double> *getIntX()           const;
    inline const double               getIntX(int i)      const;
    inline const double               getIntXmin()        const;
    inline const double               getIntXmax()        const;
    inline const int                  getIntN()           const;

    inline virtual bool isInt()    const;
    inline virtual bool isDouble() const;
    inline virtual bool isFloat()  const;

    virtual Base *clone() const = 0;

  protected:
    inline const double  getObsVal() const;
    inline void copy(const Base & other);
    //
    std::string m_name;
    std::string m_description;
    //
    // PDF def + rnd gen
    //
    PDF::Base *m_pdf;
    double m_mean;
    double m_sigma;
    PDF::DISTYPE m_dist;
    bool m_locked;
    bool m_valid;
    RND::Random *m_rndGen;
    double m_obsVal; //! copy of BaseType<T>::m_observedValue
    //
    // Integral def.
    //
    double              m_intScale;  //! range is defined by x+-scale*sigma
    int                 m_intNpts;   //! requested number of points in integral - NOTE not == actual number of points
    double              m_intTotal;  //! sum of all f(x)dx
    bool                m_intFilled; //! true if integral is filled
    std::vector<double> m_intWeight; //! array containing the weights f(x)dx
    std::vector<double> m_intX;      //! array of x, always double although X may be from an integer PDF
    Range<double>       m_intXRange; //! range of x

  };

  template <typename T>
  class BaseType: public Base {
  public:
    inline BaseType();
    inline BaseType(const char *name,const char *desc=0);
    inline BaseType(PDF::BaseType<T> *pdf, RND::Random *rndgen, const char *name, const char *description=0);

    inline BaseType(const BaseType<T> & other);
    inline virtual ~BaseType();
    //
    inline virtual T rnd();
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
    inline ObservableGauss(PDF::Gauss *pdf, RND::Random *rndGen, const char *name, const char *desc=0);
    inline ObservableGauss(const ObservableGauss & other);
    inline virtual ~ObservableGauss();
    //
    inline ObservableGauss const & operator=(ObservableGauss const & rh);
    //
    inline double rnd();

    inline ObservableGauss *clone() const;
  };

  class ObservableLogN : public BaseType<double> {
  public:
    inline ObservableLogN();
    inline ObservableLogN(PDF::LogNormal *pdf, RND::Random *rndGen, const char *name, const char *desc=0);
    inline ObservableLogN(const ObservableLogN & other);
    inline virtual ~ObservableLogN();
    //
    inline ObservableLogN const & operator=(ObservableLogN const & rh);
    inline double rnd();
    inline ObservableLogN *clone() const;
  protected:
    inline double transIntX(double x);
  };

  class ObservablePois : public BaseType<int> {
  public:
    inline ObservablePois();
    inline ObservablePois(PDF::Poisson *pdf, RND::Random *rndGen, const char *name, const char *desc=0);
    inline ObservablePois(PDF::PoisTab *pdf, RND::Random *rndGen, const char *name, const char *desc=0);
    inline ObservablePois(const ObservablePois & other);
    inline virtual ~ObservablePois();

    inline ObservablePois const & operator=(ObservablePois const & rh);
    inline void setExcludeZero();
    inline void setIncludeZero();
    inline bool getExcludeZeroFlag() const;
    inline double aprioriProb( double x );
    inline void setPdfMean(double m);
    inline void setPdfSigma(double m);
    inline int rnd();

    inline ObservablePois *clone() const;

  protected:
    bool m_excludeZero;
  };

  class ObservableFlat : public BaseType<double> {
  public:
    inline ObservableFlat();
    inline ObservableFlat(PDF::Flat *pdf, RND::Random *rndGen, const char *name, const char *desc=0);
    inline ObservableFlat(const ObservableFlat & other);

    inline virtual ~ObservableFlat();
    //
    inline ObservableFlat const & operator=(ObservableFlat const & rh);
    //
    inline double rnd();

    inline ObservableFlat *clone() const;
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
    inline void rnd();

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
