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
    Base() {
      m_pdf = 0; m_rndGen=0; m_valid=false; m_locked=false;
      m_mean = 0; m_sigma = 0; m_dist = PDF::DIST_UNDEF;
    }
    Base(const char *name, const char *description=0) {
      m_pdf=0; m_rndGen=0; m_valid=false; m_locked=false;
      m_mean = 0; m_sigma = 0; m_dist = PDF::DIST_UNDEF;
      if (name) m_name=name;
      if (description) m_description=description;
    }
    Base(const Base & other) { copy(other);}
    //
    virtual ~Base() {}
    //
    /////////////////////////////////////////////////////
    // VIRTUAL MEMBERS - some are instances of typed template BaseType<T> with T=int,double
    // TODO: Move these to separate file....
    /////////////////////////////////////////////////////

    virtual void   setObservedRnd(void)           { std::cout << "1MUST BE OVERLOADED!" << std::endl; }
    //    virtual void   setObservedValue(double v)     { std::cout << "2MUST BE OVERLOADED!" << std::endl; }
    //    virtual void   setObservedValue(int v)        { std::cout << "3MUST BE OVERLOADED!" << std::endl; }
    virtual void   setObservedValue()             { std::cout << "4MUST BE OVERLOADED!" << std::endl; }
    //
//     virtual const void   getObservedValue(int & val)    const { std::cout << "5MUST BE OVERLOADED!" << std::endl; val=0; }
//     virtual const void   getObservedValue(double & val) const { std::cout << "6MUST BE OVERLOADED! "
//                                                                          << getName() << std::endl; val=0.0; }
    /////////////////////////////////////////////////////
    /////////////////////////////////////////////////////

    inline void   setObservedValue( double v );
    inline void   setObservedValue( int v );
    inline void   getObservedValue( double & v ) const { v = m_obsVal; }
    inline void   getObservedValue( int & v ) const;
    const  double getObservedValue() const { return m_obsVal;}

    // setting name, pdf and rdn generator
    void         setName(const char *name)               { m_name=name;}
    void         setDescription(const char *description) { m_description=description;}
    void         lock()                                  { m_locked = true; }
    void         unlock()                                { m_locked = false; }
    void         setPdf(PDF::Base *pdf)                  { m_pdf = pdf; m_dist = ((pdf==0) ? PDF::DIST_NONE:pdf->getDist());}
    void         setPdfDist(const PDF::DISTYPE dist)     { if (m_pdf==0) m_dist = dist; }
    virtual void setPdfMean(double m)                    { m_mean = m; }
    virtual void setPdfSigma(double m)                   { m_sigma = m; } //if (m_pdf) m_pdf->setSigma(m); }
    void         setRndGen(RND::Random *rndgen)          { m_rndGen = rndgen;}
    virtual void validate()                              { m_valid = ((m_pdf!=0) && (m_rndGen!=0)); }

    // setting integral related
    void setIntNpts(int n)     { m_intNpts = n; }
    void setIntScale(double s) { m_intScale = s; }
    void setIntXRange(double xmin, double xmax, double step, int n=0) {
      m_intFilled = false;
      m_intXRange.setRange(xmin,xmax,step,n);
    }

    virtual double transIntX(double x) { return x; }

    void calcIntRange(double & low, double & high, double mean, double sigma, bool positive=true)
    {
      double dx;
      const double maxp = 0.9999;
      int n=0, nprev, nlow, nhigh;
      double p=0.0;
      //  PDF::DISTYPE dist = pdf->getDist();
      //
      if (m_dist==PDF::DIST_NONE) {
        low  = mean;
        high = mean;
      } else {
        switch (m_dist) {
        case PDF::DIST_GAUS2D:
        case PDF::DIST_GAUS:
        case PDF::DIST_LOGN:
          low  = mean - m_intScale*sigma;
          high = mean + m_intScale*sigma;
          break;
        case PDF::DIST_FLAT:
          dx=sigma*1.73205081; // == sqrt(12)*0.5; ignore scale - always use full range
          low  = mean-dx;
          high = mean+dx;
          break;
        case PDF::DIST_POIS:
          nlow  = -1;
          nhigh = -1;
          // find min and max range of poisson
          // this is defined by maxp above
          // low  : max N for wich sum( p(n) ) < 1.0-maxp
          // high : min N for wich sum( p(n) ) > maxp
          while (nhigh<0) {
            nprev=n;
            p += PDF::gPoisTab.getVal( n, mean );
            //          p += pdf->getVal( n, mean );
            //p += m_pdf->getVal( n, mean );
            if ((n==0) || (p<(1.0-maxp))) nlow  = n;
            if (p>maxp)       nhigh = n;
            n++;
            if (nprev>n) { // just a STUPID test; can be done better...
              std::cerr << "Infinite loop caugh in OBS::calcIntRange() for Poisson - brutal exit" << std::endl;
              exit(-1);
            }
          }
          //
          m_intNpts = nhigh-nlow+1;
          //
          low  = static_cast<double>(nlow);
          high = static_cast<double>(nhigh);
          break;
        default: // ERROR STATE
          low  = 0;
          high = 0;
          std::cerr << "OBS::calcIntRange() -> Unknown pdf type = " << m_dist << std::endl;
          break;
        }
      }
      if (positive && (low<0)) {
        high = high-low;
        low = 0;
      }
    }

    void initInt() {
      if (m_intXRange.n()<1) return;
      int np = m_intXRange.n();
      m_intWeight.resize(np,0.0);
      m_intX.resize(np,0.0);
      //
      double x;
      for (int i=0; i<np; i++) {
        x = m_intXRange.getVal(i);
        x = transIntX(x);
        m_intX[i] = x; //! for LOGN, the Xrange is in log(x)
      }

      m_intTotal = 0;
      m_intFilled = false;
    }

    void initInt(double xmin, double xmax, double step, int n=0) {
      setIntXRange(xmin,xmax,step,n);
      initInt();
    }

    //! init with zero range
    void initIntConst() {
      initInt(m_obsVal,m_obsVal,0,1);
    }

    //! init with default range, given by observed value, sigma and scale
    void initIntegral() {
      if (constant()) { // takes care of DIST_NONE, DIST_UNDEF
	initIntConst();
      } else {
	double low;
	double high;
	bool pos=true;
	double mean, sigma;
        switch ( this->m_dist ) {
        case PDF::DIST_LOGN:
 	  mean = PDF::calcLogMean(m_obsVal,double(m_sigma));
 	  sigma = PDF::calcLogSigma(m_obsVal,double(m_sigma));
	  pos = false;
          break;
        case PDF::DIST_POIS:
          mean = m_obsVal;
          sigma = sqrt((mean>0 ? mean:0.0));
          pos = true;
          break;
        default:
	  mean = m_obsVal;
	  sigma = m_sigma;
	  pos = true;
	}
        calcIntRange(low,high, mean , sigma, pos);
	initInt(low,high,0,m_intNpts);
      }
    }

    void fillInt() {
      if (m_intFilled) return;
      int np = int(m_intX.size());
      double f;
      //      double dx = static_cast<double>(m_intXRange.step());
      double dx;
      m_intTotal=0;
      if (np<1) return;
      if (np==1) { // Dirac spike - integral and weight == 1.0
	m_intWeight[0] = 1.0;
	m_intTotal = 1.0;
      } else {
	for (int i=0; i<np-1; i++) {
	  f = getPdfVal(m_intX[i]);
	  dx = static_cast<double>(m_intX[i+1]-m_intX[i]); // might be lognormal!
	  m_intWeight[i] = f*dx;
	  m_intTotal += m_intWeight[i];
	}
      }
      m_intFilled = true;
    }

    // print out
    virtual void dump() const {
      //      std::cout << "OBS::Base::dump(): observable dump not yet implemented" << std::endl;
      std::cout << "--------------------------------------------" << std::endl;
      std::cout << "OBS::Base:dump()" << std::endl;
      std::cout << "Name : " << getName() << std::endl;
      std::cout << "Mean : " << getPdfMean() << std::endl;
      std::cout << "Sigma: " << getPdfSigma() << std::endl;
      std::cout << "Dist : " << PDF::distTypeStr(getPdfDist()) << std::endl;
      std::cout << "Obs  : " << getObsVal() << std::endl;
      std::cout << "--------------------------------------------" << std::endl;
    }

    const double getPdfVal(double val) { return (this->m_pdf ? (this->m_pdf)->getVal(val,this->m_mean,this->m_sigma):0.0); }
    // accessors
    const std::string & getName()        const { return this->m_name;}
    const std::string & getDescription() const { return this->m_description;}
    const double        getPdfMean()     const { return this->m_mean; } //(m_pdf ? m_pdf->getMean():0); }
    const double        getPdfSigma()    const { return this->m_sigma; } //(m_pdf ? m_pdf->getSigma():0); }
    const PDF::DISTYPE  getPdfDist()     const { return this->m_dist; }
    PDF::Base          *getPdf()         const { return this->m_pdf;}
    RND::Random        *getRndGen()      const { return this->m_rndGen;}
    
    // status
    const bool constant() const { return (m_locked || (m_dist==PDF::DIST_NONE) || (m_dist==PDF::DIST_UNDEF)); }
    const bool locked()   const { return m_locked;}
    const bool valid()    const { return m_valid;}

    // integral accessors
    const bool                 isIntFilled()       const { return m_intFilled; }
    const double               getIntScale()       const { return m_intScale; }
    const std::vector<double> *getIntWeight()      const { return &m_intWeight; }
    const double               getIntWeight(int i) const { return m_intWeight[i]; }
    const double               getIntegral()       const { return m_intTotal; }
    const Range<double>       *getIntXRange()      const { return &m_intXRange; } //! Integration range - log(x) if DIST_LOGN
    const double               getIntdX()          const { return m_intXRange.step(); }
    const std::vector<double> *getIntX()           const { return &m_intX; } //! Always x (even for LOGN)
    const double               getIntX(int i)      const { return m_intX[i]; }
    const double               getIntXmin()        const { return m_intX.front(); }
    const double               getIntXmax()        const { return m_intX.back(); }
    const int                  getIntN()           const { return static_cast<int>(m_intX.size()); }


    Base *clone() const { // keep it as debug
      return new Base( *this );
    }
    //
  protected:
    virtual const bool isInt()    const { return false; }
    virtual const bool isDouble() const { return false; }
    virtual const bool isFloat()  const { return false; }
    const double  getObsVal() const { return m_obsVal; } // users should not use this one - only set in BaseType<T>
    void copy(const Base & other) {
      if (this != &other) {
	m_pdf    = other.getPdf();
	m_mean   = other.getPdfMean();
	m_sigma  = other.getPdfSigma();
	m_dist   = other.getPdfDist();
        m_obsVal = other.getObsVal();
	m_rndGen = other.getRndGen();
	m_valid  = other.valid();
	m_locked = other.locked();
	m_name   = "Copy of " + other.getName();
	m_description = other.getDescription();
        // integral stuff
	m_intFilled     = other.isIntFilled();
	m_intWeight     = *(other.getIntWeight());
	m_intXRange     = *(other.getIntXRange());
	m_intX          = *(other.getIntX());
	m_intTotal      = other.getIntegral();
	m_intScale      = other.getIntScale();
	m_intNpts       = other.getIntN();
      }
    }
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
    BaseType():Base() {}
    BaseType(const char *name,const char *desc=0):Base(name,desc) {}
    BaseType(PDF::BaseType<T> *pdf, RND::Random *rndgen, const char *name, const char *description=0):Base(name,description) {
      this->m_rndGen = rndgen;
      this->m_valid=((pdf!=0)&&(rndgen!=0));
      setPdf(pdf);
      setObsVal(0);
      if (pdf) setObsVal(static_cast<T>(this->m_mean));
      m_intFilled = false;
      m_intScale = 5.0;
      m_intNpts  = 20;
    }

    BaseType(const BaseType<T> & other):Base() { copy(other);}
    virtual ~BaseType() {}
    //
    virtual T rnd() { return m_observedValue; }
    //
    BaseType<T> const & operator=(BaseType<T> const & rh) {
      copy(rh);
      return *this;
    }

    BaseType<T> *clone() const {
      BaseType<T> *obj = new BaseType<T>(*this);
      return obj;
    }

    const double getPdfVal(T val) { return (this->m_pdf ? static_cast< PDF::BaseType<T> * >(this->m_pdf)->getVal(val,this->m_mean,this->m_sigma):0.0); }

    T       operator()()       { return ( this->m_locked ? m_observedValue:rnd()); }
    double  operator()(T val)  { return getPdfVal(val); }
    //
    void setObservedRnd()          { if (!this->m_locked) setObsVal(rnd());}
    void setObservedValue(T val)   { setObsVal(val); }
    void setObservedValue()        { setObsVal(static_cast<T>(getPdfMean())); }
    //
    virtual const double getObservedValue()        const { return m_obsVal; }
    virtual const void   getObservedValue(T & val) const { val = m_observedValue; }

    virtual void dump() const {
      //      std::cout << "OBS::Base::dump(): observable dump not yet implemented" << std::endl;
      std::cout << "--------------------------------------------" << std::endl;
      std::cout << "OBS::BaseType<T>:dump()" << std::endl;
      std::cout << "Name : " << getName() << std::endl;
      std::cout << "Mean : " << getPdfMean() << std::endl;
      std::cout << "Sigma: " << getPdfSigma() << std::endl;
      std::cout << "Dist : " << PDF::distTypeStr(getPdfDist()) << std::endl;
      std::cout << "Obs  : " << this->getObservedValue() << std::endl;
      std::cout << "--------------------------------------------" << std::endl;
    }

  protected:
    virtual const bool isInt()    const { return false; }
    virtual const bool isDouble() const { return false; }
    virtual const bool isFloat()  const { return false; }

    void setObsVal( T val) {
      m_obsVal = static_cast<double>(val);
      m_observedValue = val;
    }
    void copy(const BaseType<T> & other) {
      if (this != &other) {
	Base::copy(other);
	other.getObservedValue(m_observedValue);
      }
    }
    //
    T m_observedValue; //! the observed value
    //
  };

  template<> inline const bool BaseType<int>::isInt()       const { return true; }
  template<> inline const bool BaseType<double>::isDouble()    const { return true; }
  template<> inline const bool BaseType<float>::isFloat()     const { return true; }

  //
  // TODO: Maybe Observable<T> is NOT needed, might just rename BaseType to Observable
  //
//   template <typename T>
//   class Observable: public BaseType<T> {
//   public:
//     Observable():BaseType<T>() {}
//     Observable(const char *name, const char *description=0):BaseType<T>(name,description) {}
//     Observable(PDF::BaseType<T> *pdf, RND::Random *rndgen, const char *name, const char *description=0):BaseType<T>(pdf,rndgen,name,description) {}
//     Observable(const Observable<T> & other):BaseType<T>(other) {}

//     virtual ~Observable() {};

//     Observable<T> const & operator=(Observable<T> const & rh) {
//       copy(rh);
//       return *this;
//     }

//     Observable<T> *clone() const {
//       Observable<T> *obj = new Observable<T>(*this);
//       return obj;
//     }
//     void copy(const Observable<T> & other) {
//       if (this != &other) {
// 	BaseType<T>::copy(other);
//       }
//     }
//     //
//   };

  class ObservableGauss : public BaseType<double> {
  public:
    ObservableGauss():BaseType<double>("gauss","Gaussian observable") {this->m_dist=PDF::DIST_GAUS;};
    ObservableGauss(PDF::Gauss *pdf, RND::Random *rndGen, const char *name, const char *desc=0):BaseType<double>(pdf,rndGen,name,desc) {};
    ObservableGauss(const ObservableGauss & other) {
      copy(other);
    }
    virtual ~ObservableGauss() {};
    //
    ObservableGauss const & operator=(ObservableGauss const & rh) {
      copy(rh);
      return *this;
    }
    //
    double rnd() { return (this->m_valid ? this->m_rndGen->gauss(this->m_mean,this->m_sigma):0); }

    ObservableGauss *clone() const {
      ObservableGauss *obj = new ObservableGauss(*this);
      return obj;
    }
  };

  class ObservableLogN : public BaseType<double> {
  public:
    ObservableLogN():BaseType<double>("gauss","Gaussian observable") {m_dist=PDF::DIST_LOGN;};
    ObservableLogN(PDF::LogNormal *pdf, RND::Random *rndGen, const char *name, const char *desc=0):BaseType<double>(pdf,rndGen,name,desc) {};
    ObservableLogN(const ObservableLogN & other) {
      copy(other);
    }
    virtual ~ObservableLogN() {};
    //
    ObservableLogN const & operator=(ObservableLogN const & rh) {
      copy(rh);
      return *this;
    }
    //
    double rnd() {
      PDF::LogNormal *pdf = static_cast<PDF::LogNormal *>(m_pdf);
      return (m_valid ? m_rndGen->logNormalLN(pdf->getLogMean(),pdf->getLogSigma()):0);
    }

    ObservableLogN *clone() const {
      ObservableLogN *obj = new ObservableLogN(*this);
      return obj;
    }
  protected:
    double transIntX(double x) { return exp(x); }
  };

  class ObservablePois : public BaseType<int> {
  public:
    ObservablePois():BaseType<int>("poisson","Poisson observable") {this->m_dist=PDF::DIST_POIS;};
    ObservablePois(PDF::Poisson *pdf, RND::Random *rndGen, const char *name, const char *desc=0):BaseType<int>(pdf,rndGen,name,desc) {};
    ObservablePois(PDF::PoisTab *pdf, RND::Random *rndGen, const char *name, const char *desc=0):BaseType<int>(pdf,rndGen,name,desc) {};
    ObservablePois(const ObservablePois & other) {
      copy(other);
    }
    virtual ~ObservablePois() {};

    ObservablePois const & operator=(ObservablePois const & rh) {
      copy(rh);
      return *this;
    }
    //
    void setPdfMean(double m)  { this->m_mean = m;   this->m_sigma = (m>0 ? sqrt(m):0.0); }
    void setPdfSigma(double m) { this->m_mean = m*m; this->m_sigma=m; }
    int rnd() {return (this->m_valid ? this->m_rndGen->poisson(this->m_mean):0);}

    ObservablePois *clone() const {
      ObservablePois *obj = new ObservablePois(*this);
      return obj;
    }
  };

  class ObservableFlat : public BaseType<double> {
  public:
    ObservableFlat():BaseType<double>("flat","Flat observable") {m_dist=PDF::DIST_FLAT;};
    ObservableFlat(PDF::Flat *pdf, RND::Random *rndGen, const char *name, const char *desc=0):BaseType<double>(pdf,rndGen,name,desc) {};
    ObservableFlat(const ObservableFlat & other) {
      copy(other);
    }
    virtual ~ObservableFlat() {};
    //
    ObservableFlat const & operator=(ObservableFlat const & rh) {
      copy(rh);
      return *this;
    }
    //
    double rnd() { return (m_valid ? m_rndGen->flat(m_mean,m_sigma):0); }

    ObservableFlat *clone() const {
      ObservableFlat *obj = new ObservableFlat(*this);
      return obj;
    }
  };

  inline Base *makeObservable(PDF::DISTYPE dist) {
    Base *obs=0;
    switch (dist) {
    case PDF::DIST_UNDEF:
    case PDF::DIST_NONE:
      obs=new BaseType<double>();
      obs->setPdf(0);
      break;
    case PDF::DIST_POIS:
      obs=new ObservablePois();
      obs->setPdf(&PDF::gPoisTab);
      break;
    case PDF::DIST_GAUS:
      obs=new ObservableGauss();
      obs->setPdf(&PDF::gGauss);
      break;
    case PDF::DIST_FLAT:
      obs=new ObservableFlat();
      obs->setPdf(&PDF::gFlat);
      break;
    case PDF::DIST_LOGN:
      obs=new ObservableLogN();
      obs->setPdf(&PDF::gLogNormal);
      break;
    default:
      std::cout << "FATAL: Unknown distribution = " << distTypeStr(dist) << std::endl;
      exit(-1);
      break;
    }
    if (obs) {
      obs->setRndGen(&RND::gRandom);
      obs->validate();
    }
    return obs;
  };


  inline void Base::setObservedValue( double v ) {
    if (this->isDouble())   static_cast<BaseType<double> *>(this)->setObservedValue(v);
    else if (this->isInt()) static_cast<BaseType<int> *>(this)->setObservedValue(static_cast<int>(v));
    else {
      std::cout << "FATAL: OBS::Base::setObservedValue(double): not supported type!" << std::endl;
      exit(-1);
    }
  }
  inline void Base::setObservedValue( int v ) {
    if (this->isInt())         static_cast<BaseType<int> *>(this)->setObservedValue(v);
    else if (this->isDouble()) static_cast<BaseType<double> *>(this)->setObservedValue(static_cast<int>(v));
    else {
      std::cout << "FATAL: OBS::Base::setObservedValue(int): not supported type!" << std::endl;
      exit(-1);
    }
  }
  inline void Base::getObservedValue( int & v ) const {
    if (this->isInt())         static_cast<const BaseType<int> *>(this)->getObservedValue(v);
    else if (this->isDouble()) v = static_cast<const int>(m_obsVal);
    else {
      std::cout << "FATAL: OBS::Base::getObservedValue(int): not supported type!" << std::endl;
      exit(-1);
    }
  }

// #ifndef OBSERVABLE_CXX
//   extern Base *makeObservable(PDF::DISTYPE dist);
// #endif
};


