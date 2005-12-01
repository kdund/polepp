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
    Base(const Base & other) {
      copy(other);
    }

    virtual ~Base() {}
    //
    Base const & operator=(Base const & rh) {
      copy(rh);
      return *this;
    }
    //
    virtual void setObservedRnd() { std::cerr << "ERROR: Calling Base::setObservedRnd()" << std::endl;}
    void setPdf(PDF::Base *pdf)   { m_pdf = pdf; m_dist = ((pdf==0) ? PDF::DIST_NONE:pdf->getDist());}
    void setPdfDist(const PDF::DISTYPE dist) { if (m_pdf==0) m_dist = dist; }
    void setPdfMean(double m)  { m_mean = m; }
    void setPdfSigma(double m) { m_sigma = m; } //if (m_pdf) m_pdf->setSigma(m); }
    void setRndGen(RND::Random *rndgen)               { m_rndGen = rndgen;}
    void setName(const char *name)               { m_name=name;}
    void setDescription(const char *description) { m_description=description;}

    virtual void validate() { m_valid = ((m_pdf!=0) && (m_rndGen!=0)); }

    virtual void dump() const { std::cout << "OBS::Base::dump(): observable dump not yet implemented" << std::endl; }

    const double getPdfMean() const              { return m_mean; } //(m_pdf ? m_pdf->getMean():0); }
    const double getPdfSigma() const             { return m_sigma; } //(m_pdf ? m_pdf->getSigma():0); }
    const PDF::DISTYPE getPdfDist() const       { return m_dist; }
    PDF::Base         *getPdf() const           { return m_pdf;}
    const std::string & getName() const          { return m_name;}
    const std::string & getDescription() const   { return m_description;}
    const bool locked() const { return m_locked;}
    const bool valid()  const { return m_valid;}
    RND::Random *getRndGen() const { return m_rndGen;}

    Base *clone() const {
      Base *obj = new Base(*this);
      return obj;
    }
    //
  protected:
    void copy(const Base & other) {
      if (this != &other) {
	m_pdf    = other.getPdf();
	m_mean   = other.getPdfMean();
	m_sigma  = other.getPdfSigma();
	m_dist   = other.getPdfDist();
	m_rndGen = other.getRndGen();
	m_valid  = other.valid();
	m_locked = other.locked();
	m_name   = "Copy of " + other.getName();
	m_description = other.getDescription();
      }
    }
    //
    PDF::Base *m_pdf;
    double m_mean;
    double m_sigma;
    PDF::DISTYPE m_dist;
    std::string m_name;
    std::string m_description;
    RND::Random *m_rndGen;
    bool m_valid;
    bool m_locked;
  };

  template <typename T>
  class BaseType: public Base {
  public:
    BaseType():Base() {}
    BaseType(const char *name,const char *desc=0):Base(name,desc) {}
    BaseType(PDF::BaseType<T> *pdf, RND::Random *rndgen, const char *name, const char *description=0):Base(name,description) {
      m_rndGen = rndgen;
      m_valid=((pdf!=0)&&(rndgen!=0));
      m_lockedValue=0;
      setPdf(pdf);
      m_observedValue=0;
      if (pdf) m_observedValue = static_cast<T>(m_mean);
    }

    BaseType(const BaseType<T> & other):Base() { copy(other);}
    virtual ~BaseType() {}
    //
    virtual T rnd() { std::cout << "ERROR::Observable - EMPTY rnd() : " << m_rndGen << std::endl; return 0; }
    //
    BaseType<T> const & operator=(BaseType<T> const & rh) {
      copy(rh);
      return *this;
    }
    const double getPdfVal(T val) { return (m_pdf ? static_cast< BaseType<T> >(*m_pdf).getVal(val,m_mean,m_sigma):0.0); }

    T       operator()()       { return ( m_locked ? m_lockedValue:rnd()); }
    double  operator()(T val)  { return getPdfVal(val); }
    //
    void setObservedRnd()          { m_observedValue = ( m_locked ? m_lockedValue:rnd());}
    void setObservedValue(T val)   { m_observedValue = val; }
    void setLockedValue(T val)     { m_lockedValue = val;}
    //
    const T getLockedValue() const   { return m_lockedValue;}
    const T getObservedValue() const { return m_observedValue; }

    BaseType<T> *clone() const {
      BaseType<T> *obj = new BaseType<T>(*this);
      return obj;
    }

    void setIntXRange(T xmin, T xmax, T step, int n=0) { m_intXRange.setRange(xmin,xmax,step,n); }
    void initInt() {
      if (m_intXRange.n()<1) return;
      int np = m_intXRange.n();
      m_intWeight.resize(np,0.0);
      m_intX.resize(np,0);
      m_intTotal = 0;
    }
    void fillInt() {
      int np = int(m_intX.size());
      double f;
      double dx = static_cast<double>(m_intXRange.step());
      m_intTotal=0;
      if (np<1) return;
      for (int i=0; i<np; i++) {
	f = getPdfVal(m_intX[i]);
	m_intWeight[i] = f*dx;
	m_intTotal += m_intWeight[i];
      }
    }
    const double getIntWeight(int i) const { return m_intWeight[i]; }
    const double getIntegral() const       { return m_intTotal; }
  protected:
    void copy(const BaseType<T> & other) {
      if (this != &other) {
	Base::copy(other);
	m_lockedValue = other.getLockedValue();
	m_observedValue = other.getObservedValue();
      }
    }
    //
    T m_lockedValue;   // value return by rnd() if rndGen is disabled
    T m_observedValue; // the observed value

    std::vector<double> m_intWeight;
    std::vector<T>      m_intX;
    Range<T>            m_intXRange;
    double              m_intTotal;
    //
  };

  //
  // TODO: Maybe Observable<T> is NOT needed, might just rename BaseType to Observable
  //
  template <typename T>
  class Observable: public BaseType<T> {
  public:
    Observable():BaseType<T>() {}
    Observable(const char *name, const char *description=0):BaseType<T>(name,description) {}
    Observable(PDF::BaseType<T> *pdf, RND::Random *rndgen, const char *name, const char *description=0):BaseType<T>(pdf,rndgen,name,description) {}
    Observable(const Observable<T> & other):BaseType<T>(other) {}

    virtual ~Observable() {};

    Observable<T> const & operator=(Observable<T> const & rh) {
      copy(rh);
      return *this;
    }

    Observable<T> *clone() const {
      Observable<T> *obj = new Observable<T>(*this);
      return obj;
    }
    void copy(const Observable<T> & other) {
      if (this != &other) {
	BaseType<T>::copy(other);
      }
    }
    //
  };

  class ObservableGauss : public Observable<double> {
  public:
    ObservableGauss():Observable<double>("gauss","Gaussian observable") {m_dist=PDF::DIST_GAUS;};
    ObservableGauss(PDF::Gauss *pdf, RND::Random *rndGen, const char *name, const char *desc=0):Observable<double>(pdf,rndGen,name,desc) {};
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
    double rnd() { return (m_valid ? m_rndGen->gauss(m_mean,m_sigma):0); }

    ObservableGauss *clone() const {
      ObservableGauss *obj = new ObservableGauss(*this);
      return obj;
    }
  };

  class ObservablePois : public Observable<int> {
  public:
    ObservablePois():Observable<int>("poisson","Poisson observable") {m_dist=PDF::DIST_POIS;};
    ObservablePois(PDF::Poisson *pdf, RND::Random *rndGen, const char *name, const char *desc=0):Observable<int>(pdf,rndGen,name,desc) {};
    ObservablePois(PDF::PoisTab *pdf, RND::Random *rndGen, const char *name, const char *desc=0):Observable<int>(pdf,rndGen,name,desc) {};
    ObservablePois(const ObservablePois & other) {
      copy(other);
    }
    virtual ~ObservablePois() {};

    ObservablePois const & operator=(ObservablePois const & rh) {
      copy(rh);
      return *this;
    }
    //
    void setPdfMean(double m)  { m_mean = m; m_sigma = (m>0 ? sqrt(m):0.0); }
    void setPdfSigma(double m) { m_mean = m*m; m_sigma=m; }
    int rnd() {return (m_valid ? m_rndGen->poisson(m_mean):0);}

    ObservablePois *clone() const {
      ObservablePois *obj = new ObservablePois(*this);
      return obj;
    }
  };

  inline Base *makeObservable(PDF::DISTYPE dist) {
    Base *obs=0;
    switch (dist) {
    case PDF::DIST_UNDEF:
    case PDF::DIST_NONE:
      obs=new Observable<double>();
      obs->setPdf(0);
      break;
    case PDF::DIST_POIS:
      obs=new ObservablePois();
      obs->setPdf(&PDF::gPoisson);
      break;
    case PDF::DIST_GAUS:
      obs=new ObservableGauss();
      obs->setPdf(&PDF::gGauss);
      break;
    case PDF::DIST_FLAT:
      //      obs=new ObservableFlat();
      std::cout << "WARNING: Not yet implemented - ObservableFlat()" << std::endl;
      break;
    case PDF::DIST_LOGN:
      //      obs=new ObservableLogN();
      std::cout << "WARNING: Not yet implemented - ObservableLogN()" << std::endl;
      break;
    default:
      std::cout << "WARNING: Unknown distribution = " << distTypeStr(dist) << std::endl;
      break;
    }
    if (obs) {
      obs->setRndGen(&RND::gRandom);
      obs->validate();
    }
    return obs;
  }
// #ifndef OBSERVABLE_CXX
//   extern Base *makeObservable(PDF::DISTYPE dist);
// #endif
};


