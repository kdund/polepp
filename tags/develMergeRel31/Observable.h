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
  template <typename T>
  void getIntRange(T & low, T & high, double scale, double mean, double sigma, PDF::DISTYPE dist, bool positive=true) {
    //
    double dx;
    //
    if (dist==PDF::DIST_NONE) {
      low  = static_cast<T>(mean);
      high = static_cast<T>(mean);
    } else {
      switch (dist) {
      case PDF::DIST_GAUS2D:
      case PDF::DIST_GAUS:
      case PDF::DIST_LOGN:
	low  = static_cast<T>(mean - scale*sigma);
	high = static_cast<T>(mean + scale*sigma);
	break;
      case PDF::DIST_FLAT:
	dx=sigma*1.73205081; // == sqrt(12)*0.5; ignore scale - always use full range
	low  = static_cast<T>(mean-dx);
	high = static_cast<T>(mean+dx);
	break;
      case PDF::DIST_POIS:
	
      default: // ERROR STATE
	std::cerr << "OBS::getIntRange() -> Unknown pdf type = " << dist << std::endl;
	break;
      }
    }
    if (positive && (low<0)) {
      high = high-low;
      low = 0;
    }
  }

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
    void lock() { m_locked = true; }
    void unlock() { m_locked = false; }
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
    const bool constant() const { return (m_locked || (m_dist==PDF::DIST_NONE) || (m_dist==PDF::DIST_UNDEF)); }
    const bool locked() const { return m_locked;}
    const bool valid()  const { return m_valid;}
    RND::Random *getRndGen() const { return m_rndGen;}

    Base *clone() const {
      Base *obj = new Base(*this);
      return obj;
    }

    virtual void initInt()        {std::cerr << "ERROR: using Base::initInt()" << std::endl;}
    virtual void initIntConst()   {std::cerr << "ERROR: using Base::initIntConst()" << std::endl;}
    virtual void initIntDefault() {std::cerr << "ERROR: using Base::initIntDefault()" << std::endl;}
    virtual void fillInt()        {std::cerr << "ERROR: using Base::fillInt()" << std::endl;}
    //
    virtual const bool isIntFilled() const {std::cerr << "ERROR: using Base::isIntFilled()" << std::endl; return false;}
    virtual const int  getIntNpts()  const {std::cerr << "ERROR: using Base::getIntNpts()" << std::endl; return 0;}
    virtual const int  getIntN()     const {std::cerr << "ERROR: using Base::getIntN()" << std::endl; return 0;}
    virtual const double getIntdX()  const {std::cerr << "ERROR: using Base::getIntdX()" << std::endl; return 0;}
    virtual const double getIntWeight(int i) const {std::cerr << "ERROR: using Base::getIntWeight()" << std::endl; return 1.0;}
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
      this->m_rndGen = rndgen;
      this->m_valid=((pdf!=0)&&(rndgen!=0));
      setPdf(pdf);
      m_observedValue=0;
      if (pdf) m_observedValue = static_cast<T>(this->m_mean);
      m_intFilled = false;
      m_intScale = 5.0;
      m_intNpts  = 20;
    }

    BaseType(const BaseType<T> & other):Base() { copy(other);}
    virtual ~BaseType() {}
    //
    virtual T rnd() { std::cout << "ERROR::Observable - EMPTY rnd() : " << this->m_rndGen << std::endl; return 0; }
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
    void setObservedRnd()          { if (!this->m_locked) m_observedValue = rnd();}
    void setObservedValue(T val)   { m_observedValue = val; }
    //
    const T getObservedValue() const { return m_observedValue; }

    const bool isIntFilled()   const { return m_intFilled; }
    const int getIntNpts()     const { std::cerr << "WARNING:: Don't use Observable::getIntNpts() - use getIntN()" << std::endl; return m_intNpts; } // USE getIntN()
    const double getIntScale() const { return m_intScale; }
    const std::vector<double> * getIntWeight() const { return &m_intWeight; }
    const double getIntWeight(int i) const { return m_intWeight[i]; }
    const double getIntegral()       const { return m_intTotal; }
    const Range<T> * getIntXRange()  const { return &m_intXRange; }
    const std::vector<T> * getIntX() const { return &m_intX; }
    const T getIntX(int i)           const { return m_intX[i]; }
    const T getIntXmin()             const { return m_intX.front(); }
    const T getIntXmax()             const { return m_intX.back(); }
    const T getIntStep()             const { return m_intXRange.step(); }
    const double getIntdX()          const { return static_cast<double>(m_intXRange.step()); }
    const int getIntN()              const { return m_intXRange.n(); }

    void setIntNpts(int n)     { m_intNpts = n; }
    void setIntScale(double s) { m_intScale = s; }
    void setIntXRange(T xmin, T xmax, T step, int n=0) {
      m_intFilled = false;
      m_intXRange.setRange(xmin,xmax,step,n);
    }

    void initInt() {
      if (m_intXRange.n()<1) return;
      int np = m_intXRange.n();
      m_intWeight.resize(np,0.0);
      m_intX.resize(np,0);
      for (int i=0; i<np; i++) {
	m_intX[i] = m_intXRange.getVal(i);
      }
      m_intTotal = 0;
      m_intFilled = false;
    }

    void initInt(T xmin, T xmax, T step, int n=0) {
      setIntXRange(xmin,xmax,step,n);
      initInt();
    }

    //! init with zero range
    void initIntConst() {
      initInt(m_observedValue,m_observedValue,0,1);
    }

    //! init with default range, given by observed value, sigma and scale
    void initIntDefault() {
      if (constant()) { // takes care of DIST_NONE, DIST_UNDEF
	std::cout << "initIntDefault:: constant()" << std::endl;
	initIntConst();
      } else {
	T low;
	T high;
	bool pos;
	double mean, sigma;
	if (this->m_dist == PDF::DIST_LOGN) {
	  mean = PDF::calcLogMean(double(m_observedValue),m_sigma);
	  sigma = PDF::calcLogSigma(double(m_observedValue),m_sigma);
	  pos = false;
	} else {
	  mean = double(m_observedValue);
	  sigma = m_sigma;
	  pos = true;
	}
	getIntRange(low,high,m_intScale, mean , sigma, this->m_dist, pos);
	std::cout << "Int dist  = " << mean << " , " << sigma << std::endl;
	std::cout << "    range = " << low << " : " << high << std::endl;
	initInt(low,high,0,m_intNpts);
      }
    }

    void fillInt() {
      if (m_intFilled) return;
      int np = int(m_intX.size());
      double f;
      double dx = static_cast<double>(m_intXRange.step());
      m_intTotal=0;
      if (np<1) return;
      if (np==1) { // Dirac spike - integral and weight == 1.0
	m_intWeight[0] = 1.0;
	m_intTotal = 1.0;
      } else {
	for (int i=0; i<np; i++) {
	  f = getPdfVal(m_intX[i]);
	  m_intWeight[i] = f*dx;
	  m_intTotal += m_intWeight[i];
	}
      }
      m_intFilled = true;
    }
  protected:
    void copy(const BaseType<T> & other) {
      if (this != &other) {
	Base::copy(other);
	m_observedValue = other.getObservedValue();
	m_intFilled     = other.isIntFilled();
	m_intWeight     = *(other.getIntWeight());
	m_intXRange     = *(other.getIntXRange());
	m_intX          = *(other.getIntX());
	m_intTotal      = other.getIntegral();
	m_intScale      = other.getIntScale();
	m_intNpts       = other.getIntNpts();
      }
    }
    //
    T m_observedValue; //! the observed value

    std::vector<double> m_intWeight; //! array containing the weights f(x)dx
    std::vector<T>      m_intX;      //! array of x
    Range<T>            m_intXRange; //! range of x
    double              m_intScale;  //! range is defined by x+-scale*sigma
    int                 m_intNpts;   //! set number of points in integral
    double              m_intTotal;  //! sum of all f(x)dx
    bool                m_intFilled; //! true if integral is filled
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
    ObservableGauss():Observable<double>("gauss","Gaussian observable") {this->m_dist=PDF::DIST_GAUS;};
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
    double rnd() { return (this->m_valid ? this->m_rndGen->gauss(this->m_mean,this->m_sigma):0); }

    ObservableGauss *clone() const {
      ObservableGauss *obj = new ObservableGauss(*this);
      return obj;
    }
  };

  class ObservableLogN : public Observable<double> {
  public:
    ObservableLogN():Observable<double>("gauss","Gaussian observable") {m_dist=PDF::DIST_LOGN;};
    ObservableLogN(PDF::LogNormal *pdf, RND::Random *rndGen, const char *name, const char *desc=0):Observable<double>(pdf,rndGen,name,desc) {};
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
  };

  class ObservablePois : public Observable<int> {
  public:
    ObservablePois():Observable<int>("poisson","Poisson observable") {this->m_dist=PDF::DIST_POIS;};
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
    void setPdfMean(double m)  { this->m_mean = m;   this->m_sigma = (m>0 ? sqrt(m):0.0); }
    void setPdfSigma(double m) { this->m_mean = m*m; this->m_sigma=m; }
    int rnd() {return (this->m_valid ? this->m_rndGen->poisson(this->m_mean):0);}

    ObservablePois *clone() const {
      ObservablePois *obj = new ObservablePois(*this);
      return obj;
    }
  };

  class ObservableFlat : public Observable<double> {
  public:
    ObservableFlat():Observable<double>("flat","Flat observable") {m_dist=PDF::DIST_FLAT;};
    ObservableFlat(PDF::Flat *pdf, RND::Random *rndGen, const char *name, const char *desc=0):Observable<double>(pdf,rndGen,name,desc) {};
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
      obs=new ObservableFlat();
      break;
    case PDF::DIST_LOGN:
      obs=new ObservableLogN();
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


