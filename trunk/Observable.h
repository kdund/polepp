#include <iostream>
#include <vector>
#include <cmath>
#include "Random.h"

/*!
  
 */
template <typename T>
class Pdf {
public:
  Pdf() {}
  Pdf(const char *name) {m_name = name;}
  Pdf(const Pdf<T> & other) {m_name = other.getName();}
  virtual ~Pdf() {}
  //
  virtual double F(T val)=0;
  //
  double operator()(T val) { return F(val); }
  //
  const char *getName() { return m_name.c_str();}
private:
  std::string m_name;
};

class GaussPdf : public Pdf<double> {
public:
  GaussPdf():Pdf<double>("Gaussian") {m_mean=0.0; m_sigma=1.0;}
  GaussPdf(double mean, double sigma):Pdf<double>("Gaussian") {m_mean=mean; m_sigma=sigma;}
  //  GaussPdf(const GaussPdf & other):Pdf<double>(other) {m_mean=other.getMean(); m_sigma=other.getSigma(); m_name=other.getName();}
  virtual ~GaussPdf() {};
  //
  void setMean(double mean)   { m_mean  = mean;}
  void setSigma(double sigma) { m_sigma = sigma;}
  double getMean()  { return m_mean; }
  double getSigma() { return m_sigma; }
  //
  inline double F(double val);
  inline double phi(double mu);
private:
  double m_mean;
  double m_sigma;
};

class PoissonPdf : public Pdf<int> {
public:
  PoissonPdf():Pdf<int>("Poisson") {m_lambda=0.0;}
  PoissonPdf(double lambda):Pdf<int>("Poisson") {m_lambda=lambda;}
  //  PoissonPdf(const PoissonPdf & other):Pdf<int>(other) {m_lambda=other.getLambda(); m_name=other.getName();}
  virtual ~PoissonPdf() {};
  //
  void setLambda(double lambda)   { m_lambda  = lambda;}
  double getLambda()  { return m_lambda; }
  //
  inline double F(int val);
private:
  double m_lambda;
};


template <typename T>
class Observable {
public:
  Observable() {
    m_pdf=0; m_rndGen=0; m_valid=false; m_locked=false; m_lockedValue=0;
  }
  Observable(const char *name, const char *description=0) {
    m_pdf=0; m_rndGen=0; m_valid=false; m_locked=false;
    m_lockedValue=0;
    if (name) m_name=name;
    if (description) m_description=description;
  }
  Observable(Pdf<T> *pdf, Random *rndgen, const char *name, const char *description=0) {
    m_pdf=pdf; m_rndGen=rndgen; m_valid=((pdf!=0)&&(rndgen!=0));
    m_locked=false;
    m_lockedValue=0;
    if (name) m_name = name;
    if (description) m_description = description;
  }
  virtual ~Observable() {};
  //
  void setRndGen(Random *rndgen) {m_rndGen = rndgen;}
  //
  virtual T rnd()=0;
  void setLockedValue(T val) { m_lockedValue = val;}
  void setPDF(Pdf<T> *pdf)   { m_pdf = pdf;}
  //
  T       getLockedValue()   { return m_lockedValue;}
  Pdf<T> *getPDF()           { return m_pdf;}
  Random *getRndGen()        { return m_rndGen;}
  bool    valid()            { return m_valid;}

  T       operator()()       { return (m_locked ? m_lockedValue:rnd()); }
  double  operator()(T val)  { return (m_valid ? (*m_pdf)(val):0); }

  void setName(const char *name)               { m_name=name;}
  void setDescription(const char *description) { m_description=description;}
  const char *getName()                        { return m_name.c_str();}
  const char *getDescription()                 { return m_description.c_str();}
  //
protected:
  Pdf<T> *m_pdf;
  Random *m_rndGen;
  bool    m_valid;
  bool    m_locked;
  T       m_lockedValue;
  std::string m_name;
  std::string m_description;
};

class ObservableGauss : public Observable<double> {
public:
  ObservableGauss():Observable<double>() {};
  ObservableGauss(GaussPdf *pdf, Random *rndGen, const char *name, const char *desc=0):Observable<double>(pdf,rndGen,name,desc) {};
  ~ObservableGauss() {};
  //
  //  void setPDF(GaussPdf *pdf) {m_pdf = pdf;}
  inline double rnd() {return (m_valid ? m_rndGen->Gaus(dynamic_cast<GaussPdf *>(m_pdf)->getMean(),dynamic_cast<GaussPdf *>(m_pdf)->getSigma()):0);}
};

class ObservablePoisson : public Observable<int> {
public:
  ObservablePoisson():Observable<int>() {};
  ObservablePoisson(PoissonPdf *pdf, Random *rndGen, const char *name, const char *desc=0):Observable<int>(pdf,rndGen,name,desc) {};
  ~ObservablePoisson() {};
  //
  //  void setPDF(PoissonPdf *pdf) {m_pdf = pdf;}
  inline int rnd() {return (m_valid ? m_rndGen->Poisson(dynamic_cast<PoissonPdf *>(m_pdf)->getLambda()):0);}
};

inline double GaussPdf::phi(double mu) {
  return (1.0L/sqrt(2.0*M_PIl))*exp(-0.5L*mu*mu);
}

inline double GaussPdf::F(double val) {
  return phi((val-m_mean)/m_sigma);
}

inline double PoissonPdf::F(int val) {
  double prob;
  if(m_lambda<50.0) {
    prob = (pow(m_lambda,val)/exp(lgamma(val+1)))*exp(-m_lambda);
  } else {
    double t = double(val)-m_lambda;
    prob = (1.0L/(sqrt(2.0*M_PIl)*sqrt(m_lambda)))*exp(-0.5L*t*t/m_lambda);
  }
  return prob;
}

