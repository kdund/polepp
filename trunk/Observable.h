#include <iostream>
#include <vector>
#include <cmath>
#include "Random.h"

template <typename T>
class Pdf {
public:
  Pdf() {};
  virtual ~Pdf() {};
  //
  virtual double F(T val)=0;
  //
  double operator()(T val) { return F(val); }
};

class GaussPdf : public Pdf<double> {
public:
  GaussPdf():Pdf<double>() {m_mean=0.0; m_sigma=1.0;}
  GaussPdf(double mean, double sigma):Pdf<double>() {m_mean=mean; m_sigma=sigma;}
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
  PoissonPdf():Pdf<int>() {m_lambda=0.0;}
  PoissonPdf(double lambda):Pdf<int>() {m_lambda=lambda;}
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
  Observable() {m_pdf=0; m_rndGen=0; m_valid=false;}
  Observable(Pdf<T> *pdf, Random *rndgen) {m_pdf=pdf; m_rndGen=rndgen; m_valid=((pdf!=0)&&(rndgen!=0));}
  virtual ~Observable() {};
  //
  void setRndGen(Random *rndgen) {m_rndGen = rndgen;}
  //
  virtual T rnd()=0;
  void setPDF(Pdf<T> *pdf) { m_pdf = pdf;}
  Pdf<T> getPDF()    { return m_pdf;}
  Random getRndGen() { return m_rndGen;}
  bool valid() { return m_valid;}

  double operator()() { return rnd(); }
  double operator()(double val) { return (m_valid ? (*m_pdf)(val):0); }

protected:
  Pdf<T> *m_pdf;
  Random *m_rndGen;
  bool    m_valid;
};

class ObservableGauss : public Observable<double> {
public:
  ObservableGauss():Observable<double>() {};
  ObservableGauss(GaussPdf *pdf, Random *rndGen):Observable<double>(pdf,rndGen) {};
  ~ObservableGauss() {};
  //
  //  void setPDF(GaussPdf *pdf) {m_pdf = pdf;}
  inline double rnd() {return (m_valid ? m_rndGen->Gaus(dynamic_cast<GaussPdf *>(m_pdf)->getMean(),dynamic_cast<GaussPdf *>(m_pdf)->getSigma()):0);}
};

class ObservablePoisson : public Observable<int> {
public:
  ObservablePoisson():Observable<int>() {};
  ObservablePoisson(PoissonPdf *pdf, Random *rndGen):Observable<int>(pdf,rndGen) {};
  ~ObservablePoisson() {};
  //
  //  void setPDF(PoissonPdf *pdf) {m_pdf = pdf;}
  inline double rnd() {return (m_valid ? m_rndGen->Poisson(dynamic_cast<PoissonPdf *>(m_pdf)->getMean(),dynamic_cast<PoissonPdf *>(m_pdf)->getSigma()):0);}
};

inline double GaussPdf::phi(double mu) {
  return (1.0L/sqrt(2.0*M_PIl))*exp(-0.5L*mu*mu);
}

inline double GaussPdf::F(double val) {
  return phi((val-m_mean)/m_sigma);
}

inline double PoissonPdf::F(double val) {
  double prob;
  if(s<50.0) {
    prob = (pow(s,n)/exp(lgamma(n+1)))*exp(-s);
  } else {
    double t = double(n)-s;
    prob = (1.0L/(sqrt(2.0*M_PIl)*sqrt(s)))*exp(-0.5L*t*t/s);
  }
  return prob;
}

