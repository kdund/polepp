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

template <typename T>
class Observable {
public:
  Observable() {m_pdf=0; m_rndGen=0; m_valid=false;}
  Observable(Pdf<T> *pdf, Random *rndgen) {m_pdf=pdf; m_rndGen=rndgen; m_valid=((pdf!=0)&&(rndgen!=0));}
  virtual ~Observable() {};
  //
  virtual T rnd()=0;
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

class ObservableGauss : public Observable<double> {
public:
  ObservableGauss():Observable<double>() {};
  ObservableGauss(GaussPdf *pdf, Random *rndGen):Observable<double>(pdf,rndGen) {};
  ~ObservableGauss() {};
  //
  void setPDF(GaussPdf *pdf) {m_pdf = pdf;}
  void setRndGen(Random *rndgen) {m_rndGen = rndgen;}
  inline double rnd() {return (m_valid ? m_rndGen->Gaus(dynamic_cast<GaussPdf *>(m_pdf)->getMean(),dynamic_cast<GaussPdf *>(m_pdf)->getSigma()):0);}

};

inline double GaussPdf::phi(double mu) {
  return (1.0L/sqrt(2.0*M_PIl))*exp(-0.5L*mu*mu);
}

inline double GaussPdf::F(double val) {
  return phi((val-m_mean)/m_sigma);
}

int main(int argc, char *argv[]) {
  Random rndGen;
  GaussPdf gpdf;
  ObservableGauss myObs(&gpdf,&rndGen);
  double a;
  //
  gpdf.setMean(0.0);
  gpdf.setSigma(1.0);
  std::cout << gpdf(0.0)  << std::endl;
  std::cout << gpdf(1.0)  << std::endl;
  std::cout << gpdf(-1.0) << std::endl;
  std::cout << myObs(-1.0) << std::endl;
  for (int i=0; i<1000; i++) {
    a = myObs();
  }
  return 0;
}
