#include <iostream>
#include <vector>
#include <cmath>
#include "Random.h"

/*!
  
 */
template <typename T>
class Pdf {
public:
  Pdf() {m_table=0; m_tableNpts=0;}
  Pdf(const char *name) {m_table=0; m_tableNpts=0; m_name = name;}
//   Pdf(const Pdf<T> & other) {
//     m_name = other.getName();
//     double *tptr = other.getTable();
//     m_tableNpts =  other.getTableNpts();
//     if (tptr) {
//       m_table = new double[other.getTableNpts()];
//       for (int i=0; i<other.getTableNpts(); i++) { //Copy
// 	m_table[i] = tptr[i];
//       }
//     }
  virtual ~Pdf() {}
  //
  virtual double F(T val)=0;
  //
  double operator()(T val) { return F(val); }
  //
  const char *getName() { return m_name.c_str();}
  //
  virtual void tabulate() {}
  void clearTable() {if (m_table) delete [] m_table; m_table=0; m_tableNpts=0;}
  const double *getTable() {return m_table;}
  const double *getTableNpts() {return m_tableNpts;}
  bool isTabulated() { return (m_tableNpts>0); }

protected:
  std::string m_name;
  double *m_table;
  int     m_tableNpts;
};

class PdfGauss : public Pdf<double> {
public:
  PdfGauss():Pdf<double>("Gaussian") {m_mean=0.0; m_sigma=1.0;}
  PdfGauss(double mean, double sigma):Pdf<double>("Gaussian") {m_mean=mean; m_sigma=sigma;}
  //  PdfGauss(const PdfGauss & other):Pdf<double>(other) {m_mean=other.getMean(); m_sigma=other.getSigma(); m_name=other.getName();}
  virtual ~PdfGauss() {};
  //
  void setMean(double mean)   { m_mean  = mean;}
  void setSigma(double sigma) { m_sigma = sigma;}
  double getMean()  { return m_mean; }
  double getSigma() { return m_sigma; }
  //
  inline double F(double val);
  inline double phi(double mu);

  virtual void tabulate();
  void setTableParams(int npts, double mumax) { 
    if (isTabulated()) return;
    m_tableNpts=npts;
    m_tabMuMax = mumax;
  }
private:
  double m_mean;
  double m_sigma;
  double m_tabMuMax; // tabulate N(0,1) from 0 to m_tabMuMax
  double m_tabDmu;
};

class PdfPoisson : public Pdf<int> {
public:
  PdfPoisson():Pdf<int>("Poisson") {m_lambda=0.0; m_tabLmbN=0; m_tabLmbInd=0; m_tabNmax=0;}
  PdfPoisson(double lambda):Pdf<int>("Poisson") {m_lambda=lambda;}
  //  PdfPoisson(const PdfPoisson & other):Pdf<int>(other) {m_lambda=other.getLambda(); m_name=other.getName();}
  virtual ~PdfPoisson() {};
  //
  void setLambda(double lambda) {
    m_lambda  = lambda;
    if (isTabulated()) m_tabLmbInd = static_cast<int>(m_lambda/m_tabDlmb);
  }
  double getLambda()  { return m_lambda; }
  //
  inline double F(int val);

  virtual void tabulate();
  void setTableParams(int nlmb, double lmbmax, int nn) {
    if (isTabulated()) return;
    m_tableNpts = nlmb*nn;
    m_tabLmbMax = lmbmax;
    m_tabLmbN   = nlmb;
    m_tabDlmb = m_tabLmbMax/double(m_tabLmbN);
    m_tabNmax   = nn;
    m_tabLmbInd = static_cast<int>(m_lambda/m_tabDlmb);
  }
private:
  double m_lambda;
  double m_tabLmbMax;
  int    m_tabLmbN;
  double m_tabDlmb;
  int    m_tabLmbInd;
  int    m_tabNmax; // maximum N in table
  //
  double rawPoisson(int n, double s);
};

// class GeneralPdf : public Pdf<double> {
// public:
//   GeneralPdf():Pdf<double>("General") {}
//   GeneralPdf(std::vector<double> & x, std::vector<double> & f):Pdf<double>("General") { setData(x,f); init(); }
//   //  PdfGauss(const PdfGauss & other):Pdf<double>(other) {m_mean=other.getMean(); m_sigma=other.getSigma(); m_name=other.getName();}
//   virtual ~PdfGauss() {};
//   //
//   void setData(std::vector<double> & x, std::vector<double> & f) { m_x = x; m_f = f; }
//   void init();
//   double getMean()  { return m_mean; }
//   double getSigma() { return m_sigma; }
//   //
//   inline double F(double val);
//   inline double phi(double mu);
// private:
//   std::vector<double> m_x;
//   std::vector<double> m_f;
//   double m_mean;
//   double m_sigma;
// };


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
  Observable(Pdf<T> *pdf, RND::Random *rndgen, const char *name, const char *description=0) {
    m_pdf=pdf; m_rndGen=rndgen; m_valid=((pdf!=0)&&(rndgen!=0));
    m_locked=false;
    m_lockedValue=0;
    if (name) m_name = name;
    if (description) m_description = description;
  }
  virtual ~Observable() {};
  //
  void setRndGen(RND::Random *rndgen) {m_rndGen = rndgen;}
  //
  virtual T rnd()=0;
  void setLockedValue(T val) { m_lockedValue = val;}
  void setPDF(Pdf<T> *pdf)   { m_pdf = pdf;}
  //
  T       getLockedValue()   { return m_lockedValue;}
  Pdf<T> *getPDF()           { return m_pdf;}
  RND::Random *getRndGen()        { return m_rndGen;}
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
  RND::Random *m_rndGen;
  bool    m_valid;
  bool    m_locked;
  T       m_lockedValue;
  std::string m_name;
  std::string m_description;
};

class ObservableGauss : public Observable<double> {
public:
  ObservableGauss():Observable<double>() {};
  ObservableGauss(PdfGauss *pdf, RND::Random *rndGen, const char *name, const char *desc=0):Observable<double>(pdf,rndGen,name,desc) {};
  ~ObservableGauss() {};
  //
  //  void setPDF(PdfGauss *pdf) {m_pdf = pdf;}
  inline double rnd() {return (m_valid ? m_rndGen->gauss(dynamic_cast<PdfGauss *>(m_pdf)->getMean(),dynamic_cast<PdfGauss *>(m_pdf)->getSigma()):0);}
};

class ObservablePoisson : public Observable<int> {
public:
  ObservablePoisson():Observable<int>() {};
  ObservablePoisson(PdfPoisson *pdf, RND::Random *rndGen, const char *name, const char *desc=0):Observable<int>(pdf,rndGen,name,desc) {};
  ~ObservablePoisson() {};
  //
  //  void setPDF(PdfPoisson *pdf) {m_pdf = pdf;}
  inline int rnd() {return (m_valid ? m_rndGen->poisson(dynamic_cast<PdfPoisson *>(m_pdf)->getLambda()):0);}
};

inline double PdfGauss::phi(double mu) {
  return (1.0L/sqrt(2.0*M_PIl))*exp(-0.5L*mu*mu);
}

inline double PdfGauss::F(double val) {
  double rval;
  double mu = fabs((val-m_mean)/m_sigma); // symmetric around mu0
  if (mu>m_tabMuMax) {
    rval = phi(mu)/m_sigma;
  } else {
    int sind  = static_cast<int>(mu/m_tabDmu);
    if (sind>=m_tableNpts) {
      rval = phi(mu)/m_sigma;
    } else {
      rval = m_table[sind]/m_sigma;
    }
  }
  return rval;
}

inline double PdfPoisson::F(int val) {
  double rval;
  int index=0;
  if (m_tableNpts>0) index = val+m_tabLmbInd*m_tabNmax;
  if (index<m_tableNpts) {
    rval = m_table[index];
  } else {
    rval = rawPoisson(val,m_lambda);
  }
  return rval;
}


