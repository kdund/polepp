#ifndef MEASUREMENT_H
#define  MEASUREMENT_H
//
// One measurement consists of N(observed) plus a number of nuisance parameters.
// The signal is a function of these parameters and N. It is implemented in getSignal().
//
// The idea with this is:
// A measurement consists of a N(obs) + a number of nuisance parameters
//
#include <string>
#include <list>
#include <vector>
#include "Observable.h"

//
// Note virtual functions:
//  getM(s) : m = f(s,nuisance) ; i.e the relation between measured value, signal and nuisance params
//  getSignal(m,nuisance) ...
//
template <typename T>
class Measurement {
 public:
  Measurement() { m_observable=0; }
  Measurement(const char *name, const char *desc=0) { m_observable=0; m_name = name; m_description = desc;}
  Measurement(const Measurement<T> & m) { copy(m); }

  virtual ~Measurement() { if (m_observable) delete m_observable; deleteNuisance(); }
  //
  inline Measurement const & operator=( Measurement<T> const & m ) {
    copy(m);
    return *this;
  }
  inline bool operator==( const Measurement<T> & m ) {
    bool rval = false;
    return rval;
  }

  void setTrueSignal(double s)               { m_trueSignal = s; }
  void setObservable(OBS::BaseType<T> * obs) { if (m_observable) delete m_observable; m_observable = (obs ? obs->clone():0); }
  void setName(const char *name)             { m_name = name; }
  void setDescription(const char *descr)     { m_description = descr; }
  //
  OBS::Base *addNuisance(OBS::Base * nuPar) { if (nuPar) m_nuisancePars.push_back( nuPar ); return nuPar;}
  void deleteNuisance() {
    for (std::list< OBS::Base * >::iterator it = m_nuisancePars.begin();
	 it != m_nuisancePars.end();
	 ++it) {
      if (*it) delete *it;
    }
    m_nuisancePars.clear();
  }
  void dump() const {
    std::cout << "-----------MEASUREMENT------------------------\n";
    if (m_observable) m_observable->dump();
    std::cout << "----------------------------------------------\n";
  }

  void copyNuisance(std::list< OBS::Base * > & newList ) const {
    newList.clear();
    if (m_nuisancePars.size()>0) {
      for (std::list< OBS::Base * >::const_iterator it = m_nuisancePars.begin();
	   it !=  m_nuisancePars.end();
	   ++it) {
	newList.push_back((*it)->clone());
      }
    }
  }
  //
  const double getTrueSignal() const        { return m_trueSignal; }
  const T getObsVal() const                 { return (m_observable ? m_observable->getObservedValue():0); }
  const double getObsPdfMean() const        { return (m_observable ? m_observable->getPdfMean():0); }
  const double getObsPdfSigma() const       { return (m_observable ? m_observable->getPdfSigma():0); }
  const PDF::DISTYPE getObsPdfDist() const  { return (m_observable ? m_observable->getPdfDist():0); }
  const OBS::BaseType<T> *getObservable() const   { return m_observable; }
  //
  const std::string & getName()                      const { return m_name;}
  const std::string & getDescription()               const { return m_description;}
  const std::list< OBS::Base * > & getNuisanceList() const { return m_nuisancePars; }

  const double rndObs() { OBS::BaseType<T> *p = static_cast< OBS::BaseType<T> * >(m_observable); return (*p)(); }
  //
  virtual const T getM(double s)=0;
  virtual const double getSignal()=0;
  virtual const double getSignalUnc()=0; 

  bool generateObservation() { // generates a random pseudoexperiment and stores this in the observed parts of the observables
    //
    // First set observed nuisance to random values
    //
    for (std::list< OBS::Base * >::iterator it = m_nuisancePars.begin();
	 it !=  m_nuisancePars.end();
       ++it) {
      //      obs = static_cast< OBS::Base * >(*it);
      (*it)->setObservedRnd(); // set to a random value according to pdf
    }
    //
    // get the new average to be used
    //
    double mean = getM(m_trueSignal);
    m_observable->setPdfMean(mean);
    //
    // ...and then the observable
    //
    OBS::BaseType<T> *obs = static_cast< OBS::BaseType<T> * >(m_observable);
    obs->setObservedRnd();
    //
    return true;
  }
  //
 protected:
  void copy(const Measurement<T> & other) {
    if (this != &other) {
      m_name        = other.getName();
      m_description = other.getDescription();
      m_trueSignal  = other.getTrueSignal();
      if (m_observable) delete m_observable;
      m_observable  = (other.getObservable())->clone(); // make a clone - note PDF object is NOT cloned... pointer retained (speed/mem issues)
      other.copyNuisance(m_nuisancePars);// idem
    }
  }
  OBS::Base *makeNuisance(OBS::Base *np, PDF::DISTYPE dist) {
    if (np==0) {
      np = OBS::makeObservable(dist);
      addNuisance(np);
    }
    return np;
  }

  std::string m_name;
  std::string m_description;
  //

  double                                  m_trueSignal;
  OBS::BaseType<T>                       *m_observable;
  std::list< OBS::Base * >                m_nuisancePars;
  std::vector< std::vector<double> >      m_corrMat;
};

class MeasPois : public Measurement<int> {
 public:
  MeasPois() : Measurement<int>() {};
  MeasPois(const char *name, const char *desc=0) : Measurement<int>(name,desc) {};
  MeasPois(OBS::ObservablePois * obs) : Measurement<int>(obs->getName().c_str(),obs->getDescription().c_str()) {
    m_observable = obs->clone();
  }
  MeasPois(const MeasPois & other):Measurement<int>() { copy(other);}
  virtual ~MeasPois() {}
  //
  void copy(const MeasPois & other) { Measurement<int>::copy(other);}
  void setObservable(OBS::ObservablePois * obs) { if (m_observable) delete m_observable; m_observable = (obs ? obs->clone():0); }

  virtual const int    getM(double s) { return 0;}
  virtual const double getSignal()    { return 0;}
  virtual const double getSignalUnc() { return 0;}

};

class MeasPoisEB : public MeasPois {
 public:
  MeasPoisEB() : MeasPois() { m_eff=0; m_bkg=0; }
  MeasPoisEB(const char *name, const char *desc=0) : MeasPois(name,desc) { m_eff=0; m_bkg=0; }
  MeasPoisEB(const MeasPoisEB & other):MeasPois() {copy(other);}
  virtual ~MeasPoisEB() { }
  //
  void copy(const MeasPoisEB & other) {
    if (this != &other) {
      MeasPois::copy(other);
      const OBS::BaseType<double> *eff, *bkg;
      m_eff=0;
      m_bkg=0;
      eff = other.getEff();
      if (eff) m_eff = eff->clone();
      bkg = other.getBkg();
      if (bkg) m_bkg = bkg->clone();
    }
  }

  void setEffObs(double eff) {
    if (m_eff) {
      m_eff->setObservedValue(eff);
    }
  }
  void setBkgObs(double bkg) {
    if (m_bkg) {
      m_bkg->setObservedValue(bkg);
    }
  }

  void setEffPdf(double eff, double sigma, PDF::DISTYPE dist) {
    if (m_eff==0) m_eff = static_cast< OBS::BaseType<double> *>(makeNuisance(m_eff,dist));
    m_eff->setPdfMean(eff);
    m_eff->setPdfSigma(sigma);
  }

  void setBkgPdf(double bkg, double sigma, PDF::DISTYPE dist) {
    if (m_bkg==0) m_bkg = static_cast< OBS::BaseType<double> *>(makeNuisance(m_bkg,dist));
    m_bkg->setPdfMean(bkg);
    m_bkg->setPdfSigma(sigma);
  }

  void setEffPdfMean(double m)  { if (m_eff) m_eff->setPdfMean(m); }
  void setEffPdfSigma(double m) { if (m_eff) m_eff->setPdfSigma(m); }

  void setBkgPdfMean(double m)  { if (m_bkg) m_bkg->setPdfMean(m); }
  void setBkgPdfSigma(double m) { if (m_bkg) m_bkg->setPdfSigma(m); }

  const double getEffObs() const { return (m_eff ? m_eff->getObservedValue():0); }
  const double getBkgObs() const { return (m_bkg ? m_bkg->getObservedValue():0); }

  const OBS::BaseType<double> *getEff() const {
    return m_eff;
  }
  const double getEffPdfMean() const {
    return (m_eff ? m_eff->getPdfMean():0);
  }
  const double getEffPdfSigma() const {
    return (m_eff ? m_eff->getPdfSigma():0);
  }
  const PDF::DISTYPE getEffPdfDist() const {
    return (m_eff ? m_eff->getPdfDist():PDF::DIST_UNDEF);
  }

  const OBS::BaseType<double> *getBkg() const {
    return m_bkg;
  }
  const double getBkgPdfMean() const {
    return (m_bkg ? m_bkg->getPdfMean():0);
  }
  const double getBkgPdfSigma() const {
    return (m_bkg ? m_bkg->getPdfSigma():0);
  }
  const PDF::DISTYPE getBkgPdfDist() const {
    return (m_bkg ? m_bkg->getPdfDist():PDF::DIST_UNDEF);
  }

  virtual const int    getM(double s) const {
    return static_cast<int>(m_eff->getObservedValue()*s + m_bkg->getObservedValue());
  }

  virtual const double getSignal() const {
    double dn = m_observable->getObservedValue() - m_bkg->getObservedValue();
    double e = m_eff->getObservedValue(); // >0
    return (e>0 ? dn/e : 0);
  }


 private:
  OBS::BaseType<double> *m_eff; // pointers to nuisance params in list
  OBS::BaseType<double> *m_bkg;
};

#endif
