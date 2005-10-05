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
#include "ObsNew.h"

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

  void copyNuisanceList(std::list< OBS::Base * > & newList ) const {
    newList.clear();
    if (m_nuisancePars.size()>0) {
      for (std::list< OBS::Base * >::iterator it = m_nuisancePars.begin();
	   it !=  m_nuisancePars.end();
	   ++it) {
	(*it)->rnd(); // save copies
      }
    }
  }
  //
  const double getTrueSignal() const        { return m_trueSignal; }
  const T getObsVal() const                 { return (m_observable ? m_observable->getObservedValue():0); }
  const double getObsPdfMean() const        { return (m_observable ? m_observable->getPdfMean():0); }
  const double getObsPdfSigma() const       { return (m_observable ? m_observable->getPdfSigma():0); }
  const PDFN::DISTYPE getObsPdfDist() const { return (m_observable ? m_observable->getPdfDist():0); }
  //
  OBS::BaseType<T> *getObservable() const   { return m_observable; }
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
  void copy(const Measurement & other) {
    if (this != &other) {
      m_name        = other.getName();
      m_description = other.getDescription();
      m_trueSignal  = other.getTrueSignal();
      if (m_observable) delete m_observable;
      m_observable  = other.clone(); // make a clone - note PDF object is NOT cloned... pointer retained (speed/mem issues)
      other.copyNuisanceList(m_nuisancePars);// idem
    }
  }
  OBS::Base *makeNuisance(OBS::Base *np, PDFN::DISTYPE dist) {
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
  virtual ~MeasPois() {}
  //
  void setObservable(OBS::ObservablePois * obs) { if (m_observable) delete m_observable; m_observable = (obs ? obs->clone():0); }

  virtual const int    getM(double s) { return 0;}
  virtual const double getSignal()    { return 0;}
  virtual const double getSignalUnc() { return 0;}

};

class MeasPoisEB : public MeasPois {
 public:
  MeasPoisEB() : MeasPois() { m_eff=0; m_bkg=0; }
  MeasPoisEB(const char *name, const char *desc=0) : MeasPois(name,desc) { m_eff=0; m_bkg=0; }
  virtual ~MeasPoisEB() { }
  //
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

  void setEffPdf(double eff, double sigma, PDFN::DISTYPE dist) {
    m_eff = static_cast< OBS::BaseType<double> *>(makeNuisance(m_eff,dist));
    m_eff->setPdfMean(eff);
    m_eff->setPdfSigma(sigma);
  }

  void setBkgPdf(double bkg, double sigma, PDFN::DISTYPE dist) {
    m_bkg = static_cast< OBS::BaseType<double> *>(makeNuisance(m_bkg,dist));
    m_bkg->setPdfMean(bkg);
    m_bkg->setPdfSigma(sigma);
  }

  const double getEffObs() { return (m_eff ? m_eff->getObservedValue():0); }
  const double getBkgObs() { return (m_bkg ? m_bkg->getObservedValue():0); }

  OBS::BaseType<double> *getEffPdf() {
    return m_eff;
  }
  const double getEffPdfMean() {
    return (m_eff ? m_eff->getPdfMean():0);
  }
  const double getEffPdfSigma() {
    return (m_eff ? m_eff->getPdfSigma():0);
  }
  const PDFN::DISTYPE getEffPdfDist() {
    return (m_eff ? m_eff->getPdfDist():PDFN::DIST_UNDEF);
  }

  OBS::BaseType<double> *getBkgPdf() {
    return m_bkg;
  }
  const double getBkgPdfMean() {
    return (m_bkg ? m_bkg->getPdfMean():0);
  }
  const double getBkgPdfSigma() {
    return (m_bkg ? m_bkg->getPdfSigma():0);
  }
  const PDFN::DISTYPE getBkgPdfDist() {
    return (m_bkg ? m_bkg->getPdfDist():PDFN::DIST_UNDEF);
  }

  virtual const int    getM(double s) {
    return static_cast<int>(m_eff->getObservedValue()*s + m_bkg->getObservedValue());
  }

  virtual const double getSignal() {
    double dn = m_observable->getObservedValue() - m_bkg->getObservedValue();
    double e = m_eff->getObservedValue(); // >0
    return (e>0 ? dn/e : 0);
  }

 private:
  OBS::BaseType<double> *m_eff; // pointers to nuisance params in list
  OBS::BaseType<double> *m_bkg;
};

#endif
