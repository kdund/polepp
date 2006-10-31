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
#include "Tools.h"
#include "Combination.h"
#include "Observable.h"

//
// Note virtual functions:
//  getM(s) : m = f(s,nuisance) ; i.e the relation between measured value, signal and nuisance params - always double
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
  void setObservable(const OBS::BaseType<T> * obs) { if (m_observable) delete m_observable; m_observable = (obs ? obs->clone():0);}
  void setObsVal(T val)                      { if (m_observable) m_observable->setObservedValue(val); }
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

  bool removeNuisance(const OBS::Base *nptr) {
    bool rval=false;
    for (std::list< OBS::Base * >::iterator it = m_nuisancePars.begin();
	 it != m_nuisancePars.end();
	 ++it) {
      if (*it==nptr) {
	it = m_nuisancePars.erase(it);
	rval = true;
      }
    }
    return rval;
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

  //! initialize integrals of all nuisance parameters
  void initIntNuisance() {
    for (std::list< OBS::Base * >::iterator it = m_nuisancePars.begin();
	 it !=  m_nuisancePars.end();
       ++it) {
      (*it)->initIntegral();
    }
  }
  //! Assumes that all parameters have called OBS::initInt()
  void fillIntNuisance() {
    for (std::list< OBS::Base * >::iterator it = m_nuisancePars.begin();
	 it !=  m_nuisancePars.end();
       ++it) {
      (*it)->fillInt();
    }
  }
  //! Initialize nuisance weights: f(x)g(y)...dxdy
  void initNuisanceWeights() {
    std::vector<int> jj;
    std::vector<double> dx;
    double step;
    int n, nnpar;
    m_nuisanceIndecis.clear();
    m_nuisanceIndMax.clear();
    //
    nnpar=0;
    for (std::list< OBS::Base * >::iterator it = m_nuisancePars.begin();
	 it !=  m_nuisancePars.end();
       ++it) {
      nnpar++;
      if ((*it)->isIntFilled()) {
	n = (*it)->getIntN();
	if (n<1) {
	  std::cout << "WARNING: Integral with ZERO points for nuisance parameter -> " << (*it)->getName() << std::endl;
	  n=1; //HERE: SHOULD fill integral for fixed value...
	}
      } else {
	// Could in principle call fillInt() here - but let's not...
	std::cout << "WARNING: Integral not filled for nuisance parameter -> " << (*it)->getName() << std::endl;
	n=1; //HERE: SHOULD fill integral for fixed value...
      }
      m_nuisanceIndMax.push_back(n-1);
      step = (n==1 ? 1.0:(*it)->getIntdX()); // get dx (==1 if only one point)
      dx.push_back((*it)->getIntdX());
      jj.push_back(0);
    }
    //    std::cout << "--> initNuisanceWeights(); nnpar = " << nnpar << std::endl;
    //
    // Loop over all volume elements (dxdydz...)
    // and calculate the weight per volume element: f(x)g(y)h(z)...dxdydz
    //
    //    std::cout << "--> initNuisanceWeights(): loop over all volume elements" << std::endl;
    m_nuisanceWeights.clear();
    double w;
    m_nuisanceIntNorm = 0;
    //
    m_nuisanceIndecis.push_back(jj); // save null vector
    while (Combination::next_vector(jj,m_nuisanceIndMax)) {
      m_nuisanceIndecis.push_back(jj); // save vector
    }
    //
    std::list< OBS::Base * >::iterator itnp = m_nuisancePars.begin();
    //
    //    std::cout << "--> initNuisanceWeights(); nind = " << m_nuisanceIndecis.size() << std::endl;
    double wt;
    for ( unsigned int i=0; i<m_nuisanceIndecis.size(); i++) {
      w = 1.0;
      itnp = m_nuisancePars.begin();
      //! Calculate w(i,j,...) = f(xi)*dx * g(yj)*dy * ...
      for (int j=0; j<nnpar; j++) {
	wt = (*itnp)->getIntWeight(m_nuisanceIndecis[i][j]);
	w *= wt;
        // TODO: DEBUG
        if (i<0) std::cout << "NIND0: par = " << j
                            << "  weight = " << wt
                            << ", index = " << m_nuisanceIndecis[i][j]
                            << ", x = " << (*itnp)->getIntX(m_nuisanceIndecis[i][j])
                            << std::endl;
	++itnp;
      }
      // TODO: Debug info
      if (i<0) std::cout << "NIND: " << std::setw(6) << i << ". weight = " << w << std::endl;
      m_nuisanceWeights.push_back(w);
      m_nuisanceIntNorm += w;
    }
    //    std::cout << "--> initNuisanceWeights(); norm = " << m_nuisanceIntNorm << std::endl;
  }
  //
  void dump() const {
    std::cout << "-----------MEASUREMENT------------------------\n";
    if (m_observable) m_observable->dump();
    std::cout << "----------------------------------------------\n";
  }
  //
  const double getTrueSignal() const        { return m_trueSignal; }
  //  const T getObsVal() const                 { T val=0; if (m_observable) m_observable->getObservedValue(val); return val;}
  const T getObsVal() const                 { T val=0; if (m_observable) m_observable->getObservedValue(val); return val;}
  const double getObsPdfMean() const        { return (m_observable ? m_observable->getPdfMean():0); }
  const double getObsPdfSigma() const       { return (m_observable ? m_observable->getPdfSigma():0); }
  const PDF::DISTYPE getObsPdfDist() const  { return (m_observable ? m_observable->getPdfDist():0); }
  const OBS::BaseType<T> *getObservable() const   { return m_observable; }
  //
  const std::string & getName()                      const { return m_name;}
  const std::string & getDescription()               const { return m_description;}
  const std::list< OBS::Base * > & getNuisanceList() const { return m_nuisancePars; }
  const double getNuisanceIntNorm()                  const { return m_nuisanceIntNorm; }

  const double rndObs() { OBS::BaseType<T> *p = static_cast< OBS::BaseType<T> * >(m_observable); return (*p)(); }
  //
  virtual const double getM(double s) const =0;
  virtual const double getSignal() const =0;
  virtual const double getSignalUnc() const =0;
  virtual const double calcProb(T x, double s) const =0;

 // generate a random pseudoexperiment and store this in the observed parts of the observables
  void generatePseudoExperiment() {
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
    double mean = this->getM(m_trueSignal);
    m_observable->setPdfMean(mean);
    //
    // ...and then the observable
    //
    OBS::BaseType<T> *obs = static_cast< OBS::BaseType<T> * >(m_observable);
    obs->setObservedRnd();
    //
  }
  //
 protected:
  virtual void initObservable()=0;

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

  OBS::Base *makeNuisance(PDF::DISTYPE dist) {
    OBS::Base *np = OBS::makeObservable(dist);
    if (np==0) {
      std::cerr << "ERROR: Failed creating a nuisance parameter!" << std::endl;
    } else {
      addNuisance(np);
    }
    return np;
  }
  
  const int getNuisanceIndex(const OBS::Base *nptr) {
    if (nptr==0) return -1;
    int rval=-1;
    bool found = false;
    std::list< OBS::Base * >::iterator it = m_nuisancePars.begin();
    while (!found && (it != m_nuisancePars.end())) {
      found = (*it==nptr);
      rval++;
      ++it;
    }
    return rval;
  }


  std::string m_name;
  std::string m_description;
  //

  double                                  m_trueSignal;
  OBS::BaseType<T>                       *m_observable;
  std::list< OBS::Base * >                m_nuisancePars;
  std::vector< std::vector<double> >      m_corrMat;
  std::vector< double >                   m_nuisanceWeights; //! nuisance pars weights: f(x)g(y)...dxdy...
  std::vector< std::vector<int> >         m_nuisanceIndecis; //! indecis of (x,y...) -> (i,j,....) used in vector above
  std::vector< int >                      m_nuisanceIndMax;  //! maximum index per parameter
  double				  m_nuisanceIntNorm;
};

class MeasPois : public Measurement<int> {
 public:
  MeasPois() : Measurement<int>() { initObservable(); };
  MeasPois(const char *name, const char *desc=0) : Measurement<int>(name,desc) { initObservable(); };
  MeasPois(OBS::ObservablePois * obs) : Measurement<int>(obs->getName().c_str(),obs->getDescription().c_str()) {
    m_observable = obs->clone();
  }
  MeasPois(const MeasPois & other):Measurement<int>() { copy(other);}
  virtual ~MeasPois() {}
  //
  void copy(const MeasPois & other) { Measurement<int>::copy(other);}
  void setObservable(const OBS::ObservablePois * obs) { if (m_observable) delete m_observable; m_observable = (obs ? obs->clone():0);}

//   void setNObserved(int n) {
//     if (m_observable==0) {
//       m_observable = static_cast< OBS::BaseType<int> *>(OBS::makeObservable(PDF::DIST_POIS));
//       m_observable->setPdfMean(n);
//     }
//     m_observable->setObservedValue(n);
//   }

//  const int    getNObserved() const { return getObsVal(); }

  const double getM(double s) const { return 0;}
  const double getSignal()    const { return 0;}
  const double getSignalUnc() const { return 0;}
  const double calcProb(int x, double s) const {return 0;}
 protected:
  virtual void initObservable() {
    if (m_observable!=0) {
      std::cerr << "MeasPois::FATAL - m_observable non-zero ptr in initObservable()!!! BUG!" << std::endl;
      exit(-1);
    }
    m_observable = static_cast< OBS::BaseType<int> *>(OBS::makeObservable(PDF::DIST_POIS));
  }
};

class MeasPoisEB : public MeasPois {
 public:
  MeasPoisEB() : MeasPois() { m_eff=0; m_bkg=0; m_effScale=1.0; m_bkgScale=1.0; updNuisanceIndex();}
  MeasPoisEB(const char *name, const char *desc=0) : MeasPois(name,desc) { m_eff=0; m_bkg=0; m_effScale=1.0; m_bkgScale=1.0; updNuisanceIndex();}
  MeasPoisEB(const MeasPoisEB & other):MeasPois() {copy(other);}
  virtual ~MeasPoisEB() { }
  //
  void copy(const MeasPoisEB & other) {
    if (this != &other) {
      MeasPois::copy(other);
      //      const OBS::BaseType<double> *eff, *bkg;
      const OBS::Base *eff, *bkg;
      m_eff=0;
      m_bkg=0;
      eff = other.getEff();
      if (eff) m_eff = eff->clone();
      bkg = other.getBkg();
      if (bkg) m_bkg = bkg->clone();
      m_effScale = other.getEffScale();
      m_bkgScale = other.getBkgScale();
    }
  }

  void updNuisanceIndex() { // get's the indices of the bkg and eff in the nuisance list
    m_bkgIndex = getNuisanceIndex(m_bkg);
    m_effIndex = getNuisanceIndex(m_eff);
  }

  void setEffScale( double scale ) {
    if (scale>0) m_effScale = scale;
  }

  void setBkgScale( double scale ) {
    if (scale>0) m_bkgScale = scale;
  }

  void setEffInt(double scale, int n) {
    if (m_eff) {
      m_eff->setIntNpts(n);
      m_eff->setIntScale(scale);
    }
  }

  void setBkgInt(double scale, int n) {
    if (m_bkg) {
      m_bkg->setIntNpts(n);
      m_bkg->setIntScale(scale);
    }
  }

  void setEffObs(double eff) {
    if (m_eff) {
      m_eff->setObservedValue(eff);
    }
  }

  void setEffObs() {
    if (m_eff) {
      m_eff->setObservedValue();
    }
  }

  void setBkgObs(double bkg) {
    if (m_bkg) {
      m_bkg->setObservedValue(bkg);
    }
  }

  void setBkgObs() {
    if (m_bkg) {
      m_bkg->setObservedValue();
    }
  }

  void setEffPdf(double eff, double sigma, PDF::DISTYPE dist) {
    if (m_eff!=0) {
      if (m_eff->getPdfDist()!=dist) {
	removeNuisance(m_eff);
	delete m_eff;
	m_eff=0;
      }
    }
    if (m_eff==0) 
      m_eff = static_cast<OBS::Base *>(makeNuisance(dist));
    m_eff->setPdfSigma(sigma);
    m_eff->setPdfMean(eff);
    m_eff->setName("efficiency");
  
    updNuisanceIndex();
  }

  void setBkgPdf(double bkg, double sigma, PDF::DISTYPE dist) {
    if (m_bkg!=0) {
      if (m_bkg->getPdfDist()!=dist) {
	removeNuisance(m_bkg);
	delete m_bkg;
	m_bkg=0;
      }
    }
    if (m_bkg==0) m_bkg = static_cast< OBS::Base *>(makeNuisance(dist));
    m_bkg->setPdfSigma(sigma);
    m_bkg->setPdfMean(bkg);
    m_bkg->setName("background");
    updNuisanceIndex(); // this will set m_bkgIndex (and m_effIndex if needed)
  }

  void setEffPdfMean(double m)  { if (m_eff) m_eff->setPdfMean(m); }
  void setEffPdfSigma(double m) { if (m_eff) m_eff->setPdfSigma(m); }

  void setBkgPdfMean(double m)  { if (m_bkg) m_bkg->setPdfMean(m); }
  void setBkgPdfSigma(double m) { if (m_bkg) m_bkg->setPdfSigma(m); }

  void setBEcorr(double c) {} // TODO: Need to implement
  //
  const double getEffScale() const { return m_effScale; }
  const double getBkgScale() const { return m_bkgScale; }

  const double getBEcorr() const { return 0.0; } // TODO

  const double getEffObs() const { return (m_eff ? m_eff->getObservedValue():0); }

  const double getBkgObs() const { return (m_bkg ? m_bkg->getObservedValue():0); }

//   const OBS::BaseType<double> *getEff() const {
//     return m_eff;
//   }
  const OBS::Base *getEff() const {
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
  const OBS::Base *getBkg() const {
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

  const double getM(double s) const {
    double e,b;
    e = m_eff->getObservedValue()*m_effScale;
    b = m_bkg->getObservedValue()*m_bkgScale;
    return e*s + b;
  }

  const double getM(double s, double eff, double bkg) const {
    return m_effScale*eff*s + m_bkgScale*bkg;
  }

  const double getSignal() const {
    double e,b;
    e = m_eff->getObservedValue()*m_effScale;
    b = m_bkg->getObservedValue()*m_bkgScale;;
    double dn = static_cast<OBS::Base *>(m_observable)->getObservedValue() - b;
    if (dn<0.0) dn=0.0;
    return (e>0 ? dn/e : 0);
  }

  const double calcProb(int x, double s) const {
    if ((m_effIndex<0) ||
	(m_bkgIndex<0) ||
	(m_eff==0) ||
	(m_bkg==0) ||
	(m_observable==0)) {
      std::cerr << "ERROR: MeasPoisEB() - No eff and/or bkg defined!" << std::endl;
      std::cerr << "       EffPtr = " << m_eff << std::endl;
      std::cerr << "       BkgPtr = " << m_bkg << std::endl;
      std::cerr << "       EffInd = " << m_effIndex << std::endl;
      std::cerr << "       BkgInd = " << m_bkgIndex << std::endl;
      std::cerr << "       Observ = " << m_observable << std::endl;

      return 0.0;
    }
    ////// BUG TRAP!
    if (m_observable->getPdf()==0) {
      std::cerr << "ERROR: MeasPoisEB() - Observable OK, but with no PDF!" << std::endl;
      return 0.0;
    }
    int ixeff, ixbkg;
    double p=0;
    double g;
    PDF::BaseType<int> *pdf = static_cast<  PDF::BaseType<int> * >(m_observable->getPdf());
    //////////////// TEMP CODE - FOR DEBUGGING, IN CASE OF BUGS...
    if (pdf==0) {
      std::cerr << "ERROR: wrong type of PDF in MeasPoisEB (expected BaseType<int>)!" << std::endl;
      return 0;
    }
    if (m_nuisanceWeights.size() != m_nuisanceIndecis.size()) {
      std::cerr << "ERROR: nuisance weight and index arrays are not of equal size!" << std::endl;
      return 0.0;
    }
    if (m_nuisanceWeights.size()==0) {
      std::cerr << "ERROR: nuisance weight and index arrays are not initialized (size zero)!" << std::endl;
      return 0.0;
    }
    for (unsigned int i=0; i<m_nuisanceIndecis.size(); i++) {
      ixeff = m_nuisanceIndecis[i][m_effIndex];
      ixbkg = m_nuisanceIndecis[i][m_bkgIndex];
      g = getM(s, m_eff->getIntX(ixeff), m_bkg->getIntX(ixbkg));
      p += m_nuisanceWeights[i]*pdf->getVal(x,g);
      if (x<0) { // TODO: Debug info - remove when obsolete
        TOOLS::coutFixed(4,int(i));
        TOOLS::coutFixed(" s(true), e, b, mu = ",2,s);
        TOOLS::coutFixed(", ", 4, m_eff->getIntX(ixeff));
        TOOLS::coutFixed(", ", 4, m_bkg->getIntX(ixbkg));
        TOOLS::coutFixed(", ", 4, g);
        TOOLS::coutFixed(", p = ", 4, p);
        TOOLS::coutFixed(", w = ", 4, m_nuisanceWeights[i]);
        TOOLS::coutFixed(", pdfval = ", 4, pdf->getVal(x,g));
        std::cout << std::endl;
      }
    }

    //    std::cout << "calcProb( " << x << ", " << s << " ) = " << p/m_nuisanceIntNorm << std::endl;
    return p/m_nuisanceIntNorm;
  }

 private:
//   OBS::BaseType<double> *m_eff; // pointers to nuisance params in list
//   OBS::BaseType<double> *m_bkg;
  OBS::Base *m_eff; // pointers to nuisance params in list
  OBS::Base *m_bkg;
  int m_effIndex;
  int m_bkgIndex;

  double m_effScale;
  double m_bkgScale;
};

#endif
