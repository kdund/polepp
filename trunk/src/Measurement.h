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
/*! namespace MEAS
  @brief contains the Measurement class

  A measurement consists of the following:
  -# an observable described by OBS::Observable
  -# a list of nuisance parameters (list of OBS::Observable)
  -# the true value of the signal (optional of course - used in e.g coverage studies)
  -# functions that relates the signal quantity (s) with the observable and nuisance parameters
  Measurement<T>::getM()      - given s + nuisance this gives the expected observable using observed values of nuisances
  Measurement<T>::getPdfM()   - idem but the assumed pdf mean values of nuisances are used
  Measurement<T>::getSignal() - given observed nuisances + observable this yields the expected signal
  -# the probability function P( observable | signal and observed nuisances )
     
*/

namespace MEAS {
  template <typename T>
  class Measurement {
  public:
    inline Measurement();
    inline Measurement(const char *name, const char *desc=0);
    inline Measurement(const Measurement<T> & m);

    inline virtual ~Measurement();
    //
    inline Measurement const & operator=( Measurement<T> const & m );
    inline bool operator==( const Measurement<T> & m );

    inline void setTrueSignal(double s);
    inline void setObservable(const OBS::BaseType<T> * obs);
    inline void setObsVal(T val);
    inline void setName(const char *name);
    inline void setDescription(const char *descr);
    //
    inline OBS::Base *addNuisance(OBS::Base * nuPar);
    inline void deleteNuisance();
    inline bool removeNuisance(const OBS::Base *nptr);
    inline void copyNuisance(std::list< OBS::Base * > & newList ) const;

    // REMOVE THESE!!!
//     //! initialize integrals of all nuisance parameters
//     inline void initIntNuisance();
//     //! Assumes that all parameters have called OBS::initInt()
//     inline void fillIntNuisance();
//     //! Initialize nuisance weights: f(x)g(y)...dxdy
//     inline void initNuisanceWeights();
    //
    inline virtual void dump()                     const;
    //
    inline const double getTrueSignal()            const;
    //  const T getObsVal() const;
    inline const T getObsVal()                     const;
    inline const double getObsPdfMean()            const;
    inline const double getObsPdfSigma()           const;
    inline const PDF::DISTYPE getObsPdfDist()      const;
    inline const OBS::BaseType<T> *getObservable() const;
    //
    inline const std::string & getName()                      const;
    inline const std::string & getDescription()               const;
    inline const std::list< OBS::Base * > & getNuisanceList() const;
    inline const double getNuisanceIntNorm()                  const;

    inline const double rndObs();
    //
    inline virtual const double getM(double s)          const;
    inline virtual const double getPdfM(double s)       const;
    inline virtual const double getSignal()             const;
    inline virtual const double getSignalUnc()          const;
    //    inline virtual const double calcProb(T x, double s) const;

    // generate a random pseudoexperiment and store this in the observed parts of the observables
    inline void generatePseudoExperiment();
    //
  protected:
    inline virtual void initObservable() = 0;

    inline void copy(const Measurement<T> & other);

    inline OBS::Base *makeNuisance(PDF::DISTYPE dist);
  
    inline const int getNuisanceIndex(const OBS::Base *nptr);

    ////////////////////////////////////////////////////////

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
    inline MeasPois();
    inline MeasPois(const char *name, const char *desc=0);
    inline MeasPois(OBS::ObservablePois * obs);
    inline MeasPois(const MeasPois & other);
    inline virtual ~MeasPois();
    //
    inline void copy(const MeasPois & other);
    inline void setObservable(const OBS::ObservablePois * obs);
    //
    inline virtual const double getM(double s)            const;
    inline virtual const double getPdfM(double s)         const;
    inline virtual const double getSignal()               const;
    inline virtual const double getSignalUnc()            const;
    //    inline virtual const double calcProb(int x, double s) const;
    //
  protected:
    inline virtual void initObservable();
  };

  class MeasPoisEB : public MeasPois {
  public:
    inline MeasPoisEB();
    inline MeasPoisEB(const char *name, const char *desc=0);
    inline MeasPoisEB(const MeasPoisEB & other);
    inline virtual ~MeasPoisEB();
    //
    inline void copy(const MeasPoisEB & other);

    inline void updNuisanceIndex();

    inline void setEffScale( double scale );
    inline void setBkgScale( double scale );
    inline void setEffInt(double scale, int n);
    inline void setBkgInt(double scale, int n);
    inline void setEffObs(double eff);
    inline void setEffObs();
    inline void setBkgObs(double bkg);
    inline void setBkgObs();
    inline void setEffPdf(double eff, double sigma, PDF::DISTYPE dist);
    inline void setBkgPdf(double bkg, double sigma, PDF::DISTYPE dist);
    inline void setEffPdfMean(double m);
    inline void setEffPdfSigma(double m);
    inline void setBkgPdfMean(double m);
    inline void setBkgPdfSigma(double m);

    inline void setBEcorr(double c); // TODO: Need to implement
    //
    inline const double getBEcorr() const;

    inline const double getEffScale()         const;
    inline const double getEffObs()           const;
    inline const OBS::Base *getEff()          const;
    inline const double getEffPdfMean()       const;
    inline const double getEffPdfSigma()      const;
    inline const PDF::DISTYPE getEffPdfDist() const;

    inline const double getBkgScale()         const;
    inline const double getBkgObs()           const;
    inline const OBS::Base *getBkg()          const;
    inline const double getBkgPdfMean()       const;
    inline const double getBkgPdfSigma()      const;
    inline const PDF::DISTYPE getBkgPdfDist() const;

    inline virtual const double getM(double s)                         const;
    inline virtual const double getM(double s, double eff, double bkg) const;
    inline virtual const double getPdfM(double s)                      const;
    inline virtual const double getSignal()                            const;
    //    inline virtual const double calcProb(int x, double s)              const;

  private:
    OBS::Base *m_eff; //! pointer to efficiency in nuisance params list
    OBS::Base *m_bkg; //! pointer to background in nuisance params list
    int m_effIndex;   //! efficiency index
    int m_bkgIndex;

    double m_effScale;
    double m_bkgScale;
  };

};

#include "Measurement.icc"

#endif
