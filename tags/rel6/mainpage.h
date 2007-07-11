// Mainpage for Doxygen

/** @mainpage package PoleLib
 *
 * @author Fredrik Tegenfeldt, Iowa State University (fredrik.tegenfeldt@cern.ch)
 *
 * @section intro Introduction
 *
 * This library contains a method for calculating limits using the Feldman-Cousins likelihood ratio ordering and
 * a Bayesian treatment of nuisance parameters (Highland-Cousins). Some relevant publications are listed
 * in \ref refs.
 *
 * @subsection Method
 * The method assumes that a measurement is described as follows:
 *  \f[N_{obs} = \epsilon s + b\f]
 *  where
 *  - \f$N_{obs}\f$ = number of observed events (Poisson)
 *  - \f$\epsilon\f$ = measured efficiency with a known (or assumed) distribution
 *  - \f$s\f$ = unknown signal for which the limits are to be found
 *  - \f$b\f$ = measured background with a known (or assumed) distribution
 *
 * In some cases, the assumed distribution of the nuisance parameters might be given as scaled values
 * of those entering the above equation. That is, the following situation:
 *  \f[N_{obs} = c_{\epsilon}\epsilon' s + c_{b}b'\f]
 * where the distributions are defined for \f$\epsilon'\f$ and \f$b'\f$. A possible situation is, e.g,
 * that the efficiency estimate is given by a poisson with mean 10 and that the actual efficiency is 1/10 of that value.
 *
 *  The pdf used for describing the number of observed events is:  
 *  \f[
 *  q(n)_{s+b} =  
 *  \frac{1}{2\pi\sigma_{b}\sigma_{\epsilon}} \times
 *  \intop_0^{\infty}\intop_0^{\infty}p(n)_{b+ 
 *  \epsilon's}\;\;
 *  f(b'|b,\sigma_{b})\;\;
 *  g(\epsilon'|\epsilon,\sigma_{\epsilon})\;\;
 *  db'd\epsilon'
 *  \f]
 *  where \f$f()\f$ and \f$g()\f$ are the pdf for the background and efficiency respectively.
 *
 *  In order to find the lower and upper limits, a so called confidence belt is constructed.
 *  For each value of \f$s\f$, \f$n_1\f$ and \f$n_2\f$ are found such that
 *  \f$\sum_{n=n_1}^{n_2} q(n)_{s+b} = 1 - \alpha\f$ where \f$1-\alpha\f$ is the confidence limit.
 *  The confidence belt is then given by the set \f$[n_1(s+b,\alpha),n_2(s+b,\alpha)]\f$.
 *  An upper and lower limit is found by finding the intersection of the vertical line \f$n = N_{obs}\f$
 *  and the boundaries of the belt (not exactly true due to the discreteness of a Poisson).
 *
 *  However, the confidence belt is not unambigously defined. A specific ordering scheme is required.
 *  The method used for selecting \f$n_1\f$ and
 *  \f$n_2\f$ is based on likelihood ratios (aka Feldman & Cousins). For each n, a \f$s_{best}\f$ is found that
 *  will maximise the likelihood \f$\mathcal{L}(n)_{s+b}\f$ (mathematically it's just \f$q(n)_{s+b}\f$,
 *  but here it's not used as a pdf but rather as a hypothesis). For a fixed s, a likelihood ratio is calculated
 *  \f[R(n,s)_{\mathcal{L}} = \frac{\mathcal{L}(n|\epsilon s + b)}{\mathcal{L}(n|\epsilon s_{best} + b)}\f]
 *  and the n are included in the sum starting with the highest rank (ratio) and continuing with decreasing rank until
 *  the sum (ref) equals the requested confidence.
 *
 *  A problem with the above method is that the upper limit may \it{decrease}
 *  with \it{increased} \f$\sigma_{b}\f$. There is a proposed remedy for this problem by Gary Hill.
 *  It modifies the likelihood ratio in the denominator. Instead of scanning for the best hypothesis maximising the likelihood,
 *  the ratio is modified as follows:
 *
 *  \f[
 *  R(n,s)_{\mathcal{L}} = \frac{\mathcal{L}(n|\epsilon s + b)}{\mathcal{L}(n|\max(0,n_{obs}-\hat{b}) + \hat{b})}
 *  \f]
 *  where \f$\max()\f$ and \f$\hat{b}\f$ are the maximum likelihood estimates of \f$s_{true}\f$
 *  and background, respectively, given the measurement.
 *
 *  
 *
 * @section setup Setup
 * 
 *  \b Basic
 *  - Pole::setMethod()\n
 *    Sets the Likelihood ratio method (MBT=2 or FHC2=1).
 *  - Pole::setCL()\n
 *    Confidence limit, default 0.9, the requested confidence [0.0,1.0].
 *  - Pole::setNObserved()\n
 *    Number of observed events.
 *  - Pole::setEffPdf()\n
 *    Measured efficiency and assumed distribution (mean,sigma and distribution type ( PDF::DISTYPE )).
 *  - Pole::setBkgPdf()\n
 *    Measured background and assumed distribution (mean,sigma and distribution type ( PDF::DISTYPE )).
 *  - Pole::setEffBkgPdfCorr()\n
 *    Correlation between eff and bkg (if applicable), correlation coefficient [-1.0,1.0].
 *  - Pole::setEffObs()\n
 *    Set the observed efficiency. Without argument, it will just set the observed value to the mean given by setEffPdf().
 *  - Pole::setBkgObs()\n
 *    Set the observed background. For behaviour, see previous.
 *  - Pole::setEffPdfScale()\n
 *    Scaling of efficiency.
 *  - Pole::setBkgPdfScale()\n
 *    Scaling of background.
 *
 *  \b Integration
 *  - Pole::setEffInt()\n
 *    Integration range of efficiency.\n
 *    The integration range must cover the pdf such that the tails are negligable.
 *  - Pole::setBkgInt()\n
 *    Ditto, background.
 *
 *  \b Belt construction
 *  - Pole::calcBelt()\n
 *    Calculates the confidence belt [n1(s,b),n2(s,b)].
 *
 *  \b Finding best hypothesis
 *  - Pole::setBestMuStep()\n
 *    Sets the precision in findBestMu().
 *    Default = 0.01; should normally be fine.\n
 *    Only used if method is RLMETHOD::RL_FHC2.
 *
 *  \b Hypothesis testing
 *  - Pole::setTestHyp()\n
 *    This sets the hypothesis range and step size (precision) when calculating the construction or confidence belt.\n
 *    For limit calculations, this setting has no effect.
 *
 *  \b Coverage related
 *  - Pole::setTrueSignal()\n
 *    Set the true signal.
 *  - Pole::setUseCoverage()\n
 *    Setting this flag causes the limit calculation to terminate the scan as soon as it is
 *    decided whether the true signal is inside or outside the confidence limit.
 *
 *  \b Pdf
 *  - Pole::setPoisson()\n
 *    Initialises the poisson generator - use PDF::gPoisTab.
 *  - Pole::setGauss()\n
 *    Ditto but for a gauss pdf - use PDF::gGauss.
 *  - Pole::setGauss2D()\n
 *    Ditto but for a 2d gauss pdf - use PDF::gGauss2D.
 *  - Pole::setLogNormal()\n
 *    Ditto but for a log normal pdf - use PDF::gLogNormal.
 *
 *  See the files argsPole.cxx and argsCoverage.cxx for examples of usage.
 *
 * @section running Running
 *
 *  \b Print
 *  - Pole::printSetup()\n
 *    Prints the setup to stdout.
 *  - Pole::printLimit()\n
 *    Prints the calculated limit.
 *
 *  \b General
 *  - Pole::execute()\n
 *    Main routine to call. Initialises and runs with the current setup.
 *  - Pole::analyseExperiment()\n
 *    Calculates the limit using the current setup.
 *
 *  \b Debug
 *  - Pole::setVerbose()\n
 *    Sets verbose level. Produces a lot of difficult printouts. A useful level might be -V 2.
 *    That will print out the scanning for the lower and upper limit.
 *
 * @section tools Tools
 *  - polelim : limit calculation
 *  - polecov : coverage calculator
 *  - polebelt : confidence belt calculator
 *  - poleconst : prints the full construction
 *
 * @subsection tools_polelim polelim
 *
 * @section refs References
 *<table border="0" cellpadding="5" cellspacing="0">
 * <tr> <td><b> LHR ordering   </b></td><td><em> G.J. Feldman, R.D. Cousins </em></td><td> Phys.Rev D57 3873 (1998)      </td> </tr>
 * <tr> <td><b> Pole (FHC2)    </b></td><td><em> J.Conrad et al.            </em></td><td> Phys.Rev D67 012002 (2003)    </td> </tr>
 * <tr> <td><b> MBT            </b></td><td><em> G.Hill                     </em></td><td> Phys.Rev D67 118101 (2003)    </td> </tr>
 * <tr> <td><b> Coverage study </b></td><td><em> J.Conrad, F.Tegenfeldt     </em></td><td> NIM A 539 1-2 p407-413 (2005) </td> </tr>
 *</table>
 *
 *<hr>
 * @section notes  Notes for current release
 * release.notes
 *<hr>
 */
