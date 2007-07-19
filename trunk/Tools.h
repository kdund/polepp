#ifndef POLE_TOOLS_H
#define POLE_TOOLS_H
#include <iostream>
#include <iomanip>
#include <ctime>
#include <cmath>
#include <string>
#include <sstream>

namespace OBS {
  class Base;
};

namespace TOOLS {
  inline const char *yesNo(bool var);
  inline const char *boolChar(bool var, const char *tch, const char *fch);
  inline void coutFixed(int precision, int val);
  inline void coutFixed(int precision, double val);
  inline void coutFixed(const char *msg, int precision, int val);
  inline void coutFixed(const char *msg, int precision, double val);
  void makeTimeStamp( std::string & stamp );
  void makeTimeStamp( std::string & stamp, time_t time );
  void makeTimeStamp( std::string & stamp, struct tm *time );
  inline void calcFlatRange( double mean, double sigma, double & xmin, double & xmax );
  inline void calcFlatMeanSigma( double xmin, double xmax,  double &mean, double & sigma );
  void calcIntRange(const OBS::Base & obs, double scale,  double & xmin, double & xmax );

  class Timer {
  public:
    Timer():m_runningTime(false),m_runningClock(false)  {}
    ~Timer() {}

    void start() { startTimer(); startClock(); }
    void stop() { stopTimer(); stopClock(); }
    void clear() { m_startTime=0; m_stopTime=0; m_estTime=0; m_startClock=0; m_stopClock=0; m_runningTime = false; m_runningClock = false; }
    void startTimer();
    bool checkTimer(int dt);
    void stopTimer();
    time_t calcEstimatedTime( int nloops, int ntotal );
    void printTime(const char *msg, time_t t );
    void printCurrentTime(const char *msg);
    void printUsedTime();
    void printEstimatedTime(int nloops, int ntotal );
    void startClock();
    void stopClock();
    void printUsedClock(int norm);

    clock_t getStartClock() { return m_startClock; }
    clock_t getStopClock()  { return (m_runningClock ? clock():m_stopClock); }
    double  getUsedClock( double scale=1e-3 ) { return (getStopClock() - m_startClock)*scale; }

    time_t getStartTime() { return m_startTime; }
    time_t getStopTime()  { time_t t; time(&t); return (m_runningTime ? t:m_stopTime); }
    time_t getUsedTime()  { return getStopTime() - m_startTime; }
    time_t getEstTime()   { return m_estTime; }

  private:
    bool   m_runningTime;
    bool   m_runningClock;

    time_t m_startTime;
    time_t m_stopTime;
    time_t m_estTime;

    clock_t m_startClock;
    clock_t m_stopClock;
  };
};

void TOOLS::calcFlatRange( double mean, double sigma, double & xmin, double & xmax ) {
  const double d = sqrt(12.0);
  double wh = 0.5*d*sigma;
  xmin = mean-wh;
  xmax = mean+wh;
}

void TOOLS::calcFlatMeanSigma(double xmin, double xmax, double & mean, double & sigma) {
  mean  = (xmax+xmin)/2.0;
  sigma = (xmax-xmin)/sqrt(12.0);
}


const char * TOOLS::yesNo(bool var) {
  static char yes[4] = {'Y','e','s',char(0)};
  static char no[3] = {'N','o',char(0)};
  return (var ? &yes[0]:&no[0]);
}

const char * TOOLS::boolChar(bool var, const char *tch, const char *fch) {
  return (var ? tch:fch);
}
//
// allow for gcc 2.96 fix
//
# if __GNUC__ > 2
void TOOLS::coutFixed(int precision, double val) {
  std::ios_base::fmtflags old = std::cout.flags();
  std::cout << std::fixed << std::setprecision(precision) << val;
  std::cout.flags(old);
}
void TOOLS::coutFixed(int precision, int val) {
  std::ios_base::fmtflags old = std::cout.flags();
  std::cout << std::fixed << std::setprecision(precision) << val;
  std::cout.flags(old);
}

void TOOLS::coutFixed(const char *msg, int precision, double val) {
  std::cout << msg;
  std::ios_base::fmtflags old = std::cout.flags();
  std::cout << std::fixed << std::setprecision(precision) << val;
  std::cout.flags(old);
}
void TOOLS::coutFixed(const char *msg, int precision, int val) {
  std::cout << msg;
  std::ios_base::fmtflags old = std::cout.flags();
  std::cout << std::fixed << std::setprecision(precision) << val;
  std::cout.flags(old);
}
# else
void TOOLS::coutFixed(int precision, double val) {
  static char fmt[32];
  sprintf(&fmt[0],"%%.%df",precision);
  printf(fmt,val);
}

void TOOLS::coutFixed(int precision, int val) {
  static char fmt[32];
  sprintf(&fmt[0],"%%%dd",precision);
  printf(fmt,val);
}


#endif

#endif
