#include "Tools.h"

namespace TOOLS {
  void makeTimeStamp( std::string & stamp ) {
    time_t curtime;
    time(&curtime);
    struct tm *loctime;
    loctime = localtime(&curtime);
    makeTimeStamp( stamp, loctime );
  }
  void makeTimeStamp( std::string & stamp, time_t t ) {
    struct tm *timeStruct = localtime(&t);
    makeTimeStamp( stamp, timeStruct );
  }
  void makeTimeStamp( std::string & stamp, struct tm *t ) {
    static char ts[32];
    strftime( ts, 32,"%d/%m/%Y %H:%M:%S",t);
    stamp = ts;
  }
  //
  // class Timer members
  //
  void Timer::startTimer() {
    //
    m_stopTime = 0;
    m_estTime  = 0;
    //
    time(&m_startTime);
    std::string tstamp;
    TOOLS::makeTimeStamp( tstamp, m_startTime );
    std::cout << "Start of run: " << tstamp << std::endl;
  }

  bool Timer::checkTimer(int dt) {
    time_t chk;
    int delta;
    time(&chk);
    delta = int(chk - m_startTime);
    return (delta>dt); // true if full time has passed
  }

  void Timer::stopTimer() {
    time(&m_stopTime);
  }

  void Timer::printTime(const char *msg, time_t t ) {
    std::string tstamp;
    makeTimeStamp( tstamp, t );
    std::cout << msg << tstamp << std::endl;
  }

  void Timer::printCurrentTime(const char *msg) {
    std::string tstamp;
    makeTimeStamp( tstamp );
    std::cout << msg << tstamp << std::endl;
  }

  void Timer::printUsedTime() {
    if (m_startTime<=m_stopTime) {
      time_t loopTime  = m_stopTime - m_startTime;
      int hours,mins,secs;
      hours = static_cast<int>(loopTime/3600);
      mins = (loopTime - hours*3600)/60;
      secs =  loopTime - hours*3600-mins*60;
      std::cout << "Used time: ";
      std::cout << hours << "h " << mins << "m " << secs << "s" << std::endl;
    }
  }

  void Timer::printEstimatedTime(int nloops, int ntotal ) {
    if (m_startTime<=m_stopTime) {
      time_t loopTime  = m_stopTime - m_startTime;
      time_t deltaTime = (ntotal*loopTime)/ntotal;
      //
      m_estTime = deltaTime + m_startTime;
      std::string tstamp;
      makeTimeStamp( tstamp, m_estTime );
      //
      int hours,mins,secs;
      hours = static_cast<int>(deltaTime/3600);
      mins = (deltaTime - hours*3600)/60;
      secs =  deltaTime - hours*3600-mins*60;
      std::cout << "Estimated end of run: " << tstamp;
      std::cout << " ( " << hours << "h " << mins << "m " << secs << "s" << " )\n" << std::endl;
    }
  }
  void Timer::startClock() {
    m_startClock = clock();
  }
  void Timer::stopClock() {
    m_stopClock = clock();
  }
  void Timer::printUsedClock(int norm) {
    clock_t t = m_stopClock - m_startClock;
    std::cout << "Total CPU time used (ms)     : " << t/1000.0 << std::endl;
    if (norm>0) std::cout << "Per event CPU time used (ms) : " << (t/norm)/1000.0 << std::endl;
  }

};
