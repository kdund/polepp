#ifndef RANGE_H
#define RANGE_H

#include <iostream>
#include <vector>
#include <cmath>

template <typename T>
class Range {
public:
  Range(T vmin, T vmax, T vstep, int n=0) { setRange(vmin,vmax,vstep,n); }
  Range();
  Range(const Range<T> & other) { copy(other); }
  ~Range() {}
  //
  void copy(const Range<T> & other) {
    if (this != &other) {
      m_min  = other.min();
      m_max  = other.max();
      m_step = other.step();
      m_n    = other.n();
    }
  }
  void setRange(T vmin, T vmax, T vstep, int n=0) {
    if (high<=low) {
      m_min = low;
      m_max = high;
      m_step = 0;
      m_n = 1;
      return;
    }
    //
    bool defN=(step<=0);
    m_min  = low;
    m_max  = high;
    if (defN && (n<2)) n=2;
    //
    T d = high-low;
    if (defN) { // N is defined
      m_step = d/static_cast<T>(n-1);
      m_n = n;
    } else { // step is defined
      m_n = 1+static_cast<int>(d/m_step);
      m_step = step;
    }
  }
  //
  const T getVal(int index, T def=0) const {
    T rval=def;
    if ((index<m_n) && (index>=0))
      rval = m_min + static_cast<T>(index)*m_step;
    return rval;
  }

  const T min()  const {return m_min;}
  const T max()  const {return m_max;}
  const int n()  const {return m_n;}
  const T step() const {return m_step;}
private:
  T m_min;
  T m_max;
  T m_step;
  int    m_n;
};

#endif
