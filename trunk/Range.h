#ifndef RANGE_H
#define RANGE_H

#include <iostream>
#include <vector>
#include <cmath>

class Range {
public:
  Range(double vmin, double vmax, double vstep, double stepMin=0.001);
  Range();
  ~Range() {}
  //
  void setRange(double vmin, double vmax, double vstep, double stepMin=0.001);
  void forceNpts(int n);
  inline double getVal(int index, double def=0);
  inline double min()  {return m_min;}
  inline double max()  {return m_max;}
  inline int    n()    {return m_n;}
  inline double step() {return m_step;}
private:
  double m_min;
  double m_max;
  double m_step;
  int    m_n;
};

inline double Range::getVal(int index, double def) {
  double rval=def;
  if ((index<m_n) && (index>=0))
    rval = m_min + static_cast<double>(index)*m_step;
  return rval;
}

#endif
