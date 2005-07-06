#ifndef RANGE_H
#define RANGE_H

#include <iostream>
#include <vector>
#include <cmath>

class Range {
public:
  Range(double vmin, double vmax, double vstep, int nmax=0);
  Range();
  ~Range() {}
  //
  void setRange(double vmin, double vmax, double vstep, int nmax=0);
  void forceNpts(int n);
  inline const double getVal(int index, double def=0) const;
  inline const double min() const  {return m_min;}
  inline const double max() const  {return m_max;}
  inline const int    n() const    {return m_n;}
  inline const double step() const {return m_step;}
private:
  double m_min;
  double m_max;
  double m_step;
  int    m_n;
};

inline const double Range::getVal(int index, double def) const {
  double rval=def;
  if ((index<m_n) && (index>=0))
    rval = m_min + static_cast<double>(index)*m_step;
  return rval;
}

#endif
