#include "Range.h"

/////////////////////////////////////////////////////////////////////

Range::Range(double low, double high, double step, int nmax) {
  setRange(low,high,step,nmax);
}

Range::Range() {
  setRange(0,0,0,0);
}

void Range::setRange(double low, double high, double step, int nmax) {
  int n;
  double d = high - low;
  if (d>0) {
    n = static_cast<int>((d + step)/step + 0.5);
  } else {
    high = low;
    step = 1.0;
    n = 1;
  }
  m_min  = low;
  m_max  = high;
  m_n    = n;
  m_step = step;
  if ((nmax>0)&&(n>nmax)) forceNpts(nmax);
}

void Range::forceNpts(int n) {
  if ((n>0) && (n!=m_n)) {
    double d = m_max - m_min;
    m_step = d/static_cast<double>(n);
    m_n = n;
  }
}
