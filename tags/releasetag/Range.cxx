#include "Range.h"

/////////////////////////////////////////////////////////////////////

Range::Range(double low, double high, double step, double stepMin) {
  setRange(low,high,step,stepMin);
}

Range::Range() {
  setRange(0,0,0,0);
}

void Range::setRange(double low, double high, double step, double stepMin) {
  int n;
  double d = high - low;
  if (step<stepMin) step=stepMin;
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
}

void Range::forceNpts(int n) {
  if ((n>0) && (n!=m_n)) {
    double d = m_max - m_min;
    m_step = d/static_cast<double>(n);
    m_n = n;
  }
}
