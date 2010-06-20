#include "Tabulated.h"

TabFun::TabFun() {
  m_data = 0;
}

TabFun::~TabFun() {
  if (m_data) delete [] m_data;
}

void TabFun::init(double xmin, double xmax, int n) {
  if (n>0) {
    std::cout << "Init table with n = " << n << std::endl;
    m_data = new double[n];
    m_xMin = xmin;
    m_xMax = xmax;
    m_nData = n;
    m_dx = (xmax-xmin)/double(n);
    double x;
    for (int i=0; i<n; i++) {
      x = xmin + double(i)*m_dx;
      m_data[i] = F(x);
    }
  }
}

double TabFun::getVal(double x) {
  double rval;
  if ((x<m_xMin)||(x>m_xMax)) {
    rval = F(x);
  } else {
    int sind  = static_cast<int>(x/m_dx);
    if (sind>=m_nData) {
      rval = F(x);
    } else {
      rval = m_data[sind];
    }
  }
  return rval;
}

TabTrig::TabTrig():TabFun() {
}
TabTrig::~TabTrig() {}

double TabTrig::F(double x) {
  return sin(x);
}

double TabTrig::Sin(double x) {
  double rval;
  double xm = x; //mod(x,(2.0*M_PIl));
  int ind;
  //
  if (m_data) {
    ind = int((xm-m_xMin)/m_dx);
    if ((ind<0) && (ind>m_nData)) {
      rval = -m_data[-ind];
    } else if ((ind>=0) && (ind<m_nData)) {
      rval = m_data[ind];
    } else {
      std::cout << "SinOR: " << ind << std::endl;
      rval = F(x);
    }
  } else {
    std::cout << "SinND: " << x << std::endl;
    rval = F(x);
  }
  return rval;
}

double TabTrig::Cos(double x) {
  return Sin(x+M_PI*0.5);
}


TabLog::TabLog():TabFun() {
}
TabLog::~TabLog() {}

double TabLog::F(double x) {
  return log(x);
}

double TabLog::Log(double x) {
  double rval=0;
  int ind;
  //
  if (m_data) {
    ind = int((x-m_xMin)/m_dx);
    if ((ind>=0) && (ind<m_nData)) {
      rval = m_data[ind];
    } else {
      std::cout << "LogOR: " << ind << std::endl;
      if (x<=0) {
	rval = 0;
      } else {
	rval = F(x);
      }
    }
  } else {
      std::cout << "LogND: " << x << std::endl;
    rval = F(x);
  }
  return rval;
}
