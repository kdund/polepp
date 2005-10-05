#ifndef TABULATED_H
#define TABULATED_H

#include <iostream>
#include <cmath>

template <typename T>
class TabFun {
public:
  TabFun();
  virtual ~TabFun();
  //
  void init(double xmin, double xmax, int n);
  virtual double F(double x);
  virtual double getVal(double x);
  //
protected:
  int    m_nData;
  double m_xMin;
  double m_xMax;
  double m_dx;
  double *m_data;
};

class TabTrig:public TabFun {
public:
  TabTrig();
  virtual ~TabTrig();
  //
  virtual double F(double x);
  virtual double Cos(double x);
  virtual double Sin(double x);
};

class TabLog:public TabFun {
public:
  TabLog();
  virtual ~TabLog();
  //
  virtual double F(double x);
  virtual double Log(double x);
};

#endif

