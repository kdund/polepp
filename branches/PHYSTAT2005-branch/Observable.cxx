#include "Observable.h"

void PdfGauss::tabulate() {
  if (m_tableNpts<=0) return;
  if (m_table)        return; // table already defined
  //
  double mu;
  m_tabDmu = m_tabMuMax/double(m_tableNpts);
  //
  m_table = new double[m_tableNpts];
  for (int i=0; i<m_tableNpts; i++) {
    mu = m_tabDmu*double(i);
    m_table[i] = phi(mu);
  }
}

void PdfPoisson::tabulate() {
  if (m_tableNpts<=0) return;
  if (m_table)        return; // table already defined
  //
  m_table = new double[m_tableNpts];
  //
  double lmb;
  unsigned long index;
  for (int i=0; i<m_tabLmbN; i++) {
    lmb = m_tabDlmb*double(i);
    for (int j=0; j<m_tabNmax; j++) {
      index = i*m_tabNmax +j;
      m_table[index] = rawPoisson(int(j),lmb);
    }
  }
}

double PdfPoisson::rawPoisson(int n, double s) {
  double prob;
  if(s<50.0) {
    prob = (pow(s,n)/exp(lgamma(n+1)))*exp(-s);
  } else {
    double sigma = sqrt(s); // gaussian aprox.
    double c = 1.0L/(sqrt(2.0*M_PI)*sigma);
    double t = (double(n)-s)/sigma;
    prob = c*exp(-0.5L*t*t);
  }
  return prob;
}
