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
Poisson::Poisson() {
  m_data = 0; // for the destructor
  m_nN = 0;   // makes sure if we call getVal without having done init(), it will still function
  m_nTot = 0;
  m_nLambda = 0;
}

Poisson::~Poisson() {
  if (m_data) delete [] m_data;
}

void Poisson::init(int nlambda, int nn, double lmax) {
  m_nLambda = nlambda;
  m_nN = nn;
  m_nTot = nn*nlambda;
  m_lambdaMax = lmax;
  //
  m_data = new double[m_nTot];
  //
  m_dL = m_lambdaMax/double(m_nLambda);
  //
  double lmb;
  unsigned long index;
  std::cout << "==================================" << std::endl;
  if (m_nTot==0) {
    std::cout << " Will not use tabulated poisson!" << std::endl;
  } else {
    std::cout << " Initialising poisson table with:" << std::endl;
    std::cout << "   n(lambda)  = " << m_nLambda <<std::endl;
    std::cout << "   lambda_max = " << m_lambdaMax <<std::endl;
    std::cout << "   n(N)       = " << m_nN <<std::endl;
    for (unsigned long i=0; i<m_nLambda; i++) {
      lmb = m_dL*double(i);
      for (unsigned long j=0; j<m_nN; j++) {
	index = i*m_nN +j;
	m_data[index] = rawPoisson(int(j),lmb);
      }
    }
  }
  std::cout << "==================================\n" << std::endl;
}

double Poisson::rawPoisson(int n, double s) {
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

Gauss::Gauss() {
  m_data = 0;    // for the destructor
  m_nData = 0;   // makes sure if we call getVal without having done init(), it will still function
}

Gauss::~Gauss() {
  if (m_data) delete [] m_data;
}

void Gauss::init(int ndata, double mumax) {
  double mu;

  m_nData = ndata;
  m_muMax = mumax;
  m_dMu = m_muMax/double(m_nData);
  //
  std::cout << "==================================" << std::endl;
  if (m_nData==0) {
    std::cout << " Will not use tabulated gauss!" << std::endl;
  } else {
    m_data = new double[m_nData];
    std::cout << " Initialising gauss table with:" << std::endl;
    std::cout << "   n data     = " << m_nData <<std::endl;
    std::cout << "   mu max     = " << m_muMax <<std::endl;
    for (unsigned long i=0; i<m_nData; i++) {
      mu = m_dMu*double(i);
      m_data[i] = rawGauss(mu);
    }
  }
  std::cout << "==================================\n" << std::endl;
}

double Gauss::rawGauss(double mu) {
  double prob;
  double norm = 1.0L/sqrt(2.0*M_PI);
  prob = norm*exp(-0.5L*mu*mu);
  return prob;
}

// LogNormal::LogNormal() {
//   m_nmean=1.0;
//   m_nsigma=1.0;
// }

// LogNormal::~LogNormal() {
// }

// double LogNormal::init(double mean, double sigma) {
//   m_nmean  = log(mean*mean/sqrt(sigma*sigma + mean*mean));
//   m_nsigma = sqrt(log(sigma*sigma/mean/mean+1));
// }
