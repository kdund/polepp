#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include <TROOT.h>
#include <TApplication.h>

#include <TCanvas.h>
#include <TPad.h>
#include <TH1.h>
#include <TH2.h>

#include <TFile.h>
#include <TTree.h>
#include <TLegend.h>


std::vector<double> gXvec;
std::vector<double> gYvec;
std::vector<double> gZvec;

bool getXYZ(std::string data, double & x, double & y, double & z) {
  std::istringstream sstr(data);
  sstr >> x >> y >> z;
  return sstr;
}

bool loadData(std::string filename) {
  bool rval=false;
  ifstream inpf(filename.c_str());
  std::string dataLine;
  double x,y,z;
  int nch;
  //
  if (inpf.is_open()) {
    while ((nch=inpf.peek())>-1) {
      getline(inpf,dataLine);
      if (getXYZ(dataLine,x,y,z)) {
	gXvec.push_back(x);
	gYvec.push_back(y);
	gZvec.push_back(z);
      }
    }
    rval = (gXvec.size()>0);
  }
  return rval;
}

TH1F *makeH1(const char *name, const char *title, std::vector<double> x, int nbins=0, double xmin=0, double xmax=0) {
  TH1F *rval=0;
  if (x.size()==0) return rval;
  //
  if (((xmax==0)&&(xmin==0))|| (xmax<xmin)) {
    xmin=x[0];
    xmax=x[0];
    for (unsigned i=0; i<x.size(); i++) {
      if (x[i]>xmax) xmax=x[i];
      if (x[i]<xmin) xmin=x[i];
    }
  }
  if (nbins<=0) {
    nbins = static_cast<int>(xmax-xmin+0.5);
  }
  rval = new TH1F(name, title, nbins, xmin, xmax);
  for (unsigned i=0; i<x.size(); i++) {
    rval->Fill(x[i]);
  }
  return rval;
}

TH2F *makeH2(const char *name, const char *title, std::vector<double> x, std::vector<double> y, int nbinsx=0, int nbinsy=0, double xmin=0, double xmax=0, double ymin=0, double ymax=0) {
  TH2F *rval=0;
  if ((x.size()==0)||(y.size()==0)||(y.size()!=x.size())) return rval;
  //
  if (((xmax==0)&&(xmin==0))|| (xmax<xmin)) {
    xmin=x[0];
    xmax=x[0];
    for (unsigned i=0; i<x.size(); i++) {
      if (x[i]>xmax) xmax=x[i];
      if (x[i]<xmin) xmin=x[i];
    }
  }
  if (((ymax==0)&&(ymin==0))|| (ymax<ymin)) {
    ymin=y[0];
    ymax=y[0];
    for (unsigned i=0; i<y.size(); i++) {
      if (y[i]>ymax) ymax=y[i];
      if (y[i]<ymin) ymin=y[i];
    }
  }
  if (nbinsx<=0) {
    nbinsx = static_cast<int>(xmax-xmin+0.5);
  }
  if (nbinsy<=0) {
    nbinsy = static_cast<int>(ymax-ymin+0.5);
  }
  rval = new TH2F(name, title, nbinsx, xmin, xmax, nbinsy, ymin, ymax);
  for (unsigned i=0; i<x.size(); i++) {
    rval->Fill(x[i],y[i]);
  }
  return rval;
}

void plotData() {
  TCanvas *canvas = new TCanvas("H1","H1",200,10,600,280);
  canvas->cd();
  //
  //  (0,1)------------(1,1)
  //    |                |
  //    |     Canvas     |
  //    |                |
  //  (0,0) -----------(1,0)
  //
  TPad *pad1 = new TPad("pad1","The pad with the function",0.03,0.62,0.33,0.92,21);
  TPad *pad2 = new TPad("pad2","The pad with the histogram",0.34,0.62,0.64,0.92,21);
  TPad *pad3 = new TPad("pad3","The pad with the histogram",0.65,0.62,0.95,0.92,21);
  TPad *pad4 = new TPad("pad4","The pad with the histogram",0.03,0.02,0.98,0.57);
  pad1->Draw();
  pad2->Draw();
  pad3->Draw();
  pad4->Draw();

  pad1->cd();
  TH1F *histEff = makeH1("Eff", "Efficiency",gYvec,100);
  histEff->Draw();

  pad2->cd();
  TH1F *histBkg = makeH1("Bkg", "Background",gZvec,100);
  histBkg->Draw();

  pad3->cd();
  TH2F *histEffBkg = makeH2("EffBkg", "Correlation",gYvec,gZvec,100,100);
  histEffBkg->Draw();

  pad4->cd();
  TH1F *histNobs = makeH1("Nobs", "Number of observed events",gXvec,0);
  histNobs->Draw();

  //
  canvas->Modified();
  canvas->Update();
}

int main(int argc, char *argv[]) {
  if (argc>1) {
    if (loadData(std::string(argv[1]))) {
      TApplication theApp("App", &argc, argv);
      plotData();
      theApp.Run(kTRUE);
    }
  }
}
