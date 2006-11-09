void plotcov() {

  gROOT->Reset();
  //
  TFile f("coverage.dat");
  //
  TTree *tree = new TTree("coverageTree", "coverage");
  tree->ReadFile("coverage.dat");
  Long64_t nent = tree->GetEntries();
  if (nent==0) {
    std::cout << "Error: no entries in tree! Check input file (coverage.dat)" << std::endl;
    return;
  }
  tree->Draw("shyp:effmean:bkgmean","","goff");

  Double_t smin = TMath::MinElement(nent, tree->GetV1());
  Double_t smax = TMath::MaxElement(nent, tree->GetV1());
  Double_t emin = TMath::MinElement(nent, tree->GetV2());
  Double_t emax = TMath::MaxElement(nent, tree->GetV2());
  Double_t bmin = TMath::MinElement(nent, tree->GetV3());
  Double_t bmax = TMath::MaxElement(nent, tree->GetV3());
  Double_t diffs = smax-smin;
  Double_t diffe = emax-emin;
  Double_t diffb = bmax-bmin;
  Bool_t scanshyp = (diffs>1e-15);
  Bool_t scaneff  = (diffe>1e-15);
  Bool_t scanbkg  = (diffb>1e-15);
  Double_t xmin,xmax,ymin,ymax,eymax;
  Double_t dx;

  //
  TVectorD vshyp(nent);
  TVectorD veffm(nent);
  TVectorD vbkgm(nent);
  TVectorD veffm(nent);
  TVectorD vcovm(nent);
  TVectorD vcovs(nent);
  TVectorD vnull(nent);
  //
  for (int i=0; i<nent; i++) {
    vnull(i) = 0.0;
    vshyp(i) = *(tree->GetV1()+i);
    veffm(i) = *(tree->GetV2()+i);
    vbkgm(i) = *(tree->GetV3()+i);
  }
  //
  tree->Draw("cov:coverr","","goff");
  for (int i=0; i<nent; i++) {
    vcovm(i) = *(tree->GetV1()+i);
    vcovs(i) = *(tree->GetV2()+i);
  }
  eymax = vcovs.Max();
  ymin  = vcovm.Min()-eymax;
  ymax  = vcovm.Max()+eymax;

  //
  TVectorD *vxuse=0;
  TString xtit;
  if (scanshyp) {
    dx = 1.1*(vshyp(1)-vshyp(0));
    xmin = TMath::Max(smin-dx,0.0);
    xmax = smax+dx;
    vxuse = &vshyp;
    xtit = "signal";
  } else if (scaneff) {
    dx = 1.1*(veffm(1)-veffm(0));
    xmin = TMath::Max(emin-dx,0.0);
    xmax = emax+dx;
    vxuse = &veffm;
    xtit = "efficiency";
  } else if (scanbkg) {
    dx = 1.1*(vbkgm(1)-vbkgm(0));
    xmin = TMath::Max(bmin-dx,0.0);
    xmax = bmax+dx;
    vxuse = &vbkgm;
    xtit = "background";
  } else {
    std::cout << "WARNING: None of the variables are scanned!" << std::endl;
    return;
  }
  TGraphErrors *graph = new TGraphErrors( *vxuse,vcovm,vnull,vcovs);


  TCanvas *canvas;
  canvas = new TCanvas("coveragecv","coverage",200,10,600,480);
  canvas->SetFillColor(0);
  canvas->SetFillStyle(0);
  canvas->SetGrid();
  canvas->GetFrame()->SetFillColor(21);
  canvas->GetFrame()->SetBorderSize(4);
  canvas->ToggleEventStatus();
  canvas->Range(-1.97799,0.915818,11.8912,1.01673);
  canvas->SetLeftMargin(0.142617);
  canvas->SetRightMargin(0.057047);
  canvas->SetTopMargin(0.0884956);
  canvas->SetBottomMargin(0.112832);
  gStyle->SetPaperSize(26,20);
  TH1F *frame = canvas->DrawFrame(xmin,ymin,xmax,ymax);
  //
  frame->SetXTitle(xtit);
  frame->SetYTitle("<1-\\alpha_{eff}>");
  frame->SetTitleOffset(0.70,"X");
  frame->SetTitleOffset(0.95,"Y");
  frame->SetTitleSize(0.06,"X");
  frame->SetTitleSize(0.06,"Y");

  canvas->cd();

  graph->SetMarkerStyle(22);
  graph->SetMarkerSize(1.1);
  graph->Draw("LP");
}
