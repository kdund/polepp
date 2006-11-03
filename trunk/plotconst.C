{
  gROOT->Reset();
  //
  TTree *tree = new TTree("constructTree", "construct");
  tree->ReadFile("construct.dat");
  tree->Draw("shyp:N:RL","","goff");
  Long64_t nent = tree->GetEntries();
  if (nent==0) {
    std::cout << "Error: no entries in tree! Check input file (construct.dat)" << std::endl;
    return;
  }
  Double_t smin = TMath::MinElement(nent, tree->GetV1());
  Double_t smax = TMath::MaxElement(nent, tree->GetV1());
  Double_t nmin = TMath::MinElement(nent, tree->GetV2());
  Double_t nmax = TMath::MaxElement(nent, tree->GetV2());
  Double_t wmin = TMath::MinElement(nent, tree->GetV3());
  Double_t wmax = TMath::MaxElement(nent, tree->GetV3());
  Double_t s,n,w,sp;
  Double_t ds;
  //
  Int_t i=0;
  sp = *(tree->GetV1()+0);
  while ((ds==0) && (i<nent)) {
    s = *(tree->GetV1()+i);
    ds = s-sp;
    sp = s;
    i++;
  }
  Int_t    nsbins = int((smax-smin)/ds)+1;
  TH2F *hist = new TH2F("construct","construct",int(nmax-nmin)+1,nmin,nmax,nsbins,smin,smax);
  hist->SetMinimum(wmin);
  hist->SetMaximum(wmax);
  for (int i=0; i<nent; i++) {
    s = *(tree->GetV1()+i);
    n = *(tree->GetV2()+i);
    w = *(tree->GetV3()+i);
    hist->Fill(n,s,w);
  }

  hist->GetXaxis()->SetTitle("N");
  hist->GetYaxis()->SetTitle("s_{hyp}");
  hist->Draw("colz");

}
