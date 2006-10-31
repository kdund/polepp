{
  gROOT->Reset();
  //
  TFile f("construct.root");
  //
  TTree *tree = f.Get("name");

  tree->Draw("var000:var001:var002","","goff");
  Long64_t nent = tree->GetEntries();
  Double_t smin = TMath::MinElement(nent, tree->GetV1());
  Double_t smax = TMath::MaxElement(nent, tree->GetV1());
  Double_t nmin = TMath::MinElement(nent, tree->GetV2());
  Double_t nmax = TMath::MaxElement(nent, tree->GetV2());
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
