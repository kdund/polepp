{
  gROOT->Reset();
  //
  TFile f("confbelt.root");
  //
  TTree *tree = f.Get("name");

  //
  // var000 - signal hyp
  // var001 - N1
  // var002 - N2
  // var003 - P(s|N1,N2)
  //
  //
  //
  //
  tree->Draw("var000:var001:var002:var003","","goff");
  Long64_t nent = tree->GetEntries();
  Double_t smin = TMath::MinElement(nent, tree->GetV1());
  Double_t smax = TMath::MaxElement(nent, tree->GetV1());
  Double_t nmin = TMath::MinElement(nent, tree->GetV2());
  Double_t nmax = TMath::MaxElement(nent, tree->GetV3());
  Double_t s,w;
  Double_t n,n1,n2;
  Int_t nb;
  //
  TH2F *hist = new TH2F("confbelt","confbelt",int(nmax-nmin)+1,nmin,nmax,nent,smin,smax);

  for (int i=0; i<nent; i++) {
    s  = *(tree->GetV1()+i);
    n1 = *(tree->GetV2()+i);
    n2 = *(tree->GetV3()+i);
    w  = *(tree->GetV4()+i);
    nb = int(n2-n1);
    n = n1;
    for (int in=0; in<nb; in++) {
      hist->Fill(n,s,w);
      n += 1.0;
    }
  }

  hist->GetXaxis()->SetTitle("N");
  hist->GetYaxis()->SetTitle("s_{hyp}");
  hist->Draw("colz");

}
