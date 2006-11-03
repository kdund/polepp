{
  gROOT->Reset();
  //
  TFile f("confbelt.dat");
  //
  TTree *tree = new TTree("confbeltTree", "confbelt");
  tree->ReadFile("confbelt.dat");
  tree->Draw("shyp:N1:N2:P","","goff");

  Long64_t nent = tree->GetEntries();
  if (nent==0) {
    std::cout << "Error: no entries in tree! Check input file (confbelt.dat)" << std::endl;
    return;
  }
  Double_t smin = TMath::MinElement(nent, tree->GetV1());
  Double_t smax = TMath::MaxElement(nent, tree->GetV1());
  Double_t nmin = TMath::MinElement(nent, tree->GetV2());
  Double_t nmax = TMath::MaxElement(nent, tree->GetV3());
  Double_t wmin = TMath::MinElement(nent, tree->GetV4());
  Double_t wmax = TMath::MaxElement(nent, tree->GetV4());
  Double_t s,w;
  Double_t n,n1,n2;
  Int_t nb;
  //
  TH2F *hist = new TH2F("confbelt","confbelt",int(nmax-nmin)+1,nmin,nmax,nent,smin,smax);
  hist->SetMinimum(wmin);
  hist->SetMaximum(wmax);

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
