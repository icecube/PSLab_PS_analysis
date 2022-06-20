
{ 
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);


  double pThreshold = 1e-4;

  int count = 0;

  for (int i=1; i<=100; ++i) {
    char filename[500];
    sprintf(filename,"./results/ic22_trial_%08d_hLogpDec.root",i);

    TFile *f = new TFile(filename);
    if (count == 0) {
      TH2D *h = dynamic_cast<TH2D*> (f->Get("hLogpDec"));
      gROOT->cd();
      TH2D *hMain = dynamic_cast<TH2D*> (h->Clone("hMain"));
    } else {
      TH2D *h = dynamic_cast<TH2D*> (f->Get("hLogpDec"));
      gROOT->cd();
      hMain->Add(h);
    }
    f->Close();
    ++count;
  }

  new TCanvas("c1","c1");
  hMain->Draw("box");

  int nDecBins = hMain->GetNbinsY();
  int nTrials = hMain->GetSum() / nDecBins;
  cout << nTrials << " trials in each of " << nDecBins;
  cout << " declination bands.\n";



  TH1D *hSummary = new TH1D("hSummary","hSummary",
			    nDecBins, 
			    hMain->GetYaxis()->GetXmin(),
			    hMain->GetYaxis()->GetXmax());
			    
  int nLastBin = hMain->GetXaxis()->FindBin(log10(pThreshold));

  for (int i=1; i<=nDecBins; ++i) {
    int nTotal = 0;
    for (int j=0; j<nLastBin; ++j) {
      nTotal += hMain->GetBinContent(j,i);
    }

    double actualPValue = double(nTotal)/nTrials;
    hSummary->SetBinContent(i,actualPValue);
    if (nTotal>0) {
      hSummary->SetBinError(i, actualPValue / sqrt(nTotal) );
    }
  }

  new TCanvas("c2","c2");
  hSummary->Draw("e");
  gPad->SetLogy();
  TLine tline(hSummary->GetXaxis()->GetXmin(),pThreshold,
	      hSummary->GetXaxis()->GetXmax(),pThreshold);
  tline.SetLineColor(2);
  tline.Draw();

  hSummary->Draw("same");	      
}
