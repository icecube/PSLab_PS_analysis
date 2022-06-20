// This macro plots several 1D histograms of the background
// and signal (for E^-2 and E^-3 spectra) for an already-loaded
// dataset (ark40) at the zenith angle zenithDeg. It used the entire
// zenith band, and there is a tPaveText which is set by hand at the 
// end to specify the bounds of the energy zenith band.

// A figure like this was included in the IC40 time-independent paper.

{

  figpaper1();

  double zenithDeg = 100.;
  
  ZenithEnergyProb * z2 = ark40.eProb;
  int db = z2->GetZenDegBand(zenithDeg);
  
  cout << db << endl;
  
  SimpleEnergyProb s1 = z2->GetSimpleEnergyProb(db);
  
  TH1D hbkg = s1.GetHistProbBkg();
  const TH1D * hhbkg = &hbkg;
  
  hbkg.GetXaxis()->SetTitle("log_{10} (Reconstructed Energy (GeV))");
  hbkg.GetYaxis()->SetTitle("Probability Density");
    
  TH1D * h2 = s1.GetHistProbGamma(2);
  TH1D * h3 = s1.GetHistProbGamma(3);
  
  hbkg.SetLineWidth(3);
  h2->SetLineWidth(2);
  h3->SetLineWidth(2);
  
  h3->SetLineStyle(2);
  
  h2->SetLineColor(9);
  h3->SetLineColor(4);
  
  //h2->Divide(hhbkg);
  //h3->Divide(hhbkg);
  
  hbkg->Draw();
  h2->Draw("same");
  h3->Draw("same");
  
  leg = new TLegend(0.4,0.6,0.8,0.8);
  leg->AddEntry(&hbkg,"Data","l");
  leg->AddEntry(h3,"E^{-3} Signal","l");
  leg->AddEntry(h2,"E^{-2} Signal","l");
  //leg->AddEntry(hbkg,"Data");
  leg->SetBorderSize(0);
  
  leg->Draw();
  
  TPaveText *pt = new TPaveText(2,0.2,4,2.6);
  pt->AddText("10#circ < #delta < +32#circ");
  pt->SetBorderSize(0);
  pt->Draw();
  
}  
