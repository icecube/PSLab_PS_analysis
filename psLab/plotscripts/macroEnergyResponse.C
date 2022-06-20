
{
  gStyle->SetHistLineWidth(2);  // applies to all subsequent histograms

  // DEFINE THESE ELSEWHERE BEFORE RUNNING
  TTree *tree;
  TCut cut;


  TLegend *legend = new TLegend(0.5,0.65,0.9,0.9);

  TCanvas *can = new TCanvas("canNuEREsponse","NuE Response");
  can->cd();


  TCut belowHorizon = "mcPrimary_Zenith_rad*TMath::RadToDeg() > 90.";
  TCut aboveHorizon = "mcPrimary_Zenith_rad*TMath::RadToDeg() < 90.";

  double binsPerDecade = 10.;

  tree->Draw("log10(mcPrimary_Energy_GeV)>>hAtm(70,2,9)",
	     belowHorizon*cut*"BartolFluxWeightForOneFile/1000.");
  double sum = hAtm->GetSum();
  hAtm->Scale(binsPerDecade/sum);
  hAtm->SetLineColor(kGreen+2);
  hAtm->SetLineStyle(9);
  hAtm->Smooth(2);
  legend->AddEntry(hAtm,"Atmospheric (up-going)");
  hAtm->SetTitle(";log_{10} (E_{#nu} / GeV );dP/d[log_{10}(E/GeV)]");

  tree->Draw("log10(mcPrimary_Energy_GeV)>>hEm2(70,2,9)",
	     belowHorizon*cut*"pow(mcEn,-2)*mcOneWeight","same");
  double sum = hEm2->GetSum();
  hEm2->Scale(binsPerDecade/sum);
  hEm2->SetLineColor(kBlue+2);
  hEm2->Smooth(5);
  legend->AddEntry(hEm2,"E^{-2} (up-going)");

  tree->Draw("log10(mcPrimary_Energy_GeV)>>hEm2down(70,2,9)",
	     aboveHorizon*cut*"pow(mcEn,-2)*mcOneWeight","same");
  double sum = hEm2down->GetSum();
  hEm2down->Scale(binsPerDecade/sum);
  hEm2down->SetLineColor(kViolet+1);
  hEm2down->SetLineStyle(2);
  hEm2down->Smooth(20);
  legend->AddEntry(hEm2down,"E^{-2} (down-going)");

  tree->Draw("log10(mcPrimary_Energy_GeV)>>hEm1_5(70,2,9)",
	     aboveHorizon*cut*"pow(mcEn,-1.5)*mcOneWeight","same");
  double sum = hEm1_5->GetSum();
  hEm1_5->Scale(binsPerDecade/sum);
  hEm1_5->SetLineColor(kRed+1);
  hEm1_5->Smooth(20);
  hEm1_5->SetLineStyle(3);
  hEm1_5->SetLineWidth(3);
  legend->AddEntry(hEm1_5,"E^{-1.5} (down-going)");
  
  legend->SetFillColor(0);
  legend->Draw();

  hAtm->SetTitleOffset(.9,"Y");
  hAtm->GetYaxis()->CenterTitle(true);
  hAtm->GetXaxis()->CenterTitle(true);
  gPad->SetGrid();
  gPad->SetTicks();

}
