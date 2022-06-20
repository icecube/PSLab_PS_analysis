
{
  // Set up plot style
  gROOT->SetStyle("Plain");     // Clean plotting style
  gStyle->SetOptStat(0);        // remove stats box from upper right corner
  gStyle->SetHistLineWidth(2);  // applies to all subsequent histograms


  double maxAngDeg = 10.;  
  int nBins = 100;

  TCut cut0 = "OneWeight*pow(mcEn,-2) * (panZd>90)";
  TH1D *hpsf0 = PointSpreadFunction(tree, cut0, "panZr", "panAr", 
				    maxAngDeg, nBins, "hpsf0");

  TCut cut1 = cut0 * "pbf_status==0 && pfSigmaDeg<3";
  TH1D *hpsf1 = PointSpreadFunction(tree, cut1, "panZr", "panAr", 
				    maxAngDeg, nBins, "hpsf1");

  TCut cut2 = cut1 * "pandel_rlogl < 8";
  TH1D *hpsf2 = PointSpreadFunction(tree, cut2, "panZr", "panAr", 
				    maxAngDeg, nBins, "hpsf2");

  
  // coordinates of legend x0,y0,x1,y1; e.g. whole pad would be (0,0,1,1)
  TLegend *legend = new TLegend(.4,.15,.88,.4,"IC22 Point Spread Function");  

  legend.AddEntry(hpsf0,"Trigger Level");
  legend.AddEntry(hpsf1,"Sigma < 3#circ cut");
  legend.AddEntry(hpsf2,"Sigma < 3#circ && reduced LLH < 8 cut");


  // The rest is automatic ...

  hpsf0->SetLineColor(kBlack);
  hpsf1->SetLineColor(kBlue);
  hpsf2->SetLineColor(kRed);



  // REGULAR (NOT INTEGRATED, NOT NORMALIZED)


  new TCanvas("psf_a","psf_a",20,20,640,480);
  hpsf0->DrawCopy();
  hpsf1->DrawCopy("same");
  hpsf2->DrawCopy("same");


  // INTEGRATED, BUT UN-NORMALIZED


  new TCanvas("psf_b","psf_b",40,40,640,480);

  gPad->SetGrid();              // x and y grid lines

  HistIntegrate(hpsf0)->Draw();
  HistIntegrate(hpsf1)->Draw("same");
  HistIntegrate(hpsf2)->Draw("same");

  legend->Draw();
  legend->SetFillColor(0);


  // INTEGRATED AND NORMALIZED

  new TCanvas("psf_c","psf_c",60,60,640, 480);

  gPad->SetGrid();              // x and y grid lines

  TH1 *hcopy = HistIntegrate( HistNormalize(hpsf0) );
  hcopy->Draw();
  // Can force y-axis range by controlling the first (main) histogram on pad
  hcopy->SetMaximum(1.);

  HistIntegrate( HistNormalize(hpsf1) )->Draw("same");
  HistIntegrate( HistNormalize(hpsf2) )->Draw("same");

  legend->Draw();
  legend->SetFillColor(0);

}
