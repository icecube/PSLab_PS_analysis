
{
  gStyle->SetHistLineWidth(2);  // applies to all subsequent histograms
  

  // DEFINE THESE ELSEWHERE BEFORE RUNNING
  //  TCut cut;
  //  TTree *tree;
  //  TString titleString;

  TCut cut = "mDelAng<2.";


  TString titleString = "Solid-Angle-Averaged Effective Area:   #Delta#Psi < 2#circ";

  bool smoothing = false;  // "adjust" bins with low statistics


 if (1) {  // TURN OFF TO RE-RUN PLOTS WITHOUT RECALCULATING

  EffectiveArea eff;

  eff.SetEnergyRange(28,2,9);  // nBins, logEmin, logEmax

  TH1D *hArray[100];
  double zenMin[100];
  double zenMax[100];
  int color[100];
  int marker[100];

  int s=0;

  if (1) {
    zenMin[s] = 0; zenMax[s] = 30; color[s] = kViolet; marker[s] = 23;
    hArray[s] = eff.Calculate(tree, zenMin[s],zenMax[s],cut,"h"+TStringify(s));
    ++s;

    zenMin[s] = 30; zenMax[s] = 60; color[s] = kCyan+1; marker[s] = 24;
    hArray[s] = eff.Calculate(tree, zenMin[s],zenMax[s],cut,"h"+TStringify(s));
    ++s;

    zenMin[s] = 60; zenMax[s] = 90; color[s] = kOrange-6; marker[s] = 25;
    hArray[s] = eff.Calculate(tree, zenMin[s],zenMax[s],cut,"h"+TStringify(s));
    ++s;

    zenMin[s] = 90; zenMax[s] = 120;  color[s] = kRed; marker[s] = 20;
    hArray[s] = eff.Calculate(tree, zenMin[s],zenMax[s],cut,"h"+TStringify(s));
    ++s;

    zenMin[s] = 120; zenMax[s] = 150; color[s] = kGreen+2; marker[s] = 21;
    hArray[s] = eff.Calculate(tree, zenMin[s],zenMax[s],cut,"h"+TStringify(s));
    ++s;

    zenMin[s] = 150; zenMax[s] = 180; color[s] = kBlue; marker[s] = 22;
    hArray[s] = eff.Calculate(tree, zenMin[s],zenMax[s],cut,"h"+TStringify(s));
    if (smoothing) {
      cout << "MANUAL SMOOTHING IN PROGRESS !!!...\n";
      hArray[s]->SetBinContent(24,0);
      hArray[s]->SetBinContent(17,6);
    }
    ++s;
  }

  if (0) {
    zenMin[s] = 0; zenMax[s] = 90; color[s] = kGreen+9; marker[s] = 27;
    hArray[s] = eff.Calculate(tree, zenMin[s],zenMax[s],cut,"h"+TStringify(s));
    ++s;

    zenMin[s] = 90; zenMax[s] = 180; color[s] = kBlack; marker[s] = 1;
    hArray[s] = eff.Calculate(tree, zenMin[s],zenMax[s],cut,"h"+TStringify(s));
    ++s;
  }


  if (0) {
    zenMin[s] = 0; zenMax[s] = 80; color[s] = kRed+2; marker[s] = 21;
    hArray[s] = eff.Calculate(tree, zenMin[s],zenMax[s],cut,"h"+TStringify(s));
    ++s;

    zenMin[s] = 80; zenMax[s] = 100; color[s] = kGreen+2; marker[s] = 25;
    hArray[s] = eff.Calculate(tree, zenMin[s],zenMax[s],cut,"h"+TStringify(s));
    ++s;
    
    zenMin[s] = 100; zenMax[s] = 120; color[s] = kBlue; marker[s] = 23;
    hArray[s] = eff.Calculate(tree, zenMin[s],zenMax[s],cut,"h"+TStringify(s));
    ++s;

    zenMin[s] = 140; zenMax[s] = 160; color[s] = kBlue+2; marker[s] = 23;
    hArray[s] = eff.Calculate(tree, zenMin[s],zenMax[s],cut,"h"+TStringify(s));
    ++s;
  }
   
  int sTotal = s;

 }


  // coordinates of legend x0,y0,x1,y1; e.g. whole pad would be (0,0,1,1)
  TLegend *legend = new TLegend(.12, .68, .48, .93);  
  legend->SetFillColor(0);

  //  legend->SetHeader(tree->GetTitle());

  for (int s=0; s<sTotal; ++s) {
    hArray[s]->SetLineColor(color[s]);
    hArray[s]->SetMarkerColor(color[s]);
    hArray[s]->SetMarkerStyle(marker[s]);
    legend.AddEntry(hArray[s],"zenith range ("+
		    TStringify(zenMin[s])+"#circ, "+
		    TStringify(zenMax[s])+"#circ)");
  }


  new TCanvas("EffArea","EffArea");

  gPad->SetLogy(1);              // Log scale on y axis
  gPad->SetGrid(1);              // x and y grid lines
  gPad->SetTicks(1);             // Tick marks on all sides


  hArray[0]->Draw("PH");
  for (int s=1; s<sTotal; ++s) { hArray[s]->Draw("PHsame"); }

  legend->Draw();

  // Can force y-axis range by controlling the first (main) histogram on pad
  hArray[0]->SetMaximum(3e4.);
  hArray[0]->SetMinimum(3e-4);

  hArray[0]->SetTitle("Average Effective Area;"
		      "log_{10} Primary Energy / GeV;m^{2}");  
  hArray[0]->SetTitle("Average Effective Area;"
		      "log_{10} ( E_{#nu} / GeV );"
		      "#nu_{#mu} + #bar{#nu}_{#mu} Effective Area [m^{2}]");  
  hArray[0]->SetTitle(titleString);

  hArray[0]->GetXaxis()->CenterTitle(1);
  hArray[0]->GetYaxis()->CenterTitle(1);
  hArray[0]->GetYaxis()->SetTitleOffset(1.2);


  gPad->SetRightMargin(0.05);
  gPad->SetTopMargin(0.05);
  //  hArray[0]->SetTitle("");
  hArray[0]->GetXaxis()->SetNdivisions(-407);
}
