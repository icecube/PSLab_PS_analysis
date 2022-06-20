
{
  gStyle->SetHistLineWidth(2);  // applies to all subsequent histograms
  

  // DEFINE THESE ELSEWHERE BEFORE RUNNING
  //  TCut cut;
  //  TTree *tree;
  //  TString titleString;

  TCut cut = "mDelAng<2.";


  TString titleString = "Solid-Angle-Averaged Effective Area:   #Delta#Psi < 2#circ";

  bool smoothing = true;  // "adjust" bins with low statistics


 if (1) {  // TURN OFF TO RE-RUN PLOTS WITHOUT RECALCULATING

  EffectiveArea eff;

  eff.SetEnergyRange(28,2,9);  // nBins, logEmin, logEmax

  TH1D *hArray[100];
  double zenMin[100];
  double zenMax[100];
  int color[100];
  int marker[100];
  int style[100];

  // always
  tree->SetAlias("NuMuFluxScaleFactor", "1");

  int s=0;
  if (1) {
    style[s] = 1; tree->SetAlias("NuTauFluxScaleFactor","1");  // turn tau on
    zenMin[s] = 0; zenMax[s] = 30; color[s] = kViolet; marker[s] = 23;
    hArray[s] = eff.Calculate(tree, zenMin[s],zenMax[s],cut,"h"+TStringify(s));
    ++s;

    style[s] = 2; tree->SetAlias("NuTauFluxScaleFactor","0");  // turn tau off
    zenMin[s] = 0; zenMax[s] = 30; color[s] = kViolet; marker[s] = 23;
    hArray[s] = eff.Calculate(tree, zenMin[s],zenMax[s],cut,"h"+TStringify(s));
    ++s;


    style[s] = 1; tree->SetAlias("NuTauFluxScaleFactor","1");  // turn tau on
    zenMin[s] = 30; zenMax[s] = 60; color[s] = kCyan+1; marker[s] = 24;
    hArray[s] = eff.Calculate(tree, zenMin[s],zenMax[s],cut,"h"+TStringify(s));
    ++s;

    style[s] = 2; tree->SetAlias("NuTauFluxScaleFactor","0");  // turn tau off
    zenMin[s] = 30; zenMax[s] = 60; color[s] = kCyan+1; marker[s] = 24;
    hArray[s] = eff.Calculate(tree, zenMin[s],zenMax[s],cut,"h"+TStringify(s));
    ++s;


    style[s] = 1; tree->SetAlias("NuTauFluxScaleFactor","1");  // turn tau on
    zenMin[s] = 60; zenMax[s] = 90; color[s] = kOrange-6; marker[s] = 25;
    hArray[s] = eff.Calculate(tree, zenMin[s],zenMax[s],cut,"h"+TStringify(s));
    ++s;

    style[s] = 2; tree->SetAlias("NuTauFluxScaleFactor","0");  // turn tau off
    zenMin[s] = 60; zenMax[s] = 90; color[s] = kOrange-6; marker[s] = 25;
    hArray[s] = eff.Calculate(tree, zenMin[s],zenMax[s],cut,"h"+TStringify(s));
    ++s;


    style[s] = 1; tree->SetAlias("NuTauFluxScaleFactor","1");  // turn tau on
    zenMin[s] = 90; zenMax[s] = 120;  color[s] = kRed; marker[s] = 20;
    hArray[s] = eff.Calculate(tree, zenMin[s],zenMax[s],cut,"h"+TStringify(s));
    ++s;

    style[s] = 2; tree->SetAlias("NuTauFluxScaleFactor","0");  // turn tau off
    zenMin[s] = 90; zenMax[s] = 120;  color[s] = kRed; marker[s] = 20;
    hArray[s] = eff.Calculate(tree, zenMin[s],zenMax[s],cut,"h"+TStringify(s));
    ++s;


    style[s] = 1; tree->SetAlias("NuTauFluxScaleFactor","1");  // turn tau on
    zenMin[s] = 120; zenMax[s] = 150; color[s] = kGreen+2; marker[s] = 21;
    hArray[s] = eff.Calculate(tree, zenMin[s],zenMax[s],cut,"h"+TStringify(s));
    ++s;

    style[s] = 2; tree->SetAlias("NuTauFluxScaleFactor","0");  // turn tau off
    zenMin[s] = 120; zenMax[s] = 150; color[s] = kGreen+2; marker[s] = 21;
    hArray[s] = eff.Calculate(tree, zenMin[s],zenMax[s],cut,"h"+TStringify(s));
    ++s;


    style[s] = 1; tree->SetAlias("NuTauFluxScaleFactor","1");  // turn tau on
    zenMin[s] = 150; zenMax[s] = 180; color[s] = kBlue; marker[s] = 22;
    hArray[s] = eff.Calculate(tree, zenMin[s],zenMax[s],cut,"h"+TStringify(s));
    ++s;

    style[s] = 2; tree->SetAlias("NuTauFluxScaleFactor","0");  // turn tau off
    zenMin[s] = 150; zenMax[s] = 180; color[s] = kBlue; marker[s] = 22;
    hArray[s] = eff.Calculate(tree, zenMin[s],zenMax[s],cut,"h"+TStringify(s));
    if (smoothing) {
      cout << "MANUAL SMOOTHING IN PROGRESS !!!...\n";
      hArray[s]->SetBinContent(24,0);
      hArray[s]->SetBinContent(17,6);
    }
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
    hArray[s]->SetLineStyle(style[s]);
    hArray[s]->SetMarkerColor(color[s]);
    hArray[s]->SetMarkerStyle(marker[s]);
    if (s%2==0) {
      legend.AddEntry(hArray[s],"zenith range ("+
		      TStringify(zenMin[s])+"#circ, "+
		      TStringify(zenMax[s])+"#circ)");
    }
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
		      "Phi_{#nu_{#mu}+#bar{#nu}_{#mu}} = Phi_{#nu_{#tau}+#bar{#nu}_{#tau}} Effective Area [m^{2}]");  
  hArray[0]->SetTitle(titleString);

  hArray[0]->GetXaxis()->CenterTitle(1);
  hArray[0]->GetYaxis()->CenterTitle(1);
  hArray[0]->GetYaxis()->SetTitleOffset(1.2);


  gPad->SetRightMargin(0.05);
  gPad->SetTopMargin(0.05);
  //  hArray[0]->SetTitle("");
  hArray[0]->GetXaxis()->SetNdivisions(-407);
}
