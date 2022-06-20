
{
  //  EXAMPLE SESSION:
  //  root [0] .L FastEffectiveArea.C
  //  root [1] new TFile("myfile.root")
  //  root [2] TTree *tree = myTreeFromFile
  //  root [3] .x macroEffectiveExample.C



  TCut cut = "";
  TString titleString = "Solid-Angle-Averaged Effective Area";

  int eBins = 28;
  double logEMin = 2;
  double logEMax = 9;


  // Set Up

  gROOT->SetStyle("Plain");
  gStyle->SetHistLineWidth(2);  // applies to all subsequent histograms
  gStyle->SetOptStat(0);

  // coordinates of legend x0,y0,x1,y1; e.g. whole pad would be (0,0,1,1)
  TLegend *legend = new TLegend(.15, .75, .55, .89);  
  legend->SetFillColor(0);

  TH1D *hArray[100];
  double zenMin[100];
  double zenMax[100];
  int color[100];
  int marker[100];


  // Specify Different Effective Area Zenith Ranges

  int s=0;
  zenMin[s] = 90; zenMax[s] = 120;  color[s] = kRed; marker[s] = 20;
  ++s;
  zenMin[s] = 120; zenMax[s] = 150; color[s] = kGreen+2; marker[s] = 21;
  ++s;
  zenMin[s] = 150; zenMax[s] = 180; color[s] = kBlue; marker[s] = 22;
  ++s;


  for (int i=0; i<s; ++i) {
    hArray[i] = FastEffectiveArea(tree, cut, zenMin[i], zenMax[i],
				  eBins, logEMin, logEMax,
				  "h"+TStringify(i));
    hArray[i]->SetLineColor(color[i]);
    hArray[i]->SetMarkerColor(color[i]);
    hArray[i]->SetMarkerStyle(marker[i]);
    legend.AddEntry(hArray[i],"zenith range ("+
		    TStringify(zenMin[i])+"#circ, "+
		    TStringify(zenMax[i])+"#circ)");
  }


  new TCanvas("EffArea","EffArea");

  gPad->SetLogy(1);              // Log scale on y axis
  gPad->SetGrid(1);              // x and y grid lines
  gPad->SetTicks(1);             // Tick marks on all sides


  hArray[0]->Draw("PH");
  for (int i=1; i<s; ++i) { hArray[i]->Draw("PHsame"); }
  legend->Draw();


  // Can force y-axis range by controlling the first (main) histogram on pad
  hArray[0]->SetMaximum(3e4.);
  hArray[0]->SetMinimum(3e-4);

  hArray[0]->SetTitle("Average Effective Area;"
		      "log_{10} ( E_{#nu} / GeV );"
		      "#nu_{#mu} + #bar{#nu}_{#mu} Effective Area [m^{2}]");  
  hArray[0]->SetTitle(titleString);

  hArray[0]->GetXaxis()->CenterTitle(1);
  hArray[0]->GetYaxis()->CenterTitle(1);
  hArray[0]->GetYaxis()->SetTitleOffset(1.2);

  gPad->SetRightMargin(0.05);
  gPad->SetTopMargin(0.07);

  hArray[0]->GetXaxis()->SetNdivisions(-407);
}
