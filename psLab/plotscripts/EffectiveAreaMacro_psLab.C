
{

gROOT->ProcessLine(".x $LAB_MAIN_DIR/llhTimeDep/loadlibs.C");

figpaper1();

bool upgoing = kTRUE;
bool downgoing = kFALSE;

gROOT->ProcessLine(".L standAlone/FastEffectiveArea.C");

gROOT->ProcessLine(".L $LAB_MAIN_DIR/macro_llh/ic59/TreeLoader_IC59_Final.C"); // For IC59

gROOT->ProcessLine(".L $LAB_MAIN_DIR/macro_llh/ic59/TreeLoader_IC40Full_CutA6v1_final.C"); // For IC40

TTree * tree = LoadTree_IC59_nugen_numu_6471_small(); // For IC59 Bugfixed
//TTree * tree1 = LoadTree_IC59_nugen_numu_4175_small(); // For IC59
TTree * tree1 = LoadTree_IC40Full_CutA6v1_nugen3311_small(); // For IC40


double mcTotalGeneratedEvents  = 9600*5000.;
double mcTotalGeneratedEvents1 = 9300*5000.;

//tree->SetAlias("mcTotalGeneratedEvents", "9600*5000. * (NCh/NCh)");
//tree1->SetAlias("mcTotalGeneratedEvents", "9300*5000. * (NCh/NCh)");

TCut selection = "mDelAng < 3";

//FastEffectiveArea(tree,0,180,selection,"hFastEffAreaAll");

//FastEffectiveArea(tree,0,90,selection,"hFastEffAreaDown");

TLegend *legend = new TLegend(.4, .15, .88, .4);

if (downgoing) {

  FastEffectiveArea(tree,selection,0, 30,40,1,9,"hFastEffAreaD1");
  FastEffectiveArea(tree,selection,30,60,40,1,9,"hFastEffAreaD2");
  FastEffectiveArea(tree,selection,60,90,40,1,9,"hFastEffAreaD3");

  FastEffectiveArea(tree1,selection,0, 30,40,1,9,"hFastEffAreaD11");
  FastEffectiveArea(tree1,selection,30,60,40,1,9,"hFastEffAreaD21");
  FastEffectiveArea(tree1,selection,60,90,40,1,9,"hFastEffAreaD31");

  //hFastEffAreaDown->SetLineColor(5);
  hFastEffAreaD1->SetLineColor(6);
  hFastEffAreaD2->SetLineColor(7);
  hFastEffAreaD3->SetLineColor(42);

  hFastEffAreaD11->SetLineColor(6);
  hFastEffAreaD21->SetLineColor(7);
  hFastEffAreaD31->SetLineColor(42);

  hFastEffAreaD11->SetLineStyle(2);
  hFastEffAreaD21->SetLineStyle(2);
  hFastEffAreaD31->SetLineStyle(2);


  //hFastEffAreaDown->SetLineWidth(2);
  hFastEffAreaD1->SetLineWidth(2);
  hFastEffAreaD2->SetLineWidth(2);
  hFastEffAreaD3->SetLineWidth(2);
  
  hFastEffAreaD11->SetLineWidth(2);
  hFastEffAreaD21->SetLineWidth(2);
  hFastEffAreaD31->SetLineWidth(2);
  
  //legend.AddEntry(hFastEffAreaDown,"zenith range (0#circ, 90#circ)");
  legend.AddEntry(hFastEffAreaD1,"#delta = (-90#circ, -60#circ)");
  legend.AddEntry(hFastEffAreaD2,"#delta = (-60#circ, -30#circ)");
  legend.AddEntry(hFastEffAreaD3,"#delta = (-30#circ, 0#circ)");
  
}

if (upgoing) {

  //FastEffectiveArea(tree,selection, 90,180,selection,"hFastEffAreaAll");
  FastEffectiveArea(tree,selection, 90,120,40,1,9,"hFastEffAreaB1");
  FastEffectiveArea(tree,selection,120,150,40,1,9,"hFastEffAreaB2");
  FastEffectiveArea(tree,selection,150,180,40,1,9,"hFastEffAreaB3");

  //FastEffectiveArea(tree1,selection, 90,180,selection,"hFastEffAreaAll1");
  FastEffectiveArea(tree1,selection, 90,120,40,1,9,"hFastEffAreaB11");
  FastEffectiveArea(tree1,selection,120,150,40,1,9,"hFastEffAreaB21");
  FastEffectiveArea(tree1,selection,150,180,40,1,9,"hFastEffAreaB31");

  //hFastEffAreaAll->SetLineColor(1);
  hFastEffAreaB1->SetLineColor(2);
  hFastEffAreaB2->SetLineColor(4);
  hFastEffAreaB3->SetLineColor(8);

  hFastEffAreaB11->SetLineColor(2);
  hFastEffAreaB21->SetLineColor(4);
  hFastEffAreaB31->SetLineColor(8);

  hFastEffAreaB11->SetLineStyle(2);
  hFastEffAreaB21->SetLineStyle(2);
  hFastEffAreaB31->SetLineStyle(2);

  //hFastEffAreaAll->SetLineWidth(2);
  hFastEffAreaB1->SetLineWidth(2);
  hFastEffAreaB2->SetLineWidth(2);
  hFastEffAreaB3->SetLineWidth(2);

  //hFastEffAreaAll->SetLineWidth(2);
  hFastEffAreaB11->SetLineWidth(2);
  hFastEffAreaB21->SetLineWidth(2);
  hFastEffAreaB31->SetLineWidth(2);
  
  //legend.AddEntry(hFastEffAreaAll,"#delta = (90#circ, 180#circ)");
  legend.AddEntry(hFastEffAreaB1,"#delta = (0#circ, 30#circ)");
  legend.AddEntry(hFastEffAreaB2,"#delta = (30#circ, 60#circ)");
  legend.AddEntry(hFastEffAreaB3,"#delta = (60#circ, 90#circ)");

}

TCanvas * cc = new TCanvas("cc","cc",800,600);
//cc->Divide(2,1);

TH1D * hRatio1 = new TH1D("hRatio1","hRatio1",40,1,9);
TH1D * hRatio2 = new TH1D("hRatio2","hRatio2",40,1,9);
TH1D * hRatio3 = new TH1D("hRatio3","hRatio3",40,1,9);


if (upgoing) {

  hFastEffAreaB1->GetYaxis()->SetTitle("Effective Area (m^{2})");
  hFastEffAreaB1->GetXaxis()->SetTitle("log_{10} (Energy / GeV)");

  hFastEffAreaB1->Draw();
  //hFastEffAreaAll->Draw("same");
  hFastEffAreaB2->Draw("same");
  hFastEffAreaB3->Draw("same");

  hFastEffAreaB21->Draw("same");
  hFastEffAreaB31->Draw("same");
  hFastEffAreaB11->Draw("same");

  
}

if (downgoing) {
  
  if (!upgoing) { 
    hFastEffAreaD3->GetYaxis()->SetTitle("Effective Area (m^{2})");
    hFastEffAreaD3->GetXaxis()->SetTitle("log_{10} (Energy / GeV)");
    hFastEffAreaD3->Draw("");
  }
  else { hFastEffAreaD1->Draw("same"); }
  hFastEffAreaD2->Draw("same"); 
  hFastEffAreaD1->Draw("same");

  hFastEffAreaD11->Draw("same");
  hFastEffAreaD21->Draw("same"); 
  hFastEffAreaD31->Draw("same");
  
  //hFastEffAreaAll->Draw("same");
  
}

legend->Draw();
cc->SetLogy();
cc->SetGrid(1,1);

}
