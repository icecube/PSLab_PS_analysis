{

  gROOT->ProcessLine(".x ../llhTimeDep/loadlibs.C");
  gROOT->ProcessLine(".x loadlibs.C");

  figpaper1();

  TChain * chMap = new TChain("tAllSky");
  //chMap->Add("/net/user/mfbaker/lab/IC40/flares_g/maps/tAllSkyBatch_UnblindFine*.root");
  //chMap->Add("/net/user/mfbaker/psLab/trunk/macro_llh/ic59/maps/flare/tAllSkyBatch0.root");
  chMap->Add("/net/user/mfbaker/psLab/trunk/macro_llh/ic59/maps/flare/tAllSkyBatch0_unb_new.root");
  
  int var=0; //0 for estp
             //1 for mean
             //2 for sigma

  if (var==0) chMap->Draw("dec:ra>>h1(720,0,360,340,-85,85)","-log10(estp)*(estp<0.3)","goff");
  if (var==1) chMap->Draw("dec:ra>>h1(720,0,360,340,-85,85)","(mean-54900)*(estp<0.3)","goff");
  if (var==2) chMap->Draw("dec:ra>>h1(720,0,360,340,-85,85)","(sigma)*(estp<0.3)","goff");
  if (var==3) chMap->Draw("dec:ra>>h1(720,0,360,340,-85,85)","(ns)*(estp<0.3)","goff");
  if (var==4) chMap->Draw("dec:ra>>h1(720,0,360,340,-85,85)","(gamma)*(estp<0.3)","goff");

  TH2D * hInput = (TH2D*) h1;

  Projection * pro = new HammerAitoffProjection();

//  gROOT->ProcessLine(".x SimpleSkyPlot.C"); 
  
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  TCanvas *can = new TCanvas("can","can",20,20,1550,775);
  // TCanvas *can = new TCanvas("can","can",20,20,1000,500);
  can->SetLeftMargin(0.05);
  can->SetTopMargin(0.06);
  can->SetBottomMargin(0.06);
  can->SetRightMargin(0.15);
  can->SetFrameLineColor(0); // turns off histogram frame

  SetRootPalette(1);
  if (var==0) { CreatePalette(5); }
  else { gStyle->SetPalette(1); }
  //SetRootPalette(1);

  TH2D *hProjection = ProjectHistogram(hInput, pro);

  if (var==0) hProjection->SetMinimum(-log10(0.5));
  //if (var==1)  hProjection->SetMinimum(570.);
  if (var==2)  hProjection->SetMinimum(1e-6);
//  hProjection->SetMinimum(1.);
  TAxis *zaxis = hProjection->GetZaxis();
  zaxis->CenterTitle();
//  zaxis->SetLabelOffset(0.004);
  if (var==0)  zaxis->SetTitle("-log_{10} p");
  if (var==2)  zaxis->SetTitle("#sigma_{T} (days)");
  if (var==3)  zaxis->SetTitle("# source events");
  if (var==4)  zaxis->SetTitle("Spectral Index");
  if (var==4)  hProjection->SetMinimum(1.);
  if (var==1)  { 
    zaxis->SetTitle("T_{o} (MJD-54900)");
    //double tmin = h1->GetMinimum();
    //hProjection->SetMinimum(0);
  }
    
//  zaxis->SetTitle("");
  zaxis->SetTitleOffset(0.7);
  zaxis->SetTitleSize(0.05);
  zaxis->SetTitleFont(42);
  zaxis->SetLabelFont(42);  

  hProjection->Draw("AHcolz");  // AH removes axes
  gPad->Update();
  TPaletteAxis* paxis = GetHistogramPalette(hProjection);
  paxis->SetX1NDC(0.887);
  paxis->SetX2NDC(0.927);


  // Model for graph style for longitude and latitude lines

  TGraph *gModel = new TGraph(2000);
  gModel->SetMarkerStyle(1);
  gModel->SetMarkerColor(kGray+2);
  gModel->SetLineColor(kGray+2);
  char *gOpt = "P";

  // Lines of Longitude

  for (int iLon=0; iLon <= 6; ++iLon) {
    double lon = 60. * iLon;
    ProjectGraphLon(lon, pro, -85., 85., gModel)->Draw(gOpt);
  }

  // Lines of Latitude

  for (int iLat=1; iLat <= 5; ++iLat) {
    double lat = -90. + 30. * iLat;
    ProjectGraphLat(lat, pro, 0., 360., gModel)->Draw(gOpt);
  }
  ProjectGraphLat(-85, pro, 0., 360., gModel)->Draw(gOpt);
  ProjectGraphLat(+85, pro, 0., 360., gModel)->Draw(gOpt);

  // Galactic Plane

  TGraph *galPlane = ProjectGalacticPlane(pro);
  galPlane->SetMarkerStyle(1);
  galPlane->Draw("P");

  // Coordinate Labels

  double xMin = pro->GetXMin();
  double xMax = pro->GetXMax();
  double xLeft  = xMin+(xMax-xMin)*1e-7;
  double xRight = xMax-(xMax-xMin)*1e-7;
  int hourLeft  = int(pro->Lon(xLeft, 0)/15. + 0.5);
  int hourRight = int(pro->Lon(xRight,0)/15. + 0.5);
  
  const char hL[10], hR[10];
  sprintf(hL,"%ih",hourLeft );
  sprintf(hR,"%ih",hourRight);
  
  TString textLeft =  TString(hL);
  TString textRight = TString(hR);
  
//  TString textLeft = hourLeft+TString("h");
//  TString textRight = hourRight+TString("h");

  TPaveLabel *p1 = new TPaveLabel(0.019,0.5-0.015,0.033,0.5+0.015,
				  textLeft,"NDC");
  p1->SetBorderSize(0);
  p1->SetTextFont(42);
  p1->SetLineColor(0);
  p1->SetFillColor(0);
  p1->SetTextSize(1.5);
  p1->Draw();
  gPad->Update();

  TPaveLabel *p2 = new TPaveLabel(*p1);
  p2->SetLabel(textRight);
  p2->SetTextFont(42);
  p2->SetX1NDC(0.86);
  p2->SetX2NDC(0.875);
  p2->Draw();
  gPad->Update();
  
  TPaveText tp30(-2.85,0.67,-2.50,0.82);
  TPaveText tp60(-1.85,1.23,-1.50,1.35);
  
  TPaveText tpm30(2.5,-0.82,2.85,-0.67);
  TPaveText tpm60(1.50,-1.35,1.85,-1.23);

  tp30.AddText("#delta = +30#circ");
  tp60.AddText("#delta = +60#circ");
  tpm30.AddText("#delta = -30#circ");
  tpm60.AddText("#delta = -60#circ");
  
  tp30.SetBorderSize(0);
  tp60.SetBorderSize(0);
  tpm30.SetBorderSize(0);
  tpm60.SetBorderSize(0);

  tp30.SetTextFont(42);
  tp60.SetTextFont(42);
  tpm30.SetTextFont(42);
  tpm60.SetTextFont(42);
  
  tp30.SetFillColor(0);
  tp60.SetFillColor(0);
  tpm30.SetFillColor(0);
  tpm60.SetFillColor(0);

  tp30.Draw();
  tp60.Draw();
  tpm30.Draw();
  tpm60.Draw();
  
 
  


}
