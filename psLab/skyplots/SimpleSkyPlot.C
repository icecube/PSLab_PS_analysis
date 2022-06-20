
{
  TH2D *hInput;     // if already set, this will not override value
  Projection* pro;  // if already set, this will not override value
  if (!hInput) { cout << "Error: hInput not defined.\n"; return; }
  if (!pro) { cout << "Error: pro not defined.\n"; return; }


  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  TCanvas *can = new TCanvas("can","can",20,20,1550,775);
  // TCanvas *can = new TCanvas("can","can",20,20,1000,500);
  can->SetLeftMargin(0.05);
  can->SetTopMargin(0.06);
  can->SetBottomMargin(0.06);
  can->SetRightMargin(0.15);
  can->SetFrameLineColor(0); // turns off histogram frame

  SetRootPalette(17);
  

  TH2D *hProjection = ProjectHistogram(hInput, pro);

  hProjection->SetMinimum(-log10(0.5));
  TAxis *zaxis = hProjection->GetZaxis();
  zaxis->CenterTitle();
  zaxis->SetTitle("-log_{10} p");
  zaxis->SetTitleOffset(0.7);
  zaxis->SetTitleSize(0.05);

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
  int hourRight = int(pro->Lon(xRight,0)/15.+ 0.5);
  TString textLeft = hourLeft+TString("h");
  TString textRight = hourRight+TString("h");

  TPaveLabel *p1 = new TPaveLabel(0.019,0.5-0.015,0.033,0.5+0.015,
				  textLeft,"NDC");
  p1->SetBorderSize(0);
  p1->SetLineColor(0);
  p1->SetFillColor(0);
  p1->SetTextSize(1.5);
  p1->Draw();
  gPad->Update();

  TPaveLabel *p2 = new TPaveLabel(*p1);
  p2->SetLabel(textRight);
  p2->SetX1NDC(0.86);
  p2->SetX2NDC(0.875);
  p2->Draw();
  gPad->Update();
}
