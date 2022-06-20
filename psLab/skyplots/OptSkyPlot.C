{

  /* Version with many optional overlays. 
   / This is a helper script for macro_opt.C
  */

  TH2D *hInput;     // if already set, this will not override value
  Projection* pro;  // if already set, this will not override value
  bool OPT_PROBMAP;
  bool OPT_EVENTS;
  bool OPT_GRID;
  bool OPT_GALPLANE;
  if (!hInput) { cout << "Error: hInput not defined.\n"; return; }
  if (OPT_EVENTS && !gEvents) { cout << "Error: gEvents not defined.\n"; return; }
  if (!pro) { cout << "Error: pro not defined.\n"; return; }


  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  TCanvas *can = new TCanvas("can","can",20,20,1550,775);
  //TCanvas *can = new TCanvas("can","can",20,20,455,377);
  // TCanvas *can = new TCanvas("can","can",20,20,1000,500);
/*
  can->SetLeftMargin(0.00);
  can->SetTopMargin(0.00);
  can->SetBottomMargin(0.00);
  can->SetRightMargin(0.00);
*/
  can->SetLeftMargin(0.05);
  can->SetTopMargin(0.06);
  can->SetBottomMargin(0.06);
  can->SetRightMargin(0.15);
  can->SetFrameLineColor(0); // turns off histogram frame

  SetRootPalette(17);


  TH2D *hProjection = ProjectHistogram(hInput, pro);
  if (OPT_EVENTS) TGraph *gProjectionEvents = ProjectTGraph(gEvents, pro);

  hProjection->SetMinimum(-log10(0.5));
  //hProjection->SetMaximum(6.2);
  TAxis *zaxis = hProjection->GetZaxis();
  zaxis->CenterTitle();
  zaxis->SetTitle("-log_{10} p");
  zaxis->SetTitleOffset(0.7);
  zaxis->SetTitleSize(0.05);

  hProjection->Draw("AHcolz");  // AH removes axes
  //hProjection->Draw("AHcol");  // AH removes axes
  gPad->Update();
  TPaletteAxis* paxis = GetHistogramPalette(hProjection);
  paxis->SetX1NDC(0.887);
  paxis->SetX2NDC(0.927);

  // Events
  if (OPT_EVENTS) {
    if (OPT_BIG_EVENTS) {
      gProjectionEvents->SetMarkerStyle(21);
      gProjectionEvents->SetMarkerSize(0.3);
    }
    gProjectionEvents->Draw("P");
  }


  // Model for graph style for longitude and latitude lines

  TGraph *gModel = new TGraph(2000);
  gModel->SetMarkerStyle(1);
  gModel->SetMarkerColor(kGray+2);
  gModel->SetLineColor(kGray+2);
  char *gOpt = "P";

  if (OPT_GRID) {
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

  }

  // Galactic Plane

  if (OPT_GALPLANE) {
    TGraph *galPlane = ProjectGalacticPlane(pro);
    galPlane->SetMarkerStyle(20);
    galPlane->SetMarkerSize(0.15);
    galPlane->Draw("P");
  }

  if (OPT_COORD) {
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

    // New Dec labels:
    TLatex *p3 = new TLatex(2.47,0.71,"+30#circ");
    p3->SetTextSize(0.03);
    p3->Draw();
    gPad->Update();

    TLatex *p4 = new TLatex(1.38,1.24,"+60#circ");
    p4->SetTextSize(0.03);
    p4->Draw();
    gPad->Update();

    TLatex *p5 = new TLatex(2.47,-0.76,"-30#circ");
    p5->SetTextSize(0.03);
    p5->Draw();
    gPad->Update();

    TLatex *p6 = new TLatex(1.38,-1.31,"-60#circ");
    p6->SetTextSize(0.03);
    p6->Draw();
    gPad->Update();

  }
  if (OPT_GAL_SHIFT_COORD) {
    // IMPERFECT HACK for gal lon and lat:
    // Lines of Gal Longitude
    TLine *line;
    for (int iLon=0; iLon <= 6; ++iLon) {
      double lon = 60. * iLon;
      line = new TLine(lon,-90,lon,90);
      line->DrawClone();
    }

    // Lines of Gal Latitude

    for (int iLat=1; iLat <= 5; ++iLat) {
      double lat = -90. + 30. * iLat;
      line = new TLine(0,lat,360,lat);
      line->DrawClone();
    }
    line = new TLine(0,-90,360,-90);
    line->DrawClone();
    line = new TLine(0,90,360,90);
    line->DrawClone();
    //ProjectGraphLat(-85, pro, 0., 360., gModel)->Draw(gOpt);
    //ProjectGraphLat(+85, pro, 0., 360., gModel)->Draw(gOpt);

    // Coordinate Labels

    // Gal coords
    TLatex *p7 = new TLatex(-18,-1.5,"l=-180#circ");
    p7->SetTextSize(0.025);
    p7->Draw();
    gPad->Update();

    TLatex *p8 = new TLatex(361.,-1.5,"l=180#circ");
    p8->SetTextSize(0.025);
    p8->Draw();
    gPad->Update();

    TLatex *p9 = new TLatex(174,92,"b=90#circ");
    p9->SetTextSize(0.025);
    p9->Draw();
    gPad->Update();

    TLatex *p10 = new TLatex(174,-95,"b=-90#circ");
    p10->SetTextSize(0.025);
    p10->Draw();
    gPad->Update();

  }

  if (!OPT_PROBMAP) hProjection->Delete(); // Now get rid of the prob map (only used for setup)
  gPad->Update();

}
