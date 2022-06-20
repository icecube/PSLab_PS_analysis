
{
  // Must already be defined (e.g. in Efficiency_Setup.C)

  // TTree *treeRaw;
  // TTree *treeCut;
  // TCut weight;
  // TCut cutRaw;
  // TCut cutCut;


  double canFactor = 1.;
  TCanvas *can = new TCanvas("can","can",20,20,900*canFactor,550*canFactor);
  can->Divide(2,2,0.005,0.005);
  gStyle->SetNumberContours(255);
  gStyle->SetPalette(1);

  //  double b[6] = {35,2,9,40,-1,1}; // nbx, x0, x1, nby, y0, y1
  //  TString plotString = "cos(mcPrimary_Zenith_rad):log10(mcPrimary_Energy_GeV)";

  if (0) {
    double b[6] = {40,-1,1,35,2,9}; // nbx, x0, x1, nby, y0, y1
    TString plotStringX = "cos(mcPrimary_Zenith_rad)";
    TString plotStringY = "log10(mcPrimary_Energy_GeV)";
  }

  if (1) {
    double b[6] = {20,-1,1,25,1,3.5}; // nbx, x0, x1, nby, y0, y1
    TString plotStringX = "cos(mcPrimary_Zenith_rad)";
    TString plotStringY = "log10(NChan)";
  }

  if (0) {
    double b[6] = {40,-1,1,30,2,8}; // nbx, x0, x1, nby, y0, y1
    TString plotStringX = "cos(mcPrimary_Zenith_rad)";
    TString plotStringY = "log10(mmueEn)";
  }

  if (0) {
    double b[6] = {35,2,9,30,2,8}; // nbx, x0, x1, nby, y0, y1
    TString plotStringX = "log10(mcPrimary_Energy_GeV)";
    TString plotStringY = "log10(mmueEn)";
  }



  can->cd(2);
  TH2D hCut("hCut","Cut;"+plotStringX+";"+plotStringY,
	    b[0],b[1],b[2],b[3],b[4],b[5]);
  treeCut->Draw(plotStringY+":"+plotStringX+">>hCut",weight*cutCut,"colz");
  hCut->GetZaxis()->SetRangeUser(1,1e10);
  gPad->SetLogz();
  gPad->Update();

  TStopwatch ts;
  can->cd(1);
  TH2D hLevel2("hLevel2","Level 2;"+plotStringX+";"+plotStringY,
	    b[0],b[1],b[2],b[3],b[4],b[5]);
  treeRaw->Draw(plotStringY+":"+plotStringX+">>hLevel2",weight*cutRaw,"colz");
  hLevel2->GetZaxis()->SetRangeUser(1,1e10);
  gPad->SetLogz();
  ts.Print();

  can->cd(4);
  TH2D hRatio("hRatio","Ratio (Cut / Level 2);"+plotStringX+";"+plotStringY,
	    b[0],b[1],b[2],b[3],b[4],b[5]);
  hRatio.Divide(&hCut,&hLevel2);
  //  hRatio.GetZaxis()->SetRangeUser(0,1);
  hRatio.Draw("colz");
}


