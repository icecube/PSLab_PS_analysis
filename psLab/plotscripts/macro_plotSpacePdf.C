// This macro plots a 2D histogram displaying the local coordinate
// correction function for a loaded dataset (here ark59).
// This is important in the untriggered flare analysis, and this figure
// was included in the IC40 flare paper.

{

  int nbinsAz = 360;
  int nbinsZn = 180;

  double zn, az, val, bProb;
 
//  psData.SetSource(*mySignalPtr);
//  psData.SetBaseEvents(baseEvents,true,16,90,true);
//  psData.UseRealData();

  ark59.decBkgProb.FixToBase();
  //ark59.lcBkgProb.SetDecProb(ark59.decBkgProb);
    
  TH2D * hsPdf = new TH2D("","",nbinsAz,0,360,nbinsZn,-1,1);
  TH2D * hsPdf1 = new TH2D("","",nbinsAz,0,360,nbinsZn,-1,1);
  TH2D * hsPdf2 = new TH2D("","",nbinsAz,0,360,nbinsZn,-1,1);
     
  EquatorialDeg coo;
  
//  Coord * c = *coo;
  I3EventParameters par;
  I3Event ev;
  
  for (int iAz=1;iAz<=nbinsAz;iAz++) {
    az = hsPdf->GetXaxis()->GetBinCenter(iAz);
      for (int iDec=1;iDec<=nbinsZn;iDec++) {
//      zn = hsPdf->GetYaxis()->GetBinCenter(iDec);
//      zn = 90.+57.3*asin(hsPdf2->GetYaxis()->GetBinCenter(iDec));
      zn = 57.3*acos(hsPdf2->GetYaxis()->GetBinCenter(iDec));
   
      coo.SetCoords(az,zn-90.);
         
      par.recoAzimuthDeg = az;
      par.recoZenithDeg  = zn;
      
      ev.SetParams(par);

      bProb = ark59.decBkgProb.GetBkgProbDensity(coo);            
      val = ark59.lcBkgProb.BackgroundLCProb(ev);
           
      //    psData.GetEnergyProb()->GetEnergyProbBkg(ev);
        
      hsPdf->SetBinContent(iAz,iDec,val);
      hsPdf1->SetBinContent(iAz,iDec,bProb*1e5);
      hsPdf2->SetBinContent(iAz,iDec,val*bProb*1e5);
    }
  }  

hsPdf->GetXaxis()->SetTitle("Azimuth (#circ)");
//hsPdf->GetYaxis()->SetTitle("Zenith (#circ)");
hsPdf1->GetXaxis()->SetTitle("Azimuth (#circ)");
//hsPdf1->GetYaxis()->SetTitle("Zenith (#circ)");
hsPdf2->GetXaxis()->SetTitle("Azimuth (#circ)");
//hsPdf2->GetYaxis()->SetTitle("Zenith (#circ)");
hsPdf->GetYaxis()->SetTitle("cos(Zenith)");
hsPdf1->GetYaxis()->SetTitle("cos(Zenith)");
hsPdf2->GetYaxis()->SetTitle("cos(Zenith)");


hsPdf2->GetZaxis()->SetTitle("B^{space}_{i} * 10^{5}");
hsPdf1->GetZaxis()->SetTitle("10^{5}/#Omega");
hsPdf->GetZaxis()->SetTitle("LocalCoordTransform");
//hsPdf2->GetZaxis()->CenterTitle();
hsPdf->GetZaxis()->SetTitleOffset(0.8);
hsPdf1->GetZaxis()->SetTitleOffset(0.7);
hsPdf2->GetZaxis()->SetTitleOffset(0.8);

TCanvas * c = new TCanvas("c","c");
  c->SetRightMargin(0.125);
  hsPdf->Draw("colz");
TCanvas * c1 = new TCanvas("c1","c1");
  c1->SetRightMargin(0.125);
  hsPdf1->Draw("colz");
TCanvas * c2 = new TCanvas("c2","c2");
  c2->SetRightMargin(0.13);
  hsPdf2->Draw("colz");
    
}
