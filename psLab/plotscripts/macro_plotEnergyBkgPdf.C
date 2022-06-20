// This macro draws the energy background PDF vs cos(zenith) for a particular
// dataset (here loaded as ark59). 

// A figure like this was used in the IC40 flare paper

{

  int nbinsEn = 44;
  int nbinsZn = 250;

  double zn, az, val, bProb;
  
  I3Ark arkTemp = ark59;
 
  arkTemp.decBkgProb.FixToBase();
    
  //TH2D * hePdf = new TH2D("","",nbinsEn,2,8,nbinsZn,-1,1); // This matches the space bkg pdf
  TH2D * hePdf = new TH2D("","",nbinsZn,-1,1,nbinsEn,2,8);

     
  EquatorialDeg coo;
  
  I3Event ev;
  I3EventParameters pa;
  
  for (int iEn=1;iEn<=nbinsEn;iEn++) {
    //pa.energyValue = hePdf->GetXaxis()->GetBinCenter(iEn);
    pa.energyValue = hePdf->GetYaxis()->GetBinCenter(iEn);
      for (int iDec=1;iDec<=nbinsZn;iDec++) {
//      pa.recoZenithDeg = hePdf->GetYaxis()->GetBinCenter(iDec);
//      pa.recoZenithDeg = 57.3*acos(hePdf->GetYaxis()->GetBinCenter(iDec));

        pa.recoZenithDeg = 57.3*acos(hePdf->GetXaxis()->GetBinCenter(iDec));
//      pa.recoZenithDeg = 90.+57.3*asin(hePdf->GetYaxis()->GetBinCenter(iDec));
      //cout << pa <<  " " << flush;
      
      ev.SetParams(pa);
      
      val = arkTemp.eProb->GetEnergyProbBkg(ev);
        
      //hePdf->SetBinContent(iEn,iDec,val); // this matches local coord bkg in zenith
      hePdf->SetBinContent(iDec,iEn,val);
    }
  }  

hePdf->GetYaxis()->SetTitle("log_{10} (Energy Proxy)");
//hePdf->GetYaxis()->SetTitle("Zenith (#circ)");
hePdf->GetXaxis()->SetTitle("cos(Zenith)");

hePdf->GetZaxis()->SetTitle("B^{Energy}_{i}");
hePdf->GetZaxis()->SetTitleOffset(1);

TCanvas * c = new TCanvas("c","c");
c->SetRightMargin(0.125);
hePdf->GetZaxis()->SetTitleOffset(0.5);
hePdf->Draw("colz");
    
}
