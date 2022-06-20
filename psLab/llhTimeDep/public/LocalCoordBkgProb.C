#include "llhTimeDep/public/LocalCoordBkgProb.h"


void LocalCoordBkgProb::FillLCBkgHisto(const vector<I3Event>& events) {
  dataAZ.Reset();
  if(!useCosZ) { 
    dataAZ.SetBins(nbAz,0.,360.,nbZn,0.,180.);
  } else {
    dataAZ.SetBins(nbAz,0.,360.,nbZn,-1,1.);
  }
  
  dataZT.Reset();
  if (nBinsT) { dataZT.SetBins(nBinsT,tMin,tMax); }

  for (unsigned int i=0; i<events.size(); i++) {
    if(!useCosZ) {
      dataAZ.Fill(events[i].GetParams().recoAzimuthDeg,events[i].GetParams().recoZenithDeg);
    } else { 
      dataAZ.Fill(events[i].GetParams().recoAzimuthDeg,
           -1.*cos(events[i].GetParams().recoZenithDeg/57.2957795) ); 
    }
    if (events[i].GetParams().recoZenithDeg < 90. && nBinsT) {
      dataZT.Fill(events[i].GetTime().GetTime());
    }
  }

  pdfAZ.Reset();
  if(!useCosZ) {
    pdfAZ.SetBins(nbAz,0.,360.,nbZn,0.,180.);
  } else { 
    pdfAZ.SetBins(nbAz,0.,360.,nbZn,-1.,1.);
  }
  
  double n[nbZn+1];
  for (int i=1;i<=nbZn;i++) {
    n[i] = dataAZ.ProjectionX("px1",i,i,"")->GetSum();
  }
  
  for (int i=0;i<=nbAz;i++) {
    for (int j=0;j<=nbZn;j++) {
      if (n[j]) {
        pdfAZ.SetBinContent( i,j, (dataAZ.GetBinContent(i,j)*nbAz/n[j]) );
      } else {
        pdfAZ.SetBinContent(i,j, 1.);
      }
    }
  }

  if (nBinsT) {
    pdfZT.Reset();
    pdfZT.SetBins(nBinsT,tMin,tMax);

    double sum = dataZT.GetEntries();
                  
    for (int i=1;i<=dataZT.GetNbinsX();i++) {
      pdfZT.SetBinContent(i, 1/(dataZT.GetBinContent(i)*nBinsT/sum) );
    }
  } 
}

void LocalCoordBkgProb::FillLCBkgHisto_folded(const vector<I3Event>& events,double minDegZenith,double maxDegZenith) {

    dataAZ.Reset();
    if(!useCosZ) { 
        dataAZ.SetBins(nbAz,0.,60.,nbZn,minDegZenith,maxDegZenith);
    } else {
        dataAZ.SetBins(nbAz,0.,60.,nbZn,-cos(minDegZenith*TMath::DegToRad()),-cos(maxDegZenith*TMath::DegToRad()));
    }

    dataZT.Reset();
    if (nBinsT) { dataZT.SetBins(nBinsT,tMin,tMax); }

    for (unsigned int i=0; i<events.size(); i++) {
        double folded_az = events[i].GetParams().recoAzimuthDeg - (int)(events[i].GetParams().recoAzimuthDeg/60.0)*60;
        if(!useCosZ) {
            dataAZ.Fill(folded_az,events[i].GetParams().recoZenithDeg);
        } else { 
            dataAZ.Fill(folded_az,-1.*cos(events[i].GetParams().recoZenithDeg/57.2957795) ); 
        }
        if (events[i].GetParams().recoZenithDeg < 90. && nBinsT) {
            dataZT.Fill(events[i].GetTime().GetTime());
        }
    }

    pdfAZ.Reset();
    if(!useCosZ) {
        pdfAZ.SetBins(nbAz,0.,60.,nbZn,minDegZenith,maxDegZenith);
    } else { 
        pdfAZ.SetBins(nbAz,0.,60.,nbZn,-cos(minDegZenith*TMath::DegToRad()),-cos(maxDegZenith*TMath::DegToRad()));
    }
  
    double n[nbZn+1];
    for (int i=1;i<=nbZn;i++) {
        n[i] = dataAZ.ProjectionX("px1",i,i,"")->GetSum();
    }
  
    for (int i=0;i<=nbAz;i++) {
        for (int j=0;j<=nbZn;j++) {
            if (n[j]) {
                pdfAZ.SetBinContent( i,j, (dataAZ.GetBinContent(i,j)*nbAz/n[j]) );
            } else {
                pdfAZ.SetBinContent(i,j, 1.);
            }
        }
    }
    
    if (nbZn==1) {
        TH1D *pdfAZ_1D_tmp=pdfAZ.ProjectionX("pdfAZ_1D",1,1,"");
        pdfAZ_1D=new TGraph(pdfAZ_1D_tmp->GetNbinsX()+4); 
        double border_value=0.5*(pdfAZ_1D_tmp->GetBinContent(pdfAZ_1D_tmp->GetNbinsX())+pdfAZ_1D_tmp->GetBinContent(1));
        pdfAZ_1D->SetPoint(0,pdfAZ_1D_tmp->GetBinCenter(1)-2*pdfAZ_1D_tmp->GetBinWidth(1),pdfAZ_1D_tmp->GetBinContent(pdfAZ_1D_tmp->GetNbinsX()-1));        
        pdfAZ_1D->SetPoint(1,pdfAZ_1D_tmp->GetBinCenter(1)-pdfAZ_1D_tmp->GetBinWidth(1),border_value);
        pdfAZ_1D->SetPoint(2,pdfAZ_1D_tmp->GetBinCenter(1),border_value);
        for (int i=1;i<pdfAZ_1D_tmp->GetNbinsX()-1;i++) {
            pdfAZ_1D->SetPoint(i+2,pdfAZ_1D_tmp->GetBinCenter(i+1),pdfAZ_1D_tmp->GetBinContent(i+1));
        }
        
        pdfAZ_1D->SetPoint(pdfAZ_1D_tmp->GetNbinsX()+1,pdfAZ_1D_tmp->GetBinCenter(pdfAZ_1D_tmp->GetNbinsX()),border_value);
        pdfAZ_1D->SetPoint(pdfAZ_1D_tmp->GetNbinsX()+2,pdfAZ_1D_tmp->GetBinCenter(pdfAZ_1D_tmp->GetNbinsX())+pdfAZ_1D_tmp->GetBinWidth(1),border_value);
        pdfAZ_1D->SetPoint(pdfAZ_1D_tmp->GetNbinsX()+3,pdfAZ_1D_tmp->GetBinCenter(pdfAZ_1D_tmp->GetNbinsX())+2*pdfAZ_1D_tmp->GetBinWidth(1),pdfAZ_1D_tmp->GetBinContent(2));
    }
    if (nBinsT) {
        pdfZT.Reset();
        pdfZT.SetBins(nBinsT,tMin,tMax);

        double sum = dataZT.GetEntries();
                  
        for (int i=1;i<=dataZT.GetNbinsX();i++) {
            pdfZT.SetBinContent(i, 1/(dataZT.GetBinContent(i)*nBinsT/sum) );
        }
    } 
}

void LocalCoordBkgProb::SetTimeParams(double min, double max, int nbins){
  tMin = min;
  tMax = max;
  nBinsT = nbins;
}


double LocalCoordBkgProb::BackgroundLCProbA_mcozZ(double azimuth,double zenith) {
    double p1;
    if(useCosZ) {
        p1 = pdfAZ.Interpolate(azimuth,zenith);
    }
   return p1;
}
    
double LocalCoordBkgProb::BackgroundLCProbA_mcozZ_folded(double azimuth,double zenith) {
    double p1;
    if(useCosZ) {
        double azimuth_folded  = azimuth - 60.*(int)(azimuth/60.0);
        if (nbZn==1) p1 = pdfAZ_1D->Eval(azimuth_folded,0, "S");
        else p1 = pdfAZ.Interpolate(azimuth_folded,zenith);
    }
   return p1;
}
    
double LocalCoordBkgProb::BackgroundLCProb(const I3Event& ev) {
  double mjd = ev.GetMJD();
  double azimuth = ev.GetParams().recoAzimuthDeg;
  double zenith = ev.GetParams().recoZenithDeg;
  double p1=1.;
  
  if(!useCosZ) {
    p1 = pdfAZ.Interpolate(azimuth,zenith);
  } else {
    p1 = pdfAZ.Interpolate(azimuth,-1*cos(zenith/57.2957795));
  }

  double p2=1.;
  int t;

  if (nBinsT) {
    t = pdfZT.GetXaxis()->FindBin(mjd); // if we want atmospheric variation
    if (zenith<90) {
      p2 = pdfZT.GetBinContent(t);
    } else {
      p2 = 1.;
    }
  }

  return p1 * p2;
}

double LocalCoordBkgProb::BackgroundLCProb_folded(const I3Event& ev) {
    double mjd = ev.GetMJD();
    double azimuth_folded = ev.GetParams().recoAzimuthDeg - (int)(ev.GetParams().recoAzimuthDeg/60.0)*60;
    double zenith = ev.GetParams().recoZenithDeg;
    double p1=1.;
   
    if (nbZn==1) {
        pdfAZ_1D->Eval(azimuth_folded,0, "S");
        /*
        float x;
        float sum=0;
        for (int i=0; i<100000; i++){
            x=60*i/100000.;
            sum=sum+pdfAZ_1D->Eval(x,0, "S");
        }
        cout << "norm check " <<  sum << " " << sum/100000. << endl;
        exit(0);
        */
    }
    else {
        if(!useCosZ) p1 = pdfAZ.Interpolate(azimuth_folded,zenith);
        else p1 = pdfAZ.Interpolate(azimuth_folded,-1*cos(zenith/57.2957795));
    }

    double p2=1.;
    int t;

    if (nBinsT) {
        t = pdfZT.GetXaxis()->FindBin(mjd); // if we want atmospheric variation
        if (zenith<90) {
            p2 = pdfZT.GetBinContent(t);
        } else {
        p2 = 1.;
        }
    }

    return p1 * p2;
}
