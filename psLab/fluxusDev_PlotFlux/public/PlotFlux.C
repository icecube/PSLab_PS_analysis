#include "fluxusDev_PlotFlux/public/PlotFlux.h"

#include "TCanvas.h"
#include "TPad.h"

#include "rootExt/public/TStringify.h"



void FluxPlotManager::Plot() {
  hFrame.SetBins(1,logEMin,logEMax);
    
  TString titlePower = "E^{"+TStringify(yPower)+"}";
  if (yPower==1) { titlePower = "E"; }
  if (yPower==0) { titlePower = ""; }

  TString titleUnits = "("+TStringify(yOutputUnits/1e9)+" GeV)";
  if (yOutputUnits==1e9) { titleUnits = "GeV"; }
  if (yOutputUnits==1e12) { titleUnits = "TeV"; }

  if (yPower==1) { titleUnits = ""; }
  else if (yPower == 2) { }  // no change
  else {
    titleUnits += "^{"+TStringify(yPower-1)+"}";
  }

  titleUnits += " cm^{-2} s^{-1}";

  hFrame.SetTitle(";log_{10} E/GeV;"+titlePower+" d#Phi/dE  [" +
		  titleUnits + "]");

  plotGraphVect = graphVect;

  double yLowest = 0;
  double yHighest = 0;
  double x, y;

  double yScaleFactor = yOutputUnits/1e9;  // GeV is original reference unit

  for (unsigned int i=0; i<plotGraphVect.size(); ++i) {
    for (int j=0; j<plotGraphVect[i].GetN(); ++j) {

      plotGraphVect[i].GetPoint(j,x,y);
      y *= yScaleFactor * pow( pow(10,x)/yScaleFactor, yPower);
      plotGraphVect[i].SetPoint(j,x,y);

      if (y>yHighest) { yHighest = y; }
      if (y>0) {
	if (yLowest==0 || y<yLowest) { yLowest = y; }
      }
    }
  }

  if (yLowest==0) {
    cout << "Error: flux was <=0 everywhere.\n";
    return;
  }

  if (yMin > 0 && yMax > 0) {
    hFrame.SetMinimum(yMin);
    hFrame.SetMaximum(yMax);
  } else {
    double log10Ratio = log10(yHighest/yLowest);
    if (log10Ratio<1) { log10Ratio = 1.; }
    hFrame.SetMinimum(yLowest/pow(10,0.1*log10Ratio));
    hFrame.SetMaximum(yHighest*pow(10,0.1*log10Ratio));
  }

  hFrame.Draw();
  gPad->SetLogy(1);
  for (unsigned int i=0; i<plotGraphVect.size(); ++i) {
    plotGraphVect[i].Draw("L");
  }
}



TGraph* GraphFlux(FluxBase &flux, double yPower, 
		  double xmin, double xmax, int nPoints)
{
  double logEMin = log10(xmin);
  double logEMax = log10(xmax);

  double log_inc = (logEMax-logEMin)/nPoints;

  TGraph *graphFn = new TGraph(nPoints);
  for (int i=0; i<nPoints; ++i) {
    double logE = logEMin+(i+0.5)*log_inc;
    double x = pow(10,logE);
    double y = flux.GetFlux(x)*pow(x,yPower);
    graphFn->SetPoint(i,logE,y);
  }
    
  return graphFn;
}



TGraph* PlotFlux(FluxBase &flux, double yPower, 
		 double xmin, double xmax, 
		 double ymin, double ymax)
{
  TCanvas *fluxCanvas = new TCanvas("fluxCanvas","fluxCanvas"); //,20,20,800,600);

  int nPoints = 250;
  double logEMin = log10(xmin);
  double logEMax = log10(xmax);

  double log_inc = (logEMax-logEMin)/nPoints;

  TGraph *graphFn = new TGraph(nPoints);
  for (int i=0; i<nPoints; ++i) {
    double logE = logEMin+(i+0.5)*log_inc;
    double x = pow(10,logE);
    double y = flux.GetFlux(x)*pow(x,yPower);
    graphFn->SetPoint(i,logE,y);

    if (y>ymax) { ymax = y; }
    if (y>0) {
      if (ymin==0 || y<ymin) { ymin = y; }
    }
  }

  if (ymin==0) {
    cout << "Error: flux was <=0 everywhere.\n";
    return NULL;
  }

  TH1D *hFrame = new TH1D("","",1,log10(xmin),log10(xmax));
  double log10Ratio = log10(ymax/ymin);
  if (log10Ratio<1) { log10Ratio = 1.; }
  hFrame->SetMinimum(ymin/pow(10,0.1*log10Ratio));
  hFrame->SetMaximum(ymax*pow(10,0.1*log10Ratio));
  hFrame->Draw();
  gPad->SetLogy(1);

  fluxCanvas->cd();
  graphFn->Draw("L");

  // Setting ymin by default
    
  return graphFn;
}
