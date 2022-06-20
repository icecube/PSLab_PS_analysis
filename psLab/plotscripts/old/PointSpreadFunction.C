
#include "TTree.h"
#include "TCut.h"
#include "TString.h"
#include "TH1D.h"
#include <iostream>
#include "TDirectory.h"
#include "MakeNewName.h"

TH1D* PointSpreadFunction(TTree *tree, TCut cut, 
			  TString recoZenRad, TString recoAziRad,
			  double maxAngDeg=10, int nBins=100, 
			  TString hname = "hPSF") {

  // Check that our tree has the required aliases defined:

  if (!tree->GetBranch("mcZr") && !tree->GetAlias("mcZr") ) {
    cout << "Error.  Must define \"mcZr\" alias for mc true zenith (Rad)\n";
    return NULL;
  }
  if (!tree->GetBranch("mcAr") && !tree->GetAlias("mcAr") ) {
    cout << "Error.  Must define \"mcAr\" alias for mc true azimuth (Rad)\n";
    return NULL;
  }

  if (gDirectory->FindObject(hname)) { delete gDirectory->FindObject(hname);}
  TH1D *hPSF = new TH1D(hname,"; #Psi [#circ]",
			nBins, 0, maxAngDeg);

  TString expression="SpaceAngleRad(mcZr,mcAr,"+recoZenRad+","+recoAziRad+")";

  tree->Draw("TMath::RadToDeg()*"+expression+">>"+hname, cut, "goff");

  return hPSF;
}



TH1* HistIntegrate(TH1 *h) {

  TH1* hnew = dynamic_cast<TH1*> (h->Clone(MakeNewName(h->GetName())));

  double sum = 0.;

  // note: 0 is the underflow bin, and nBins+1 is the overflow bin
  for (int i=0; i <= h->GetNbinsX()+1; ++i) {
    sum += h->GetBinContent(i);
    hnew->SetBinContent(i, sum);
  }

  return hnew;
}


TH1D* HistNormalize(TH1D *h) {
  h->Scale(1./h->GetSum());
  return h;
}

