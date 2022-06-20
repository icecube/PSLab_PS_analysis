#include "plotscripts/public/EffectiveArea.h"

#include "TTreeFormula.h"
#include "TCut.h"
#include "TMath.h"
#include "TDirectory.h"

#include "rootExt/public/FunctionsRoot.h"
#include "rootExt/public/TStringify.h"



void EffectiveArea::SetEnergyRange(int nBins, double logEMin, double logEMax)
{
  nBins_ = nBins;
  logEMin_ = logEMin;
  logEMax_ = logEMax;
}

TH1D* EffectiveArea::Calculate(TTree *tree, double zenMinDeg, double zenMaxDeg, 
			       TCut cut, TString hname)
{
  // Check that our tree has the required aliases defined:

  if (!VerifyExpression(tree,"mcPrimary_Zenith_rad") ||
      !VerifyExpression(tree,"mcPrimary_Energy_GeV") ||
      !VerifyExpression(tree,"mcOneWeight") ||
      !VerifyExpression(tree,"mcTotalGeneratedEvents") ) 
  {
    return NULL;
  }

  if (gDirectory->FindObject(hname)) { delete gDirectory->FindObject(hname);}
  TH1D *hist = new TH1D(hname,
			"Solid-Angle-Averaged Effective Area;"
			"log_{10} E_{#nu} / GeV;Area [m^{2}]",
			nBins_, logEMin_, logEMax_);

  double zenMinRad = zenMinDeg*TMath::DegToRad();
  double zenMaxRad = zenMaxDeg*TMath::DegToRad();

  double solidAngle = 2 * TMath::Pi() * (cos(zenMinRad)-cos(zenMaxRad));

  tree->SetAlias("SolidAngle",TStringify(solidAngle));

  TCut angleMinCut = TCut("mcPrimary_Zenith_rad>"+TStringify(zenMinRad) );
  TCut angleMaxCut = TCut("mcPrimary_Zenith_rad<"+TStringify(zenMaxRad) );

  double EBinsPerDecade = nBins_/(logEMax_-logEMin_);
  tree->SetAlias("EBinsPerDecade", TStringify(EBinsPerDecade) );
  tree->SetAlias("mcLogEBin",
		 "int(log10(mcPrimary_Energy_GeV)*EBinsPerDecade)");
  tree->SetAlias("mcEMin", "pow(10., mcLogEBin/EBinsPerDecade)");
  tree->SetAlias("mcEMax", "pow(10., (1+mcLogEBin)/EBinsPerDecade)");

  tree->Draw("log10(mcPrimary_Energy_GeV)>>"+hname, 
	     cut * angleMinCut * angleMaxCut * 
	     "1e-4*(mcOneWeight/mcTotalGeneratedEvents)/"
	     "(SolidAngle*(mcEMax-mcEMin))","goff");

  return hist;
}
