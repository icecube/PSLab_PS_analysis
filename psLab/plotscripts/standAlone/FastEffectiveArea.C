
// The histogram returned by this function gives the *solid-angle-averaged* 
// effective area, in units of m^2, as a function of energy.

// It is calculated as an average over a zenith-angle range, but it 
// represents the effective area for a point-source flux at a single zenith 
// within this range.




#ifndef ROOTEXT_FUNCTIONSROOT_H_
bool VerifyExpression(TTree *tree, const char* expression) {
  TTreeFormula formula("","1", tree);  // initialize with expression="1", (OK)
  if (formula.Compile(expression)) {   // now test expresssion
    return false;
  }
  return true;
}
#endif

#ifndef ROOTEXT_TSTRINGIFY_H_
TString TStringify(double x) {
  char temp[100];
  // formerly, this was %lg, but that didn't keep enough digits for double
  // precision.  This should work better: up to 16 digits if necessary.
  sprintf(temp,"%.16g",x);
  return TString(temp);
}
#endif



TH1D* FastEffectiveArea(TTree *tree, TCut cut="", 
			double zenMinDeg=90, double zenMaxDeg=180, 
			int EBins=40, double LogEMin=1, double LogEMax=9,
			TString hname="hFastEffArea") {

  // Check that our tree has the required aliases defined:

  if (!VerifyExpression(tree,"mcPrimary_Zenith_rad") ||
      !VerifyExpression(tree,"mcPrimary_Energy_GeV") ||
      !VerifyExpression(tree,"mcOneWeight") ) //||
      //!VerifyExpression(tree,"mcTotalGeneratedEvents") ) 
  {
    cout << "\nMake sure tree has variables or aliases for:\n";
    cout << "  mcPrimary_Zenith_rad\n";
    cout << "  mcPrimary_Energy_GeV\n";
    cout << "  mcOneWeight\n";
    //cout << "  mcTotalGeneratedEvents"
      " (no. of mc events generated per file * no. of files)\n\n";
    return NULL;
  }

  if (gDirectory->FindObject(hname)) { delete gDirectory->FindObject(hname);}
  TH1D *hFastEffArea = new TH1D(hname,"Effective Area",
				EBins, LogEMin, LogEMax);

  cout << zenMinDeg << " " << zenMaxDeg << endl;

  double zenMinRad = zenMinDeg*TMath::DegToRad();
  double zenMaxRad = zenMaxDeg*TMath::DegToRad();

  double solidAngle = 2 * TMath::Pi() * (cos(zenMinRad)-cos(zenMaxRad));

  tree->SetAlias("SolidAngle",TStringify(solidAngle));
  TCut solidAngCut = "1/"+TStringify(solidAngle);

  TCut angleMinCut = "mcPrimary_Zenith_rad>"+TStringify(zenMinRad);
  TCut angleMaxCut = "mcPrimary_Zenith_rad<"+TStringify(zenMaxRad);
  TCut normCut = "1./mcTotalGeneratedEvents";//+TStringify(mcTotalGeneratedEvents);

  double EBinsPerDecade = EBins/(LogEMax-LogEMin);
  //tree->SetAlias("EBinsPerDecade", TStringify(EBinsPerDecade) );
  tree->SetAlias("mcLogEBin",
		 "int(log10(mcPrimary_Energy_GeV)*"+TStringify(EBinsPerDecade)+")");
  tree->SetAlias("mcEMin", "pow(10., mcLogEBin/"+TStringify(EBinsPerDecade)+")");
  tree->SetAlias("mcEMax", "pow(10., (1+mcLogEBin)/"+TStringify(EBinsPerDecade)+")");		 

//  tree->SetAlias("mcEMin", "pow(10., mcLogEBin/EBinsPerDecade)");
//  tree->SetAlias("mcEMax", "pow(10., (1+mcLogEBin)/EBinsPerDecade)");

  //tree->SetAlias("mcTGE", TStringify(mcTotalGeneratedEvents * (LogEMax/LogEMax)) );

  // factor of 1e-4 below converts from cm2 (intrinsinc units of OneWeight)
  // to m2, which is traditional effective area unit for IceCube plots

//  tree->Draw("log10(mcPrimary_Energy_GeV)>>"+hname, 
//	     cut * angleMinCut * angleMaxCut * 
//	     "1e-4*(mcOneWeight/mcTotalGeneratedEvents)/"
//	     "(SolidAngle*(mcEMax-mcEMin))","goff");

  tree->Draw("log10(mcPrimary_Energy_GeV)>>"+hname, 
	     cut * angleMinCut * angleMaxCut * solidAngCut * normCut *
	     "1e-4*(mcOneWeight)/(mcEMax-mcEMin)","goff");

  return hFastEffArea;
}
