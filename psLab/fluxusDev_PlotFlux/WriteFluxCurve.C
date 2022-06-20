
{
  gROOT->ProcessLine(".L FluxTools.C");
  gROOT->ProcessLine(".L PlotFlux.C+");
  gStyle->SetPadTopMargin(0.05);


  // Assume these are already set:
  //    srcTree = LoadTree_mayIC22_nugen651_extuple();
  //    sourceZenWidthDeg = 4.;
  //    evLoader.SetSourceTree(srcTree);
  //    evLoader.SetSourceZenWidthDeg(sourceZenWidthDeg);
  //
  //    llhPtr (configured w/ or w/o energy, etc.)
  //    disco (loops, disc or upper limit, est. median, etc.)

  double srcDecDeg = 6.00;
  EquatorialDeg srcLocation(90.,srcDecDeg);
  vector<I3Event> sourceEvents;
  evLoader.LoadSourceEvents(sourceEvents, srcLocation);

  double binsPerDecade = 2;

  double spectralIndex = -2.0;

  double baseThinning = 10./20.; // src livetime will be scaled by this factor

  psData.SetBaseThinningProb(baseThinning);

  char *filename = "results/fluxDecades_ic40_10years_Em2_disc_E_dec6.txt";

  // Need to make General
  //  TString pullRangeString = "g32SpaceAngleDeg/pf32SigmaDeg<2";
  TString pullRangeString = "mDelAng/mpfSigmaDeg<2";
  // COULD SWITCH THIS TO MPFSIGMA CORRECTED... THIS IS JUST ESTIMATOR ANYWAY

  int thresholdEvents = 4;


  disco.SetLogLikelihood(*llhPtr);

  vector<double> eMinVect;
  vector<double> eMaxVect;
  vector<double> indexVect;
  vector<double> fluxVect;
  vector<double> nSigVect;

  vector<TString> fluxStringVect;

  double logEMin = 2.;
  double logEMax = 9.;

  TH1D hEnergyRange("hEnergyRange","hEnergyRange",
		    (logEMax-logEMin)*binsPerDecade,logEMin,logEMax);
  TString zenRangeString = "abs(mcPrimary_Zenith_rad*TMath::RadToDeg()-";
  zenRangeString += TStringify(90+srcDecDeg);
  zenRangeString += ")<"+TStringify(sourceZenWidthDeg);
  TCut zenRangeCut = zenRangeString;
  TCut pullRangeCut = pullRangeString;

  TCanvas *canER = new TCanvas("canER");
  srcTree->Draw("log10(mcPrimary_Energy_GeV)>>hEnergyRange",
		evLoader.GetCuts()*zenRangeCut*pullRangeCut);
  canER->Update();

  double logEinc = 1./binsPerDecade;

  // DECADES
  if (1) {
    for (double logE=logEMin; logE<logEMax; logE += logEinc) {
      eMinVect.push_back(pow(10,logE)) ; 
      eMaxVect.push_back(pow(10,logE+logEinc));
      indexVect.push_back(spectralIndex);
      fluxVect.push_back(0);
      nSigVect.push_back(0);
      fluxStringVect.push_back("");
    }
  }

  // CUTOFFS
  if (0) {
    for (double logE=logEMin+logEinc; logE<=logEMax; logE += logEinc) {
      eMinVect.push_back(pow(10,logEMin)) ; 
      eMaxVect.push_back(pow(10,logE));
      indexVect.push_back(spectralIndex);
      fluxVect.push_back(0);
      nSigVect.push_back(0);
      fluxStringVect.push_back("");
    }
  }
  
  TCanvas *can = new TCanvas("can","can",20,20,800,800);
  can->Divide(1,2,0.001,0.001);

  FluxPlotManager fpm;
  fpm.yPower = fabs(spectralIndex);

  can->cd(1);
  TH1D hFrameEvents;
  hFrameEvents.SetBins(1,logEMin,logEMax);
  hFrameEvents.SetTitle(";log_{10} E/GeV;Events");
  hFrameEvents.SetMaximum(30);
  hFrameEvents.Draw();

  can->cd(2);
  TH1D hFrameFlux;
  hFrameFlux.SetBins(1,logEMin,logEMax);
  hFrameFlux.SetTitle(";log_{10} E/GeV;E^{#gamma} TeV^{-1} cm^{-2} s^{-1}");
  hFrameFlux.Draw();
  hFrameFlux.SetMinimum(1e-12);
  gPad->SetLogy(1);


  for (int i=0; i<eMinVect.size(); ++i) {

    // x is in GeV
    TString fluxString = "pow(x,";
    fluxString += TStringify(indexVect[i]);
    fluxString += ")*(x > ";
    fluxString += TStringify(eMinVect[i]);
    fluxString += " && x < ";
    fluxString += TStringify(eMaxVect[i]);
    fluxString += ")";

    FormulaFlux thisRangeFlux(fluxString);
    I3PointGenerator thisSignal(sourceEvents, thisRangeFlux, srcLocation,
				livetime * baseThinning);

    psData.SetSource(thisSignal);

    double meanSrcEv_ForDetection;
    double meanFlux_ForDetection;

    // are there any signal events at all in this range?
    TString eString = "mcPrimary_Energy_GeV > " + TStringify(eMinVect[i]);
    eString += " && mcPrimary_Energy_GeV < " + TStringify(eMaxVect[i]);
    TCut eCut = eString;
    int evPass = 
      srcTree->Draw("1", evLoader.GetCuts()*zenRangeCut*pullRangeCut*eCut);
    cout << "\n" << evPass << " good events in range:  " << eString << endl;

    if ( evPass > thresholdEvents ) {

      disco.AnalyzeDiscoveryPotential(&psData,srcLocation);
      meanSrcEv_ForDetection = disco.MeanSrcEv_ForDetection_;
      meanFlux_ForDetection = 
	thisSignal.GetFluxScaleForNev(meanSrcEv_ForDetection);
    }
    else {
      // No signal events, so can't estimate discovery potential!
      meanSrcEv_ForDetection = 0.;
      meanFlux_ForDetection = 0.;
    }


    fluxVect[i] = meanFlux_ForDetection;
    nSigVect[i] = meanSrcEv_ForDetection;
    fluxStringVect[i] = TStringify(meanFlux_ForDetection)+"*"+fluxString;


    can->cd(1);
    TLine *eventLine = new TLine(log10(eMinVect[i]),meanSrcEv_ForDetection,
				 log10(eMaxVect[i]),meanSrcEv_ForDetection);
    eventLine->Draw();

    can->cd(2);
    TGraph *g = GraphFlux(fluxStringVect[i], 0, eMinVect[i], eMaxVect[i]);
    fpm.AddGraph(*g);
    fpm.Plot();
    can->Update();
  }

  FILE *fp = fopen(filename,"w");

  for (int i=0; i<eMinVect.size(); ++i) {
    fprintf(fp, "%lg %lg %lg %lg %lg\n",
	  eMinVect[i], eMaxVect[i], indexVect[i],
	  fluxVect[i], nSigVect[i]);
    fprintf(fp, "%s\n",dynamic_cast<char*>(fluxStringVect[i]));
  }
  
  fclose(fp);

}
