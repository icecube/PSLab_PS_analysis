
{
  gROOT->ProcessLine(".L FluxTools.C");
  gROOT->ProcessLine(".L PlotFlux.C+");

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.10);


  // Assume these are already set:
  //    srcTree = LoadTree_mayIC22_nugen651_extuple();
  //    sourceZenWidthDeg = 4.;
  //    evLoader.SetSourceTree(srcTree);
  //    evLoader.SetSourceZenWidthDeg(sourceZenWidthDeg);

  double srcDecDeg = 6.00;
  EquatorialDeg srcLocation(90.,srcDecDeg);

  sourceZenWidthDeg = 1.;
  evLoader.SetSourceZenWidthDeg(sourceZenWidthDeg);

  vector<I3Event> sourceEvents;
  evLoader.LoadSourceEvents(sourceEvents, srcLocation);

  double binsPerDecade = 2;

  double spectralIndex = -2.0;

  bool discovery = false;

  bool zenUseEnergy = true;

  double baseThinning = 1.00;  // src livetime will be scaled by this factor
  // DONT SET FOR NOW, RUNNING WITH JONS VERSION
  //  double baseThinning = 0.003;  // src livetime will be scaled by this factor
  //  psData.SetBaseThinningProb(baseThinning);

  char *filename = "results/fluxDecades_ic40_full_Em2_sens_E_dec06.txt";

  // Need to make General
  //  TString pullRangeString = "g32SpaceAngleDeg/pf32SigmaDeg<2";
  TString pullRangeString = "mDelAng/mpfSigmaDeg<2";
  // COULD SWITCH THIS TO MPFSIGMA CORRECTED... THIS IS JUST ESTIMATOR ANYWAY


  vector<double> eMinVect;
  vector<double> eMaxVect;
  vector<double> indexVect;
  vector<double> fluxVect;
  vector<double> nSigVect;

  vector<TString> fluxStringVect;

  double logEMin = 2.5;
  double logEMax = 8.5.;




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

  int thresholdEvents = 1;


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

    psData.SetBaseEvents(baseEvents);
    psData.SetSource(thisSignal);
    psData.SetRandomizeSrc(false);  // in case true was set in another script

    double meanSrcEv_ForDetection;
    double meanFlux_ForDetection;

    // are there any signal events at all in this range?
    TString eString = "mcPrimary_Energy_GeV > " + TStringify(eMinVect[i]);
    eString += " && mcPrimary_Energy_GeV < " + TStringify(eMaxVect[i]);
    TCut eCut = eString;
    int evPass = 
      srcTree->Draw("1", evLoader.GetCuts()*zenRangeCut*pullRangeCut*eCut);
    cout << evPass << " good events in range:  " << eString << endl;

    if ( evPass > thresholdEvents ) {

      DiscoveryPotential disco;
      disco.monitor_ = true;
      disco.method_ = 2;
      disco.loops_ = 20;
      cout << "LOOPS = " << disco.loops_ << endl;
      
      LlhEnergy llhEnergyFn;
      llhEnergyFn.SetUseEnergy(zenUseEnergy);
      llhEnergyFn.SetOptimizeTolerance(.01);
      llhEnergyFn->SetMonitorLevel(0);
      disco.SetLogLikelihood(llhEnergyFn);
      
      if (discovery == true) {
	cout << "Setting for discovery\n";
	disco.SetForDiscovery();
	disco.SetDetectionSignificance(2.87e-7); // one-sided p-value for 5sigma
      } else {
	cout << "\nSetting for Neyman upper limit... check script\n\n";
	disco.SetForMedianUpperLimit(1000); // ntrials
	cout << "Calculating Median Bkg Prob (set sensitivity p-threshold):\n";
      }

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


    TLine *eventLine = new TLine(log10(eMinVect[i]),meanSrcEv_ForDetection,
				   log10(eMaxVect[i]),meanSrcEv_ForDetection);
    can->cd(1);
    eventLine->Draw();


    // covert E^-2 flux from GeV to TeV
    //    double scale = ScaleFlux(1e9, 1e12, 1e9, 1e12, indexVect[i]);
    //    double yValue = scale*meanFlux_ForDetection;
    //    TLine *fluxLine = new TLine(log10(eMinVect[i]),yValue,
    //				log10(eMaxVect[i]),yValue);


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
