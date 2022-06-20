
{
  /* ASSUME WE HAVE ALREADY SET UP:

  I3Analysis psData;
  // Smoothing radius has to already be set
  // vector<I3Event> baseEvents  already defined

  EventLoader evLoader;
  // Set up to repeatedly load source events

  FOR SIGNAL:
  double spectralIndex;
  double livetime;

  */

  gROOT->ProcessLine(".x fc/loadlibs.C");


  // Some Options

  double bkgBinSize = 2.0; // degrees

  bool MACRO_DISCOVERY_ZENITH_SinDec = true;

  double zenithDegLT = 175.; // evaluate only for zenith less than this


  // 1k NuMu, no NuTaus
  //char fcDirectory[1000] = "/net/user/jdumm/testlab/fc/ic40_sigUpgrade_Em2/results/";
  // 10k NuMu, 10k NuTau
  //char fcDirectory[1000] = "/net/user/jdumm/testlab/fc/ic40_10kNuMu_5kNuTau_Em2/results/";
  // 10k NuMu, No NuTau
  //char fcDirectory[1000] = "/net/user/jdumm/testlab/fc/ic40_10kNuMu_NoNuTau_Em2/results/";
  char fcDirectory[1000] = "$LAB_MAIN_DIR/fc/results/ic4059_10kNuMu_NoNuTau_Em2";
  double spectralIndex = -2.0;
  bool optMedian = true;

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat("");
  TCanvas *can = new TCanvas("CanZenith","Discovery Potential in Zenith",
			     20,20,800,700);
  can->Divide(2,2,0.005,0.005);

  double ListDec[1000];
  int NDec;


  // Declination
  if (!MACRO_DISCOVERY_ZENITH_SinDec) {
    double decMin = 0.;
    double decMax = 85.;
    NDec=17;

    for (int iDec=0; iDec<NDec; ++iDec) {
      ListDec[iDec] = (iDec+0.5)*(decMax-decMin)/NDec + decMin;
    }
    TH1D *hBkgEvents = new TH1D("hBkgEvents","Bkg Events in "+TStringify(bkgBinSize)+"#circ radius bin ; Declination [#circ]",NDec,decMin,decMax);
    TH1D *hFixedFluxEvents = new TH1D("hFixedFluxEvents","Mean Source Events for Fixed Flux ; Declination [#circ]",NDec,decMin,decMax);
    TH1D *hDetectEvents = new TH1D("hDetectEvents","Mean Source Events for Detectable Flux ; Declination [#circ]",NDec,decMin,decMax);
    TH1D *hDetectFlux = new TH1D("hDetectFlux","Detectable Flux ; Declination [#circ]",NDec,decMin,decMax);
  }


  if (MACRO_DISCOVERY_ZENITH_SinDec) {
    double sinDecMin = -1;
    double sinDecMax = 1.;
    NDec=40;
    //NDec=200;

    for (int iDec=0; iDec<NDec; ++iDec) {
      ListDec[iDec] = asin( (iDec+0.5)*(sinDecMax-sinDecMin)/NDec + sinDecMin)
	*TMath::RadToDeg();
    }
    TH1D *hBkgEvents = new TH1D("hBkgEvents","Bkg Events in "+TStringify(bkgBinSize)+"#circ radius bin ; sin(#delta)",NDec,sinDecMin,sinDecMax);
    TH1D *hFixedFluxEvents = new TH1D("hFixedFluxEvents","Mean Source Events for Fixed Flux ; sin(#delta)",NDec,sinDecMin,sinDecMax);
    TH1D *hDetectEvents = new TH1D("hDetectEvents","Mean Source Events for Detectable Flux ; sin(#delta)",NDec,sinDecMin, sinDecMax);
    TH1D *hDetectFlux = new TH1D("hDetectFlux","Detectable Flux ; sin(#delta)",NDec,sinDecMin,sinDecMax);
  }


  //
  // DECLINATION LOOP
  //


  for (int iDec=0; iDec<NDec && ListDec[iDec]+90.<zenithDegLT; ++iDec) {

    double srcDec = ListDec[iDec];

    EquatorialDeg scanSrcLocation(0.,srcDec);
    cout << "Source position: " << scanSrcLocation.GetRa() << "deg r.a.,  ";
    cout << scanSrcLocation.GetDec() << "deg dec\n";

    // LOAD SOURCE EVENTS

    vector<I3Event> scanSrcEvents;
    // MERGED FUNCTIONS... //   evLoader.SetSourceEqDeg(scanSrcLocation);
    evLoader.LoadSourceEvents(scanSrcEvents, scanSrcLocation);
    if (baseEvents.size() == 0 || scanSrcEvents.size() == 0) {
      cout << "ERROR: ZERO EVENTS!  That's going to be a problem...\n";
      cout << "Forcing macro to stop now...\n\n";
      // force root to stop macro here:
      cout << 1/0 << endl;
    }

    // OLD SOURCE
    //    mySource.SetSourceParams(flux, spectralIndex, livetime);
    //    mySource.StoreCandidateEvents(scanSrcLocation, scanSrcEvents);

    // NEW SIGNAL
    PowerLawFlux scanFlux(1,spectralIndex);  //  1 GeV^-1 cm^-2 s^-1 , index
    I3PointGenerator scanSignal(scanSrcEvents, scanFlux, 
				scanSrcLocation, livetime);
    cout << "Power Law Flux with spectralIndex = " << spectralIndex << endl;

    double meanSrcEv_ForInitialFlux = scanSignal.GetMeanSrcNev();
    cout << "  Mean Src Events for Init. Flux: ";
    cout << meanSrcEv_ForInitialFlux << endl;


    psData.SetBaseEvents(baseEvents);
    psData.SetSource(scanSignal);
    psData.SetRandomizeSrc(false);  // in case true was set in another script

    // Have to do this before calculating number density,
    // since number density includes all events, even if src is added
    // (hmm, maybe should remove the term 'bkg' from that function call)
    psData.GenerateDataSet_with_nSrcEvents(0); 

    double bkgIn2Bin = psData.BkgNumberDensity(scanSrcLocation)
      * bkgBinSize * bkgBinSize * TMath::Pi();

    if (MACRO_DISCOVERY_ZENITH_SinDec) {
      hBkgEvents->Fill(sin(srcDec*TMath::DegToRad()),bkgIn2Bin);
      hFixedFluxEvents->Fill(sin(srcDec*TMath::DegToRad()),
			     meanSrcEv_ForInitialFlux);
    } else {
      hBkgEvents->Fill(srcDec,bkgIn2Bin);
      hFixedFluxEvents->Fill(srcDec,meanSrcEv_ForInitialFlux);
    }

    can->cd(1);
    hBkgEvents->Draw();
    can->cd(2);
    hFixedFluxEvents->Draw();
    can->Update();

    // FC STUFF
    double sigUnc = 0.15;  // systematic unc. on signal

    cout << "SIGNED SQRT FC\n";
    SignedSqrtFC signedfc;
    signedfc.SetParams(800,-10.,30.,410, sigUnc);
    GeneralFC &gfc = signedfc;


/*
    DiscoveryPotential disco;
    disco.monitor_ = true;

    disco.method_ = 2;
    disco.loops_ = 50;
    //disco.loops_ = 20;
*/

    NewLlhEnergy llhEnergyFn;
    bool zenUseEnergy = true;
    llhEnergyFn.SetUseEnergy(zenUseEnergy);
    llhEnergyFn.SetOptimizeAngleDeg(10);
    llhEnergyFn.SetOptimizeTolerance(.01);
    //llhEnergyFn.SetOptimizeTolerance(0); // No optimization done!
    llhEnergyFn->SetMonitorLevel(0);
    //disco.SetLogLikelihood(llhEnergyFn);

    /*
    BinnedAnalysis binFn;
    binFn.SetBinSize(bkgBinSize);
    TProfile profileLogProb(MakeNewName("profileLogProb"),
			    "Binned Analysis;Bin Size [#circ];Log_{10} Prob",
			    50,0.,10.,
			    "S"); // RMS spread option
    binFn.SetTProfile(&profileLogProb);
    disco.SetLogLikelihood(binFn);
    */


    // DISCOVERY  OR  UPPER-LIMIT

/*
    if (SUPER_MACRO_DISCOVERY) { // Select for Discovery Potential    
      disco.SetForDiscovery();

      //disco.detectionSignificance_ = 2.87e-7; // one-sided p-value for 5sigma
      disco.SetDetectionSignificance(2.87e-7); // one-sided p-value for 5sigma

      disco.SetDetectionPower(0.5);
    } else { // Select for Sensitivity
      //disco.SetForUpperLimit();
      disco.SetForMedianUpperLimit(1000);
      disco.SetDetectionPower(0.90);
      cout << "Calculating Median Bkg Prob (set sensitivity p-threshold):\n";
      double medBkgProb = disco.MedianProbForBkg(&psData, srcLocation, 1000);
      //disco.detectionSignificance_ = medBkgProb;
      disco.SetDetectionSignificance(medBkgProb);
      cout << " Median Bkg Prob: " << medBkgProb << endl;
    }


    disco.AnalyzeDiscoveryPotential(&psData,scanSrcLocation);
*/

    // FC sensitivities
    SetDecGeneralFC(gfc,fcDirectory, 50, 30., srcDec);

    double confLev = 0.90;
    signedfc.BuildConfidenceBand(confLev);

    //double nSrcUpLimit = signedfc.GetUpperLimit(sqrt(2*fabs(maxLlh)));
    double meanSrcEv_ForDetection;

    //optMedian = true;
    if (optMedian) {
      meanSrcEv_ForDetection = gfc.GetMedianUpperLimit();
    } else {
      meanSrcEv_ForDetection = gfc.GetAverageUpperLimit();
    }

    double meanFlux_ForDetection =
      scanSignal.GetFluxScaleForNev(meanSrcEv_ForDetection);


    vector<I3Event> sourceEvents;
    evLoader.SetMonitor(false);
    evLoader.LoadSourceEvents(sourceEvents, scanSrcLocation);

    if (MACRO_DISCOVERY_ZENITH_SinDec) {
      double sinSrcDec = sin(srcDec*TMath::DegToRad());
      hDetectEvents->Fill(sinSrcDec, meanSrcEv_ForDetection);
      hDetectFlux->Fill(sinSrcDec, meanFlux_ForDetection);
    } else {
      hDetectEvents->Fill(srcDec, meanSrcEv_ForDetection);
      hDetectFlux->Fill(srcDec, meanFlux_ForDetection);
    }

    can->cd(3);
    hDetectEvents->Draw();
    can->Update();

    can->cd(4);
    hDetectFlux->Draw();
    can->Update();

    if (optMedian) {
      cout << confLev*100 << "% median FC upper limits\n";
    } else {
      cout << confLev*100 << "% average FC upper limits\n";
    }
    cout << meanSrcEv_ForDetection << " events, ";
    cout << meanFlux_ForDetection << " flux\n";
    cout << "Results for spectral index= "<<scanFlux.GetSpectralIndex()<<endl;
    cout << "Livetime (days) = " << scanSignal.GetLivetime()/86400. << endl;
    //    cout << "\nUsing energy info = " << zenUseEnergy << endl;
    cout << "Running average flux: " << hDetectFlux->GetSum()/(iDec+1) << endl;



/*
    cout << "Detection Power " << disco.GetDetectionPower()*100 << "%\n"; 
    cout << "Detection Significance " << disco.GetDetectionSignificance() << endl;
    cout << meanSrcEv_ForDetection << " events, ";
    cout << meanFlux_ForDetection << " flux\n";
    cout << "Results for spectral index= "<<scanFlux.GetSpectralIndex()<<endl;
    cout << "Livetime (days) = " << scanSignal.GetLivetime()/86400. << endl;
    cout << "\nUsing energy info = " << zenUseEnergy << endl;
    cout << "Running average flux: " << hDetectFlux->GetSum()/(iDec+1) << endl;
*/
  }


}
