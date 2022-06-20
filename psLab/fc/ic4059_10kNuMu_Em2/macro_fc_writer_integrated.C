 
{
  
  //nSigMin, nSigMax, nTrials and resultsDir set in rootSubmitMacro, which calls this.
  
  bool optUseEnergy = true;

  gSystem->cd("$LAB_MAIN_DIR/macro_llh/ic59");

  gROOT->ProcessLine(".x $LAB_MAIN_DIR/macro_llh/ic59/Track_loadmulti_4059.C");
  // this load the Ark script, the IC40+59 multiArk (mark) and the multiAnalysisFn (maf)
  //initialize_ran1(ranSeedI); // set up to initialize with -55

  gSystem->cd(startDir);

  //int nSigMin = 0;
  //int nSigMax = 50;

  //int nTrials = 10000;//000;

  double spectralIndex = -2;


  char *mainDirName = resultsDir;  // SET BY rootSubmitMacro.C

  // double srcDecDeg = 11.375;
  // srcDecDeg IS SET BY rootSubmitMacro.C

  //I3Analysis &aSet = psData;
  //EventLoader &evLoader = ns; // Already the correct name

  // LOAD SOURCE EVENTS

  //vector<I3Event> srcEvents = mark.psData.GetSource()sourceEvents;
  
  EquatorialDeg srcLocation(0.,decDeg);
  //evLoader.SetSourceEqDeg(srcLocation); // Old
  //evLoader.LoadSourceEvents(srcEvents, srcLocation);
  
  //if (baseEvents.size() == 0 || srcEvents.size() == 0) {
  //  cout << "ERROR: ZERO EVENTS!  That's going to be a problem...\n";
  //  cout << "Forcing macro to stop now...\n\n";
    // force root to stop macro here:
  //  cout << 1/0 << endl;
  //}

  // OLD SOURCE
  //NugenSource fcSource;
  //fcSource.SetSourceParams(flux, spectralIndex, livetime);

  //TimePdf * tPdf = new PeriodicGaussianTimePdf(0.,1.,0.5,sigma,1.);

  mark.SetPointSource(srcLocation, PowerLawFlux(1.,spectralIndex));
  maf.SetSearchCoord(srcLocation);

  // NEW SIGNAL
  //PowerLawFlux fcFlux(1,spectralIndex);  //  1 GeV^-1 cm^-2 s^-1 , index
  //I3PointGenerator fcSignal(srcEvents, fcFlux,
  //            srcLocation, tPdf, livetime);
  //cout << "Power Law Flux with spectralIndex = " << spectralIndex << endl;


  //aSet.SetBaseEvents(bkgEvents);
  //aSet.SetBaseEvents(baseEvents);

  //fcSignal.StoreCandidateEvents(srcLocation, srcEvents); # Old
  //aSet.SetSource(fcSignal);
  //aSet.SetRandomizeSrc(false);  // in case true was set in another script

  

  char decDirName[500], notesFileName[500], sigFileName[500];
  sprintf(decDirName,"%s/zenith_%d_%d_integrated",mainDirName,decDeg,-1.0*spectralIndex);
  sprintf(notesFileName,"%s/notes.txt",decDirName);

  cout << decDirName << endl;

  gSystem->MakeDirectory(mainDirName);  // okay if already exists
  gSystem->MakeDirectory(decDirName);   // okay if already exists

  FILE *fp = fopen(notesFileName,"w");
  if (fp) { cout << "Writing notes to " << notesFileName << endl;
  }
  else { 
    cout << "Could not open " << notesFileName << endl;
    goto END;
  }

  cout << nSigMin << " " << nSigMax << endl;

  fprintf(fp,"nSigMin , nSigMax: %d , %d\n",nSigMin,nSigMax);
  //fprintf(fp,"src , sigma: %s, %d\n",srcName[nsrc],sigma);
  fprintf(fp,"nTrials: %d\n",nTrials);
  fprintf(fp,"SpectralIndex: %f\n",spectralIndex);
  //fprintf(fp,"Zenith Range (+/- deg): %f\n",sourceZenWidthDeg);
  //fprintf(fp,"Signal events loaded (passed cuts): %d\n", srcEvents.size());
  fprintf(fp,"Cuts:  %s\n", ark59.evLoader.GetCuts().GetTitle());
  fclose(fp);

  //if (srcEvents.size()==0) { 
  //  cout << "ERROR: no signal events loaded.\n";
  //  goto END;
  //}

/*
  int TimePdfType = 0;

  bool useE = true;

  NewLlhEnergy llhEnergyFn;
  llhEnergyFn.SetUseEnergy(useE);
  llhEnergyFn.SetOptimizeTolerance(0.01);
  llhEnergyFn->SetMonitorLevel(0);
  llhEnergyFn.SpectralPenalty = false;
  llhEnergyFn.ndof = 2.;
  
  AnalysisLlh &llhFn = llhEnergyFn; */


  double maxLlh;
  double bestFitNs;
  double bestFitGamma = 0.;
//  double bestFitPeriod;
//  double bestFitSigma;

  cout << "nSig = ";

  for (int nSig = nSigMin; nSig<=nSigMax; ++nSig) {
    cout << nSig << ", " << flush;
    sprintf(sigFileName,"%s/nSig_%03d.txt",decDirName,nSig);
    FILE *fpSig = fopen(sigFileName,"w");
    
    for (int nt=0; nt<nTrials; ++nt) {
     
      mark.psData->GenerateDataSet_with_nSrcEvents(nSig);
      maf.MaximizeLlh();
      
      maxLlh = maf.Get_logLambdaBest();
      bestFitNs = maf.GetPar(0);
      if (optUseEnergy) {
	bestFitGamma = maf.GetPar(1);
      }
      
      //  if (bestFitNs<0.) { maxLlh = -maxLlh; }
      fprintf(fpSig,"%lg %.3f %.3f\n",maxLlh,bestFitNs,bestFitGamma);
    }

    fclose(fpSig);
  }

  cout << "done\n";

 END:
}
