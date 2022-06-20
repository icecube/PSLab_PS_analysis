 
{
  bool optUseEnergy = true;

  int nSigMin = 0;
  int nSigMax = 50;

  int nTrials = 10000;

  double spectralIndex = -2;


  char *mainDirName = resultsDir;  // SET BY rootSubmitMacro.C

  // double srcDecDeg = 11.375;
  // srcDecDeg IS SET BY rootSubmitMacro.C


  I3Analysis &aSet = psData;
  //EventLoader &evLoader = ns; // Already the correct name

  // LOAD SOURCE EVENTS

  vector<I3Event> srcEvents;
  EquatorialDeg srcLocation(0.,srcDecDeg);
  //evLoader.SetSourceEqDeg(srcLocation); // Old
  evLoader.LoadSourceEvents(srcEvents, srcLocation);

  if (baseEvents.size() == 0 || srcEvents.size() == 0) {
    cout << "ERROR: ZERO EVENTS!  That's going to be a problem...\n";
    cout << "Forcing macro to stop now...\n\n";
    // force root to stop macro here:
    cout << 1/0 << endl;
  }

  // OLD SOURCE
  //NugenSource fcSource;
  //fcSource.SetSourceParams(flux, spectralIndex, livetime);

  // NEW SIGNAL
  PowerLawFlux fcFlux(1,spectralIndex);  //  1 GeV^-1 cm^-2 s^-1 , index
  I3PointGenerator fcSignal(srcEvents, fcFlux,
              srcLocation, livetime);
  cout << "Power Law Flux with spectralIndex = " << spectralIndex << endl;


  //aSet.SetBaseEvents(bkgEvents);
  aSet.SetBaseEvents(baseEvents);

  //fcSignal.StoreCandidateEvents(srcLocation, srcEvents); # Old
  aSet.SetSource(fcSignal);
  aSet.SetRandomizeSrc(false);  // in case true was set in another script


  char decDirName[500], notesFileName[500], sigFileName[500];
  sprintf(decDirName,"%s/dec_%08.4f/",mainDirName,srcDecDeg);
  sprintf(notesFileName,"%s/notes.txt",decDirName);

  gSystem->MakeDirectory(mainDirName);  // okay if already exists
  gSystem->MakeDirectory(decDirName);   // okay if already exists

  FILE *fp = fopen(notesFileName,"w");
  if (fp) { cout << "Writing notes to " << notesFileName << endl;
  }
  else { 
    cout << "Could not open " << notesFileName << endl;
    goto END;
  }

  fprintf(fp,"nSigMin , nSigMax: %d , %d\n",nSigMin,nSigMax);
  fprintf(fp,"nTrials: %d\n",nTrials);
  fprintf(fp,"SpectralIndex: %f\n",spectralIndex);
  fprintf(fp,"Zenith Range (+/- deg): %f\n",sourceZenWidthDeg);
  fprintf(fp,"Signal events loaded (passed cuts): %d\n", srcEvents.size());
  fprintf(fp,"Cuts:  %s\n", evLoader.GetCuts().GetTitle());
  fclose(fp);

  if (srcEvents.size()==0) { 
    cout << "ERROR: no signal events loaded.\n";
    goto END;
  }


  LlhEnergy llhEnergyFn;
  llhEnergyFn.SetUseEnergy(optUseEnergy);
  llhEnergyFn.SetOptimizeAngleDeg(10);
  llhEnergyFn.SetOptimizeTolerance(.01);
  llhEnergyFn.SetMonitorLevel(0);

  AnalysisLlh &llhFn = llhEnergyFn;


  double maxLlh;
  double bestFitNs;
  double bestFitGamma = 0.;

  cout << "nSig = ";

  for (int nSig = nSigMin; nSig<=nSigMax; ++nSig) {
    cout << nSig << ", " << flush;
    sprintf(sigFileName,"%s/nSig_%03d.txt",decDirName,nSig);
    FILE *fpSig = fopen(sigFileName,"w");
 
    for (int nt=0; nt<nTrials; ++nt) {

      aSet.GenerateDataSet_with_nSrcEvents(nSig);
      llhFn.SetAnalysis(aSet, srcLocation);

      llhFn.MaximizeLlh();
      maxLlh = llhFn.Get_logLambdaBest();
      bestFitNs = llhFn.GetPar(0);
      if (optUseEnergy) {
	bestFitGamma = llhFn.GetPar(1);
      }
      //  if (bestFitNs<0.) { maxLlh = -maxLlh; }
      fprintf(fpSig,"%lg %.3f %.3f\n",maxLlh,bestFitNs,bestFitGamma);
    }

    fclose(fpSig);
  }

  cout << "done\n";

 END:
}
