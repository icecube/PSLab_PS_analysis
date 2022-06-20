 
{
  bool optUseEnergy = true;

  // These now set by rootSubmitMacro_stack.C:
  int nSigMin;
  int nSigMax;
  int nTrials;
  long ranSeedGlobal;

  char *mainDirName = resultsDir;  // SET BY rootSubmitMacro.C

  //double spectralIndex = -2;
  cout << "Setting up analysis for: " << srcName << endl; 
  //gROOT->ProcessLine(".x ic40_fix_llh6_Ext/Milagro/macro_SetupMilagroAnalysis.C");
/*
  if (srcName == "GP_FermiGalDiffuse") {
    gROOT->ProcessLine(".x $LAB_MAIN_DIR/macro_llh/ic59/galplane/GP_FermiGalDiffuse/Track_multi_4059.C");
  }
*/

  // SET UP TO WORK WITH MULTI-ANALYSIS

  char resultDirName[500], notesFileName[500], sigFileName[500];
  //sprintf(resultDirName,"%s/Poisson_%s/meanNSig_%03d/",mainDirName,srcName,nSigMin);
  sprintf(resultDirName,"Poisson_%s/meanNSig_%03d/",srcName,nSigMin);
  //sprintf(notesFileName,"%s/notes.txt",resultDirName);

  //gSystem->MakeDirectory(mainDirName);  // okay if already exists
  //gSystem->MakeDirectory(resultDirName);   // okay if already exists, recursive
  gSystem->mkdir(resultDirName,1);   // okay if already exists, recursive

  double maxLlh;
  double bestFitNs;
  double bestFitGamma = 0.;

  cout << "nSig = ";

  int nTest=0;
  for (int nSig = nSigMin; nSig<=nSigMax; ++nSig) {
    cout << nSig << ", " << flush;
    sprintf(sigFileName,"%s/08%d.txt",resultDirName,-1*ranSeedGlobal);
    FILE *fpSig = fopen(sigFileName,"w");
 
    for (int nt=0; nt<nTrials; ++nt) {
      nTest = random_poisson(nSig);
      mas.GenerateDataSet_with_nSrcEvents(nTest);

      maf.MaximizeLlh();
      maxLlh = maf.Get_logLambdaBest();
      bestFitNs = maf.GetPar(0);
      if (optUseEnergy) {
      	bestFitGamma = maf.GetPar(1);
      }
      //  if (bestFitNs<0.) { maxLlh = -maxLlh; }
      fprintf(fpSig,"%d %d %lg %.3f %.3f\n",nSig,nTest,maxLlh,bestFitNs,bestFitGamma);
    }

    fclose(fpSig);
  }

  cout << "done\n";

 END:
}
