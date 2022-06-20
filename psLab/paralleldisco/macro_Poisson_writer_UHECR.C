 
{
  bool optUseEnergy = true;

  // These now set by rootSubmitMacro_stack.C:
  int nSigMin;
  int nSigMax;
  int nTrials;
  long ranSeedGlobal;
  bool isUHECR;

  char *mainDirName = resultsDir;  // SET BY rootSubmitMacro.C

  cout << "Setting up analysis for: " << srcName << endl; 
  cout << "is it UHECR? "<<isUHECR<<endl;

  // SET UP TO WORK WITH MULTI-ANALYSIS

  char resultDirName[500], notesFileName[500], sigFileName[500];
  sprintf(resultDirName,"Poisson_%s/meanNSig_%03d/",srcName,nSigMin);

  gSystem->mkdir(resultDirName,1);   // okay if already exists, recursive

  double maxLlh;
  double bestFitNs;
  double bestFitGamma = 0.;

  cout << "nSig = ";

  int nTest=0;
  for (int nSig = nSigMin; nSig<=nSigMax; ++nSig) {
    cout << nSig << ", " << flush;
    cout<<endl;
    sprintf(sigFileName,"%s/08%d.txt",resultDirName,-1*ranSeedGlobal);
    FILE *fpSig = fopen(sigFileName,"w");
 
    for (int nt=0; nt<nTrials; ++nt) {
      nTest = random_poisson(nSig);
  	double bestTS=0.0;
  	double bestPar0=0.0;
  	double bestPar1=0.0;
	double bestSig=-1.0;

        for (double mysig=0.0; mysig<5.1; mysig=mysig+0.5){
              //first setup the sigma vector used on the llh side as mysig
              vector<double> mySigmas;
              mySigmas.clear();
              for (int j=0; j<srcLocations.size(); j++)  mySigmas.push_back(mysig);
   
              // NEW FOR EXT 
              newllhstack40.SetSourceSigmas(mySigmas);
              newllhstack59.SetSourceSigmas(mySigmas);
   
              MultiAnalysisFn * maf2 = new MultiAnalysisFn;
              maf2->AddAnalysisFn(&newllhstack40);
              maf2->AddAnalysisFn(&newllhstack59);
              maf2->SetParTranslator(&pt);
              mas.GenerateDataSet_with_nSrcEvents(nTest);
   
              vector<MinuitParDef> pdv;
              pdv.push_back( MinuitParDef("nSrc",2,0.1, 0.,100.) );
              pdv.push_back( MinuitParDef("gamma",2.5,0.5, 1., 4.) );
              maf2->SetParDefs(pdv);
   
              maf2->MaximizeLlh();
              if(maf2->GetTestStatistic()>bestTS) {
		//GetTestStatistic() is equivallent to Get_logLambdaBest()
                      bestTS = maf2->GetTestStatistic();
                      bestPar0 = maf2->GetPar(0);
                      bestPar1 = maf2->GetPar(1);
                      bestSig  = mysig;
              }
              delete maf2;
   
         }//end of sigma loop
	 fprintf(fpSig,"%d %d %lg %.3f %.3f %.3f\n",nSig,nTest,bestTS,bestPar0,bestPar1, bestSig);
	 cout<<nSig<<" "<<nTest<<" "<<bestTS<<" "<<bestPar0<<" "<<bestPar1<<" "<<bestSig<<endl;
    }

    fclose(fpSig);
  }

  cout << "done\n";

 END:
}
