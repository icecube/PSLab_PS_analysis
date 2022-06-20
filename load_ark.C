void load_ark(I3Ark& ark, bool opt_UseRealData, TString period, TString RecoName, bool useAngFloor, bool useNewDataSet=true) {

  Printf("loading period %s",period.Data());

  //
  // CONFIGURE EVENTLOADER WITH BACKGROUND / REAL DATA SAMPLE
  //

  vector<double> startVect, stopVect;

  cout<<"Using new data structure? "<<useNewDataSet<<endl;
  ark.evLoader.SetNewDataStructure(useNewDataSet);
  ark.evLoader.SetApplyAngFloor(useAngFloor);
  ark.evLoader.SetRndOnlyTime(false);
  ark.evLoader.SetBkgTree( LoadTree_GoodRuns_Full(period, startVect, stopVect) );
  cout << "Loaded tree" << endl;

  ark.livetime = GetValueFromTree(ark.evLoader.GetBkgTree(),"livetimeTotal");
  ark.tmin = GetValueFromTree(ark.evLoader.GetBkgTree(),"tmin");
  ark.tmax = GetValueFromTree(ark.evLoader.GetBkgTree(),"tmax"); 
  ark.startMissRuns = startVect;
  ark.stopMissRuns = stopVect;
 
  ark.evLoader.SetBkgLoadMethod_Exact();
    
  if (opt_UseRealData) {
    ark.evLoader.SetTimeMethod_Actual("timeMJD"); 
    cout << "Will load  * * * * U N B L I N D E D * * * *  Base Data!\n";
  }
  else {
    
    ark.evLoader.SetTimeMethod_Scramble();  // for blindness!
    cout << "Will load scrambled times for Base Data.\n";
  }
  
  cout << "Livetime (Days): " << ark.livetime/86400. << "\n";
    
  //
  // CONFIGURE EVENTLOADER WITH SIGNAL SIMULATION DATA SAMPLE
  //
  // evLoader will load signal events for any specified declination

  cout << "Configuring Source Event Sample...\n";

  ark.sourceZenWidthDeg = 0.5;
  ark.evLoader.SetSourceZenWidthDeg(ark.sourceZenWidthDeg);
 
  
  //
  // SPECIFY CUTS, if any
  //
  //TCut IC86_Cut = "";//no cut
  // Cuts can be added or reset any time, and events re-loaded with new cuts
  // (cuts are also evaluated in order, for possible efficiency gain)
  //ark.evLoader.AddCut(IC86_Cut);
  

  //
  // SET NAMES OF VARIABLES TO BE USED FROM ROOT FILES
  //

  TString recoZenRadName = "zenRad"; 
  TString recoAziRadName = "aziRad"; 
  TString energyVar      = "logMuexEn";  
  ark.evLoader.SetName_recoZenith_rad(recoZenRadName);
  ark.evLoader.SetName_recoAzimuth_rad(recoAziRadName);
  ark.evLoader.SetName_energyValue(energyVar);
  ark.evLoader.SetName_eventID("EventID");
  ark.evLoader.SetName_runID("RunID");
  ark.evLoader.SetName_sigmaDeg("sigmaDeg");

  //
  // LOAD EVENTS, SET UP ENERGY PDFs
  // 
  cout << "Loading Background Events: " << ark.evLoader.GetBkgTree()->GetTitle() << endl;
  ark.evLoader.SetMonitor(true);
  ark.evLoader.LoadBkgEvents(ark.baseEvents);
  
  ZenithEnergyProb* zen_eProb = new ZenithEnergyProb();
  
  cout << "Filling Energy PDFs:\n";
  zen_eProb->SetSourceZenWidthDeg( ark.sourceZenWidthDeg);
  zen_eProb->SetName_recoZenith_rad(recoZenRadName);

  vector<double> zenMinDegVect;
      
  // Fill vector with bin edges that match "CutDMS" bin edges
  // 0 - 90  zen added
  double tempBot;
  if(period=="IC40") {
    tempBot=1;
    for(int i=0; i<21;i++){
      zenMinDegVect.push_back(acos(tempBot)*TMath::RadToDeg());
      tempBot-=0.05;
    }
    zenMinDegVect.push_back(110.00); // ~acos(-0.33)
    zenMinDegVect.push_back(132.00); // ~acos(-0.67)
    zenMinDegVect.push_back(180.);
  }
  else if(period=="IC59") {
    tempBot=1;
    for(int i=0; i<21;i++){
      zenMinDegVect.push_back(acos(tempBot)*TMath::RadToDeg());
      tempBot-=0.05;
    }
    zenMinDegVect.push_back(110.00); // ~acos(-0.33)
    zenMinDegVect.push_back(132.00); // ~acos(-0.67)
    zenMinDegVect.push_back(180.);
  } 
  else if(period=="IC79") {
    tempBot=1;
    for(int i=0; i < 19;i++){
      zenMinDegVect.push_back(acos(tempBot)*TMath::RadToDeg());
      tempBot-=0.05;
    }
    zenMinDegVect.push_back(86.00); 
    zenMinDegVect.push_back(88.00); 
    zenMinDegVect.push_back(90.00);
    zenMinDegVect.push_back(97.); // Follow BDT regions
    zenMinDegVect.push_back(125.);
    zenMinDegVect.push_back(180.);
  }
  else if(period=="IC86_I") {
    zenMinDegVect.push_back(acos(1.)*TMath::RadToDeg());
    zenMinDegVect.push_back(acos(0.97)*TMath::RadToDeg());
    zenMinDegVect.push_back(acos(0.94)*TMath::RadToDeg());
    tempBot=0.9;
    for(int i=0; i<17;i++){
      zenMinDegVect.push_back(acos(tempBot)*TMath::RadToDeg());
      tempBot-=0.05;
    }
    zenMinDegVect.push_back(86.00); 
    zenMinDegVect.push_back(88.00); 
    tempBot-=0.05;
    for(int i=0; i<5;i++){
      zenMinDegVect.push_back(acos(tempBot)*TMath::RadToDeg());
      tempBot-=0.2;
    }
    zenMinDegVect.push_back(180.00);

  }
  else if(period=="IC86_II_VII" || period == "IC86_II" || period == "IC86_III" || period == "IC86_IV" || period == "IC86_V" || period == "IC86_VI" || period == "IC86_VII") {
    zenMinDegVect.push_back(acos(1.)*TMath::RadToDeg());
    zenMinDegVect.push_back(acos(0.97)*TMath::RadToDeg());
    zenMinDegVect.push_back(acos(0.94)*TMath::RadToDeg());
    tempBot=0.9;
    for(int i=0; i<17;i++){
      cout << tempBot << " " << acos(tempBot)*TMath::RadToDeg() << endl;
      zenMinDegVect.push_back(acos(tempBot)*TMath::RadToDeg());
      tempBot-=0.05;
    }
    zenMinDegVect.push_back(86.00); 
    zenMinDegVect.push_back(88.00); 
    tempBot-=0.05;
    for(int i=0; i<5;i++){
      cout << tempBot << " " << acos(tempBot)*TMath::RadToDeg() << endl;
      zenMinDegVect.push_back(acos(tempBot)*TMath::RadToDeg());
      tempBot-=0.2;
    }
    zenMinDegVect.push_back(180.00);
  }

  zen_eProb->SetZenithBandsDeg(zenMinDegVect);
  zen_eProb->SetLoadModeNew(true); // true is now faster

  if (energyVar == "logMuexEn") {
    int nBackFill = 35; // don't backfill previous bins
    cout << "backfill" << endl;
    zen_eProb->SetEnergyGammaRangeAndBackFill(40,2.,9., 30,1.,4., nBackFill);
  }

  zen_eProb->SetTableBkg(ark.baseEvents);
  TStopwatch ts;
  cout << "set table gamma" << endl;
  zen_eProb->SetTableGamma(ark.GetEffectiveAreaDistribution(), ark.GetRecoEnergyDistribution());
  ts.Print();

  ark.eProb = zen_eProb;
    
  ark.decBkgProb.Initialize(180, 1);
  ark.decBkgProb.SetBaseDecMap(ark.baseEvents);
  ark.lcBkgProb.Initialize(12.,90.,true); //so why doesn't this work together?
  ark.lcBkgProb.FillLCBkgHisto(ark.baseEvents);
   
  I3Analysis *psData = new I3Analysis();
  psData->SetRndOnlyTime(false);
 
  psData->SetBkgSpaceProb(ark.decBkgProb);
  psData->SetBaseEvents(ark.baseEvents);
  psData->SetEnergyProb(*(ark.eProb));
  
  EventTimeModule * ic_times = new EventTimeModuleDiscrete();

  if (opt_UseRealData) {
    psData->UseRealData();  // the event set is now exactly equal
    // to the data set (i.e. no scrambling, no fake signal added.)
  } else {
    psData->SetRandomizeBase(true); 
    if(period=="IC86_II_VII") ic_times->SetTimesFromMJDFile("ListOfTimes/HugeListOfTimes_IC86-II-III-IV-V-VI-VII_GRL_version-003-p02.txt");
    else if(period=="IC40") ic_times->SetTimesFromMJDFile("ListOfTimes/HugeListOfTimes_IC40_GRL_version-003-p02.txt");
    else if(period=="IC59") ic_times->SetTimesFromMJDFile("ListOfTimes/HugeListOfTimes_IC59_GRL_version-003-p02.txt");
    else if(period=="IC79") ic_times->SetTimesFromMJDFile("ListOfTimes/HugeListOfTimes_IC79_GRL_version-003-p02.txt");
    else if(period=="IC86_I") ic_times->SetTimesFromMJDFile("ListOfTimes/HugeListOfTimes_IC86-I_GRL_version-003-p02.txt");
    else if(period=="IC86_II") ic_times->SetTimesFromMJDFile("ListOfTimes/HugeListOfTimes_IC86-II_GRL_version-003-p02.txt");
    else if(period=="IC86_III") ic_times->SetTimesFromMJDFile("ListOfTimes/HugeListOfTimes_IC86-III_GRL_version-003-p02.txt");
    else if(period=="IC86_IV") ic_times->SetTimesFromMJDFile("ListOfTimes/HugeListOfTimes_IC86-IV_GRL_version-003-p02.txt");
    else if(period=="IC86_V") ic_times->SetTimesFromMJDFile("ListOfTimes/HugeListOfTimes_IC86-V_GRL_version-003-p02.txt");
    else if(period=="IC86_VI") ic_times->SetTimesFromMJDFile("ListOfTimes/HugeListOfTimes_IC86-VI_GRL_version-003-p02.txt");
    else if(period=="IC86_VII") ic_times->SetTimesFromMJDFile(TString("ListOfTimes/HugeListOfTimes_IC86-VII_GRL_version-003-p02.txt"));
    psData->SetEventTimeModulePtr(ic_times);
    psData->GenerateDataSet_with_nSrcEvents(0); // needs an I3SignalGenerator first??
  } 
  
  ark.psData = psData;
 

}
