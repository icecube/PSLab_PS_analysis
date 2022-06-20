#include "llh/public/EventLoader.h"

#include "rootExt/public/FunctionsRoot.h" // for GetValueFromTree function
#include "rootExt/public/log_report.h"  // for log_warns, etc.
#include "rootExt/public/randomfunctions.h"
#include "rootExt/public/TreeConverter.h"

#include "llh/public/CoordinateTransform.h"


void EventLoader::InstallEvents(vector<I3Event>& eventVect, 
				TTree* tree, bool installSrc) {
  eventVect.clear();

  cout << endl;
  cout << "New data structure " << GetDataStructure() << endl;
  cout << endl;
  
  cout << endl;
  cout << "Apply ang floor " << GetApplyAngFloor() << endl;
  cout << endl;
  
  if (installSrc) {
    srcMCZenRadVect_.clear();
    srcMCAziRadVect_.clear();
    srcMCOneWeightVect_.clear();
  }
  else {   // installing base data
    bkgWeightVect_.clear();

    if ( bkgTimeMethod_ == SCRAMBLE ) {
      // check that we actually have an EventTimeModule to point to
      if ( !evTimeModulePtr_ ) {
	evTimeModulePtr_ = new EventTimeModule();
	// if not, make a default one
      }
      // (very minor memory leak here if never deleted, but it only happens
      //  once so not very worrisome)
      
      // times from generic module with default constructor 
      // are only generated over one sidereal day... this is *supposed* to look
      // wierd if someone is expecting a realistic distribution of times
      // during the year...
      
      // This maintains backward compatibility for time-independent scripts...
      // eventually we may think to remove this and force _all_ scripts to
      // specify how time is handled.

      // Clear list (if any) of used times which the module tracks internally
      evTimeModulePtr_->ResetUsedTimes();
    }
  }

  TreeConverter tc;
  tc.SetMonitor(monitor_);

  for (unsigned int i=0; i<cutStringVect_.size(); ++i) {
    tc.AddCut( cutStringVect_[i] );
  }

  tc.AddVar("ZenRad",name_recoZenith_rad_);
  tc.AddVar("AziRad",name_recoAzimuth_rad_);
  tc.AddVar("SigmaDeg",name_sigmaDeg_);
  tc.AddVar("energyValue",name_energyValue_);
  tc.AddVar("runID",name_runID_);
  tc.AddVar("eventID",name_eventID_);

  if (installSrc) {
    
    tc.AddVar("mcEnergyGeV",name_mcEnergy_GeV_);
    
    if (newDataStructure_){
        cout << "Setting vars for new data structure " << endl;
        tc.AddVar("timeMJD",name_timeMJD_);
        tc.AddVar("mcRaRad",name_mcRa_rad_);
        tc.AddVar("mcDecRad",name_mcDec_rad_);
    }
    
    else {
        cout << "Setting vars for old data structure " << endl;
        tc.AddVar("timeMJD","0.");  // ignore for now, will be set for source
        tc.AddVar("mcZenRad",name_mcZenith_rad_);
        tc.AddVar("mcAziRad",name_mcAzimuth_rad_);
    }
    
    tc.AddVar("mcOneWeight",name_mcOneWeight_);
  } else {
    tc.AddVar("timeMJD",name_timeMJD_);
    tc.AddVar("bkgWeight",bkgWeight_);
  }

  tc.SetTree(tree);

  const map<TString,double>* resultMap;
  
  // this action automatically skips to next entry which passes cuts
  while( tc.GetNextEntry() ) {

    resultMap = tc.GetResultMap();

    I3Event e;
    I3EventParameters params;

    double ZenRad = resultMap->find("ZenRad")->second;
    double AziRad = resultMap->find("AziRad")->second;

    params.recoZenithDeg = TMath::RadToDeg() * ZenRad;
    params.recoAzimuthDeg = TMath::RadToDeg() * AziRad;
    params.parafitSigmaDeg = resultMap->find("SigmaDeg")->second;
    //temp. fix to avoid too small sigma -> 0.2deg floor
    if(applyAngFloor_ && params.parafitSigmaDeg<0.2) {
        
        //cout << "Applying angular floor" << endl;
        params.parafitSigmaDeg = 0.2;
    }
    params.energyValue = resultMap->find("energyValue")->second;
    params.runID = int(resultMap->find("runID")->second);
    params.eventID = int(resultMap->find("eventID")->second);

    e.SetParams(params);

    double timeMJD = resultMap->find("timeMJD")->second;
    e.SetMJD(timeMJD);

    // NOTE THAT Equatorial Coordinates are not set during Install...
    //    e.SetCoord(EquatorialDeg(raDeg,decDeg));

    if (installSrc) {
      I3MCParameters mcparams;

      mcparams.mcEnergy = resultMap->find("mcEnergyGeV")->second;

      // this gets set later, once the method of source generation is specified
      mcparams.PS_FlatSpectrumRate = 0.;

      // this gets set later during analysis, when spectralIndex is chosen
      mcparams.srcWeight = 0.;   

      e.SetMCParams(mcparams);

      // Other mc information
      if (newDataStructure_){
          double mcRaDeg     = (resultMap->find("mcRaRad")->second)*TMath::RadToDeg();
          double mcDecDeg    = (resultMap->find("mcDecRad")->second)*TMath::RadToDeg();
          double mcZenDeg=0, mcAziDeg=0;
          EqToLocal(mcRaDeg,mcDecDeg,timeMJD,mcZenDeg, mcAziDeg);
          double mcZenRad=mcZenDeg*TMath::DegToRad();
          double mcAziRad=mcAziDeg*TMath::DegToRad();
          srcMCZenRadVect_.push_back(mcZenRad);
          srcMCAziRadVect_.push_back(mcAziRad);
      }
      
      else {
          double mcZenRad    = resultMap->find("mcZenRad")->second;
          double mcAziRad    = resultMap->find("mcAziRad")->second;
          srcMCZenRadVect_.push_back(mcZenRad);
          srcMCAziRadVect_.push_back(mcAziRad);
      }
      
      double mcOneWeight = resultMap->find("mcOneWeight")->second;
      
      
      srcMCOneWeightVect_.push_back(mcOneWeight);
    }
    else {
      double bkgWeight = resultMap->find("bkgWeight")->second;
      bkgWeightVect_.push_back(bkgWeight);
    }

    eventVect.push_back(e);
  }

  if (monitor_) {
    cout << "(" << eventVect.size() << " entries in tree passed cuts)\n";
  }

}


void EventLoader::LoadBkgEvents(vector<I3Event>& bkgEvents) {
  if (rawBkgEventVect_.size() == 0) {
    cout << "Installing Events " << endl;

    InstallEvents(rawBkgEventVect_, bkgTree_, false);
  }

  bkgEvents.clear();

  // Sanity check if we are going to over-sample:
  if (bkgLoadMethod_ != EXACT ) {
    if (minZenithRad_ >= maxZenithRad_) {
      log_fatal("\nEventLoader: Bad range for spreading bkg events:  "
		"minZenithDeg >= maxZenithDeg\n");
      exit(1);
    }
    if (spreadZenithRad_ <= 0.) {
      log_fatal("\nEventLoader: Bad range for spreading bkg events:  "
		"spreadZenithDeg <= 0.\n");
      exit(1);
    }
  }

  // For Monitoring Purposes
  unsigned int nMultiSampled[5] = {0};
  int maxSampled = 0;
  double sumWeights = 0.;

  for (int n=0; n<int(rawBkgEventVect_.size()); ++n) {

    // Figure out whether to add this event 0 times, once, or more...

    int nLoadEvents = 0;   // if some undefined method shows up, load nothing

    if (bkgLoadMethod_ == EXACT) {
      nLoadEvents = 1;
    } 
    else if (bkgLoadMethod_ == SELECTIVE_SAMPLE) 
    {
      double thisWeight = bkgWeightVect_[n];
      sumWeights += thisWeight;
      double r = random_uniform(0.,1.);
      nLoadEvents = ( r < thisWeight );  // either 0 or 1
    } 
    else if (bkgLoadMethod_ == POISSON_SAMPLE_PLUS)
    {
      double thisWeight = bkgWeightVect_[n];
      sumWeights += thisWeight;
      nLoadEvents = random_poisson(thisWeight);  // can be 0,1,2,...
    }
    else if (bkgLoadMethod_ == FIXED_UPSAMPLE)
    {
      double thisWeight = bkgWeightVect_[n];
      sumWeights += thisWeight;
      nLoadEvents = int(thisWeight);  // can be 0,1,2,...
    }

    
    // For Monitoring Purposes
    if (nLoadEvents<4) {
      ++nMultiSampled[nLoadEvents];
    } else { 
      ++nMultiSampled[4]; 
    }
    if (nLoadEvents>maxSampled) {
      maxSampled = nLoadEvents;
    }


    for (int i=0; i<nLoadEvents; ++i) {

      // copy the nth raw event (i.e. the current event in rawBkgEventVect_)
      I3Event e = rawBkgEventVect_[n];

      // PSUEDO-EVENTS for over-sampling:
      // If we're sampling this event more than once,
      // Modify ZenRad to create a semi-random new zenith
      // (but: near original, still passing zen cut, and still same params)
      if (i>0) 
      {
	// we will be modifying the zenith param:
	I3EventParameters params = e.GetParams();

	double ZenRad = params.recoZenithDeg  * TMath::DegToRad();

	double minBoundRad = ZenRad - spreadZenithRad_;
	if (minBoundRad < minZenithRad_) {
	  minBoundRad = minZenithRad_;
	}
	double maxBoundRad = ZenRad + spreadZenithRad_;
	if (maxBoundRad > maxZenithRad_) {
	  maxBoundRad = maxZenithRad_;
	}
	// assign a new random zenith, within range, according to phase space
	double cosNewZen = 
	  random_uniform(cos(maxBoundRad), cos(minBoundRad));

	ZenRad = acos(cosNewZen);

	// Modify event here
	params.recoZenithDeg = TMath::RadToDeg() * ZenRad;
	e.SetParams(params);
      }

      // Scramble event time, if requested
      // If randomizeOnlyTime_ is true, true RA and dec are assigned
      // to the events, and then the time is scrambled.
      // If it is false, time is scrambled and then RA is corrected (hence scrambled itself)
      if(randomizeOnlyTime_){
        double raDeg=0, decDeg=0;
        LocalToEq(e.GetParams().recoZenithDeg, e.GetParams().recoAzimuthDeg,
    e.GetTime().GetMJD(),
    raDeg, decDeg);
        e.SetCoord(EquatorialDeg(raDeg,decDeg));

      }

      if (bkgTimeMethod_ == SCRAMBLE) {
	e.SetTime( evTimeModulePtr_->GetRandomTime() );
      }

      if(!randomizeOnlyTime_){
        double raDeg=0, decDeg=0;
        LocalToEq(e.GetParams().recoZenithDeg, e.GetParams().recoAzimuthDeg,
                e.GetTime().GetMJD(),
                raDeg, decDeg);
        e.SetCoord(EquatorialDeg(raDeg,decDeg));
      }

      bkgEvents.push_back(e);
    }

  }


  // Give some output if user is monitoring 

  if (monitor_) {
    cout << "Background set has " << bkgEvents.size() << " events.\n";

    if (bkgLoadMethod_ == SELECTIVE_SAMPLE) {
      cout << "SelectiveSample based on weights: \n"; 
      cout << "  " << nMultiSampled[0] << " events skipped\n";
      cout << "  " << nMultiSampled[1] << " events loaded once.\n";
      cout << "Number of events expected from sum of Weights: ";
      cout << sumWeights << endl;
    }

    if (bkgLoadMethod_ == POISSON_SAMPLE_PLUS ||
	bkgLoadMethod_ == FIXED_UPSAMPLE) {
      cout << "PoissonSamplePlus or FixedUpsample based on weights: \n"; 
      cout << "  " << nMultiSampled[0] << " events skipped.\n";
      cout << "  " << nMultiSampled[1] << " events loaded once.\n";
      cout << "  " << nMultiSampled[2] << " events sampled twice.\n";
      cout << "  " << nMultiSampled[3] << " events sampled thrice.\n";
      cout << "  " << nMultiSampled[4] << " events sampled >=4  times.\n";
      if (maxSampled>=4) {
	cout << "  " << maxSampled;
	cout << " is the maximum times a single event was sampled.\n";
      }
      cout << "Number of events expected from sum of Weights: ";
      cout << sumWeights << endl;
    }
  }

}



void EventLoader::SetSourceTree(TTree *tree) {
    
  if (newDataStructure_) {
      name_mcRa_rad_  = "mcPrimary_Ra_rad";
      name_mcDec_rad_ = "mcPrimary_Dec_rad";
  }
  
  else {
      name_mcZenith_rad_  = "mcPrimary_Zenith_rad";
      name_mcAzimuth_rad_ = "mcPrimary_Azimuth_rad";
  }
  
  name_mcEnergy_GeV_  = "mcPrimary_Energy_GeV";
  name_mcOneWeight_   = "mcOneWeight";
  TString name_TotGenEvents = "mcTotalGeneratedEvents";

  // Stop if any of these names cannot be evaluated for the tree
  
  if (newDataStructure_) {
      assert( VerifyExpression(tree, name_mcRa_rad_) );
      assert( VerifyExpression(tree, name_mcDec_rad_) );
  }
  else {
    assert( VerifyExpression(tree, name_mcZenith_rad_) );
    assert( VerifyExpression(tree, name_mcAzimuth_rad_) );
  }
  assert( VerifyExpression(tree, name_mcEnergy_GeV_) );
  assert( VerifyExpression(tree, name_mcOneWeight_) );
  assert( VerifyExpression(tree, name_TotGenEvents) );
  
  totalGeneratedEvents_ = GetValueFromTree(tree, name_TotGenEvents);

  sourceTree_ = tree;
}


void EventLoader::LoadSourceEvents(vector<I3Event>& srcCandidateEvents,
				   EquatorialDeg srcEqDeg, bool rotate) {
  if (rawSrcEventVect_.size() == 0) {
    InstallEvents(rawSrcEventVect_, sourceTree_, true);
  }

  srcCandidateEvents.clear();

  // time is arbitrary here, we just care about the whole zenith band
  double arbitrary_timeMJD = 51544.5;


  double bandSolidAngle=0;
  double srcMinRad=0, srcMaxRad=0;

  // Set up things related to grabbing only mc events near source zenith
  {
    double srcZenDeg=0, srcAziDeg=0;

    EqToLocal(srcEqDeg.GetRa(),srcEqDeg.GetDec(), arbitrary_timeMJD,
	      srcZenDeg, srcAziDeg);

    // this assumes we are at pole, and only need one zenith band

    double srcZenRad = srcZenDeg*TMath::DegToRad();
    double zenWidthRad = zenWidthDeg_*TMath::DegToRad();

    srcMinRad = srcZenRad - zenWidthRad;
    srcMaxRad = srcZenRad + zenWidthRad;
    // restrict within bounds, since we can't select events outside this range
    if (srcMinRad<0.) { srcMinRad = 0.; }
    if (srcMaxRad>TMath::Pi()) { srcMaxRad = TMath::Pi(); }

    bandSolidAngle = 2.*TMath::Pi()*( cos(srcMinRad)-cos(srcMaxRad) );

    if (monitor_) {
      cout << "Src at Zenith = " << srcZenRad*TMath::RadToDeg();
      cout << " ,  selecting events with MCZenith in range ( ";
      cout << srcMinRad*TMath::RadToDeg() << " , ";
      cout << srcMaxRad*TMath::RadToDeg() << " )\n";
    }
  }


  double mcZenDeg_new=0, mcAziDeg_new=0;

  // for this timeMJD, find out where the local mc dir. should have been
  EqToLocal(srcEqDeg.GetRa(), srcEqDeg.GetDec(), 
	    arbitrary_timeMJD,
	    mcZenDeg_new, mcAziDeg_new);
  double mcZenRad_new = mcZenDeg_new*TMath::DegToRad();
  double mcAziRad_new = mcAziDeg_new*TMath::DegToRad();


  for (int n=0; n<int(rawSrcEventVect_.size()); ++n) {

    if (srcMCZenRadVect_[n]<srcMinRad || srcMCZenRadVect_[n]>srcMaxRad) {
      continue;
    }

    I3Event e = rawSrcEventVect_[n];
    I3EventParameters params = e.GetParams();
    I3MCParameters mcparams = e.GetMCParams();

    double ZenRad = params.recoZenithDeg  * TMath::DegToRad();
    double AziRad = params.recoAzimuthDeg  * TMath::DegToRad();

    double ZenRad1 = ZenRad;
    double AziRad1 = AziRad;
    
    // Now Shift for MC source events
    if (rotate){
      double mcZenRad = srcMCZenRadVect_[n];
      double mcAziRad = srcMCAziRadVect_[n];

      TVector3 v;
      v.SetMagThetaPhi(1.,ZenRad, AziRad);
      v.RotateZ(-mcAziRad);                 // mc true now in x-z plane
      v.RotateY(mcZenRad_new - mcZenRad);   // mc true now at point src zenith
      v.RotateZ(mcAziRad_new);              // mc true now at point src azimuth
    
      // now our reconstructed track dir should have the same deviation
      // w.r.t. the point source location as it did to the original mc true

      ZenRad = v.Theta();
      AziRad = v.Phi();
    }

    // these have been rotated from the simulation, as per above!!
    params.recoZenithDeg = TMath::RadToDeg() * ZenRad;
    params.recoAzimuthDeg = TMath::RadToDeg() * AziRad;

    // MAJOR ISSUE ESP. FOR TIME-DEPENDENT SIGNAL SIMULATION:
    // THIS OLD METHOD MOVES THE ORIGINAL AZIMUTH (BOTH TRUE AND RECO)
    // TO A COMPLETELY DIFFERENT LOCATION

    e.SetParams(params);

    mcparams.PS_FlatSpectrumRate = srcMCOneWeightVect_[n] /
      (totalGeneratedEvents_ * bandSolidAngle);
    // In the usual formula for OneWeight (IceSim 2.0), we have:
    // "MyFlux * (OneWeight/NGeneratedEvents) * LiveTime"
    // Here, MyFlux and LiveTime will be specified for the source later,
    // So we just need the middle term (OneWeight/NGeneratedEvents)
    //
    // The extra factor 'bandSolidAngle':
    // OneWeight was calculated, assuming it would be multipled by 
    // a diffuse flux per steradian.
    // Here, we are selecting MC events within the zenith band specified,
    // and concentrating them at point, creating a point flux with no
    // solid angle term.  

    e.SetMCParams(mcparams);

    // Fill Coordinate Values

    double raDeg=0, decDeg=0;
    // (the important thing is that both mc and reco are
    // transformed to equatorial using the same timeMJD, even if arbitrary)
    LocalToEq(params.recoZenithDeg, params.recoAzimuthDeg, arbitrary_timeMJD,
	      raDeg, decDeg);

    //if (rotate) { e.SetCoord(EquatorialDeg(raDeg,decDeg)); }
    //else { e.SetCoord(EquatorialDeg(TMath::RadToDeg() * AziRad,TMath::RadToDeg() * ZenRad)); }
    e.SetCoord(EquatorialDeg(raDeg,decDeg));

    e.SetMJD(arbitrary_timeMJD);
    
    params.recoZenithDeg  = TMath::RadToDeg() * ZenRad1;
    params.recoAzimuthDeg = TMath::RadToDeg() * AziRad1;
    e.SetParams(params);

    srcCandidateEvents.push_back(e);

  }

  if (monitor_) {
    cout << "Selected " << srcCandidateEvents.size() << " events.\n";
  }
}
