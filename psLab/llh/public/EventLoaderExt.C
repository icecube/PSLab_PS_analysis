#include "EventLoaderExt.h"

// This has GetValueFromTree function
#include "rootExt/public/FunctionsRoot.h"
#include "rootExt/public/log_report.h"  // for log_warns, etc.
#include "rootExt/public/TreeConverter.h"
#include "rootExt/public/randomfunctions.h"
#include "llh/public/CoordinateTransform.h"

#include "TRandom.h"
#include "TF2.h"
#include "TH3D.h"
#include "TCanvas.h"

#include "llh/public/Kent.C"
#include "llh/public/BetaProf.C"

void EventLoaderExt::LoadSourcePdf_Kent(EquatorialDeg srcEqDeg, double sigmaSrc) {

  srcEqDeg_ = EquatorialDeg(srcEqDeg);

  // Implementation #3: Simplified Kent Distribution acts like Gaus2 over whole sky!
  //
  // Kent Distribution with no ellipticity (Beta=0).
  // This is a Gaussian on a sphere, angles in radians.
  // Parameters are: 0) Normalization (ignorable unless numerics go awry!)
  //                 1) Source Sigma (acts just like a Gaus Sigma!)
  //                 2) Mean Theta in Rad
  //                 3) Mean Phi in Rad
  // Beta=0; whole sky; no approx
  //TF2 *fSourcePdf = new TF2("Kent","[0]*exp( 1/([1]**2)*( sin(y)*sin([2])*cos(x-[3])+cos(y)*cos([2]) ) )",0,2*TMath::Pi(),0,TMath::Pi()); // Kent

  // Set range of function to have a reach of 3*sigma from the mean. 

  // This function handles the Kent distribution, using an approx when sigma 
  //  is small to avoid numerical problems (float limitations)
  //cerr << "Trying to load pdf from Kent...\n";
  // Input params in degrees
  hSourcePdf_ = (TH2D*)Kent(srcEqDeg.GetDec()+90., srcEqDeg.GetRa(), sigmaSrc);
  // Function will be drawn: TODO: Find a way to do this without drawing.
  // The function still works in batch mode though: "root -b" so not a problem 

}
/*
void EventLoaderExt::LoadSourcePdf_CG_CentralAGN(EquatorialDeg srcEqDeg, double z, double rCore, bool singleBeta, double n1, double rc1, double beta1, double n2, double rc2, double beta2) {
  //  Central AGN
  cout << "APPROX CENTRAL AGN MODEL AS POINT SOURCES!!!\n";
  cout << "This is strictly true up to the break energies ~1e17.5 eV...\n";
  cout << "Use LoadSourcePdf_Kent with srcSigma = 0.\n";
  exit(1);

}
*/

void EventLoaderExt::LoadSourcePdf_CG_ABI(char *model, EquatorialDeg srcEqDeg, double z, double rVirialMpc, bool singleBeta, double n1, double rc1, double beta1, double n2, double rc2, double beta2) {

  // "Model A", "Model B", and Isobaric Model can all be loaded with this function.

  TString modelSt(model);

  if (modelSt == "A") {
    // "Model A" (uniform CRs within rShock = 0.56*rVirial, 
    cout << "Loading CG \"Model A\" (uniform CR density within rShock)...\n";
  } else if (modelSt == "B") {
  // "Model B" (uniform CRs within r_virial)
    cout << "Loading CG \"Model B\" (uniform CR density within rVirial)...\n";
  } else if (modelSt == "I") {
  // Isobaric
  // Under this model, the CR distribution follows the electron distribution,
  //  so we just square the BetaProf of each cluster to get hSourcePdf_
    cout << "Loading CG \"Isobaric Model\" (CR density traces electron density)...\n";
  } else {
    cout << "Error: Must configure the CG model correctly:\n options are: \"A\", \"B\", or \"I\"\n";
    exit(1);
  }

  srcEqDeg_ = EquatorialDeg(srcEqDeg);

  // From z, calc dist and angular extent from Hubble's Law:
  const double c = 299792.458; //      km/s
  const double H = 70; // Hubble Const km/s/Mpc
  double distClusterMpc = z*c/H; // Mpc
  // Angular Extent = Radial Dist / distClusterMpc
  double rVirialRad = rVirialMpc/distClusterMpc; // virial radius in radians
  //cout << "rVirialMpc: " << rVirialMpc << endl;
  //cout << "rVirialRad: " << rVirialRad << endl;

  // Electron Density (arbitrary units) vs Radius (given in kpc)
  TF1 *fBetaProf = (TF1*)BetaProf(singleBeta, n1, rc1, beta1, n2, rc2, beta2, rVirialMpc);

  int nBins = fBetaProf->GetNpx();  // Using default: 100
  double aziClusterRad  = srcEqDeg.GetRa()*TMath::DegToRad();
  double zenClusterRad = (srcEqDeg.GetDec()+90)*TMath::DegToRad();
  //cout << "aziClusterRad: " << aziClusterRad<< endl;
  //cout << "zenClusterRad: " << zenClusterRad<< endl;

  double xminRad, xmaxRad, yminRad, ymaxRad;
  xminRad = aziClusterRad-rVirialRad/sin(zenClusterRad);
  xmaxRad = aziClusterRad+rVirialRad/sin(zenClusterRad);
  yminRad = zenClusterRad-rVirialRad;
  ymaxRad = zenClusterRad+rVirialRad;

  double zminMpc, zmaxMpc;
  zminMpc = -rVirialMpc;
  zmaxMpc = rVirialMpc;

  // Z-units are left as Mpc so we can integrate column depth
  TH3D *h3Dmodel = new TH3D("h3Dmodel","h3Dmodel",
    nBins, xminRad, xmaxRad,
    nBins, yminRad, ymaxRad,
    nBins, zminMpc, zmaxMpc);


  // To construct hSourcePdf from 3D cluster models:
  //  loop over histo, calc radius of this bin (Rad),
  //  eval fBetaProf @ this radius, integrate column depth (Mpc), fill 2D histo.
  // (This is how to convert the 1d beta-prof and 3D CR dist, to a 2d source pdf)
  double xRad=0, yRad=0, zMpc=0, rhoRad=0, rhoMpc=0, rMpc=0, value=0;
  for (int i=1; i<=nBins; i++) {
    for (int j=1; j<=nBins; j++) {
      for (int k=1; k<=nBins; k++) {
        xRad = h3Dmodel->GetXaxis()->GetBinCenter(i);
        yRad = h3Dmodel->GetYaxis()->GetBinCenter(j);
        zMpc = h3Dmodel->GetZaxis()->GetBinCenter(k);
        // polar dist ('r_perp') to center of CG in radians
        rhoRad = acos(sin(yRad)*sin(zenClusterRad)*cos(xRad-aziClusterRad) + 
                 cos(yRad)*cos(zenClusterRad));
        rhoMpc = rhoRad*distClusterMpc;
        rMpc = sqrt(rhoMpc*rhoMpc + z*z); // Mpc from cluster center
        if (modelSt == "A") {
          // CRs uniform inside rShock = 0.56*rVirial
          if (rMpc>0.56*rVirialMpc) {
            value = 0;
          } else {
            // Product of CR*electron density gives neutrinos
            value = fBetaProf->Eval(rMpc*1.e3); // takes kpc as input
          }
        } else if (modelSt == "B") {
          if (rMpc>rVirialMpc) {
            value = 0;
          } else {
            // CRs proportional to electron density for isobaric model
            value = fBetaProf->Eval(rMpc*1.e3); // takes kpc as input
          }
        } else if (modelSt == "I") {
            // CRs proportional to electron density for isobaric model
            // Product of CR*electron density gives neutrinos
          if (rMpc>rVirialMpc) {
            value = 0; // approx as 0 outside of virial radius
          } else {
            value = pow( fBetaProf->Eval(rMpc*1.e3), 2 );
          }
        } else {
          cout << "Error: Must configure the CG model correctly:\n";
          cout << " options are: \"A\", \"B\", or \"I\"\n";
          cout << "Double Error: This should have been caught by an earlier check!\n";
          cout << "You should never see this message!\n";
          exit(1);
        }
        h3Dmodel->SetBinContent(i,j,k,value);
      } // finish loop over z
    }   // finish loop over y
  }     // finish loop over x

  // We now have a histogram of the 3D model of the cluster (in Mpc).  We need
  //  to project it onto the x-y plane and convert to angular coords on sky.

  // hSourcePdf_ is in radians and represents the source distribution.
  // It is centered at the center of the cluster and extends +or- virial radius.
  // It's coords represent Theta and Phi of the source (like Kent), where
  //  theta = dec+90, phi = ra
  //TH2D *hSourcePdf = new TH2D("hSourcePdf","hSourcePdf",
  //  nBins, xminRad, xmaxRad,
  //  nBins, yminRad, ymaxRad);

  TH2D *hSourcePdf = (TH2D*)h3Dmodel->Project3D("yx")->Clone("hSourcePdf");

  // Normalize hSourcePdf (not strictly necessary)
  double sum = hSourcePdf->Integral();
  cout << sum << endl;
  hSourcePdf->Scale(1./sum);

  hSourcePdf_ = hSourcePdf; // Set member pointer to new source pdf

}


void EventLoaderExt::LoadSourceEvents(vector<I3Event>& srcCandidateEvents) {

  stopwatch_LoadSourceEvents_.Start(false);

  if (!hSourcePdf_) {
    cout << "Must load source pdf first (for flexibility)\n";
    cout << "E.g. LoadSourcePdf_Kent(EquatorialDeg srcEqDeg, double sigmaSrc)\n";
  }

  assert(hSourcePdf_);
  
  // Is there a way to assert source coord is set also?
  //assert(srcEqDeg_.GetRa());

  //nk new vector stuff
  if (rawSrcEventVect_.size() == 0) {
    InstallEvents(rawSrcEventVect_, sourceTree_, true);
  }

  srcCandidateEvents.clear();

  // time is arbitrary here, we just care about the whole zenith band
  double arbitrary_timeMJD = 54650.5;
  
  
  double bandSolidAngle;
  
  // Set up things related to grabbing only mc events near source zenith
    
  double srcZenDeg, srcAziDeg;
  EqToLocal(srcEqDeg_.GetRa(),srcEqDeg_.GetDec(), arbitrary_timeMJD,
	    srcZenDeg, srcAziDeg);
  
  // this assumes we are at pole, and only need one zenith band
  
  double srcZenRad = srcZenDeg*TMath::DegToRad();
  double zenWidthRad = zenWidthDeg_*TMath::DegToRad();
  /*if(zenWidthDegRangeVect_.size() != 0)
    {
      for(int i = 0; i < int(zenWidthDegRangeVect_.size()); i++)
	{
	  
	  if (srcZenDeg > zenWidthDegRangeVect_[i][1] && srcZenDeg < zenWidthDegRangeVect_[i][2])
	    {
	      
	      zenWidthRad = zenWidthDegRangeVect_[i][0] * TMath::DegToRad();
	      break;
	    }
	}
      
    }*/
  
  double srcMinRad = srcZenRad - zenWidthRad;
  double srcMaxRad = srcZenRad + zenWidthRad;
  bandSolidAngle = 2.*TMath::Pi()*( cos(srcMinRad)-cos(srcMaxRad) );
  gRandom->SetSeed(get_ran1_seed());
  
  // We select the different patches according to the source PDF.
  // Then we do the point-source simulation for each patch (keeping the band).
  
  for (int i = 0; i < nUpSamples_; i++)
    {
      double xRan = 0, yRan = 0, xRanDeg = 0, yRanDeg = 0;
      hSourcePdf_->GetRandom2(xRan, yRan); // Random Phi and Theta in radians
      xRanDeg = xRan*TMath::RadToDeg();
      yRanDeg = yRan*TMath::RadToDeg();
      double srcZenDeg_new, srcAziDeg_new;
      
      //New coordinates for a pseudo-point source.
      
      cout << "New source position: " << xRanDeg 
	   << " r.a.  " << yRanDeg - 90 << " declination " << endl;
      

      EqToLocal(xRanDeg, yRanDeg-90., // Convert Theta to Dec
		arbitrary_timeMJD,
		srcZenDeg_new, srcAziDeg_new);
      
      srcMinRad = srcZenDeg_new * TMath::DegToRad() - zenWidthRad;
      srcMaxRad = srcZenDeg_new * TMath::DegToRad() + zenWidthRad;
      //New MC local coordinates

      double mcZenRad_new = srcZenDeg_new * TMath::DegToRad();
      double mcAziRad_new = srcAziDeg_new * TMath::DegToRad();
      
      //Now let's loop over events the old fashion way
      for (int n = 0; n < int(rawSrcEventVect_.size()); ++n) {
	if (srcMCZenRadVect_[n] < srcMinRad || srcMCZenRadVect_[n] > srcMaxRad) {
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
	
	// these have been rotated from the simulation, as per above!!
	params.recoZenithDeg = TMath::RadToDeg() * ZenRad;
    
    
	params.recoAzimuthDeg = TMath::RadToDeg() * AziRad;
    
	// MAJOR ISSUE ESP. FOR TIME-DEPENDENT SIGNAL SIMULATION:
	// THIS OLD METHOD MOVES THE ORIGINAL AZIMUTH (BOTH TRUE AND RECO)
	// TO A COMPLETELY DIFFERENT LOCATION
	
	e.SetParams(params);
	
	mcparams.PS_FlatSpectrumRate = srcMCOneWeightVect_[n] /
	  (totalGeneratedEvents_ * bandSolidAngle) / 
	  nUpSamples_;
	  
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
    }
  
  
  if (monitor_) {
    cout << "Selected " << srcCandidateEvents.size() << " events using new ExtLoader with many dec bands.\n";
  }
  stopwatch_LoadSourceEvents_.Stop();
  stopwatch_LoadSourceEvents_.Print();
}

void EventLoaderExt::LoadSourceEvents(vector<I3Event>& srcCandidateEvents,
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
