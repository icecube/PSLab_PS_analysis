#include "llh/public/ZenithEnergyProb.h"

#include "TCut.h"
#include "TMath.h"  // for TMath::RadToDeg()

#include "rootExt/public/FunctionsRoot.h"  // for VerifyExpresion()
#include "rootExt/public/log_report.h"  // for log_warns, etc.


// CONSTRUCTOR

ZenithEnergyProb::ZenithEnergyProb() : 
  sourceZenWidthDeg_(0),
  recoZenRadName_("") 
{
  vector<double> tempVect;
  tempVect.push_back(0.);
  tempVect.push_back(180.);
  SetZenithBandsDeg(tempVect);
}


//   PUBLIC FUNCTIONS   //

//
// SET FUNCTIONS
//


void ZenithEnergyProb::SetZenithBandsDeg(const vector<double>& zenMinDegVect) {
  selectZenithBand_ = 1; // just reset it to safe value

  hZenDegBands_.Reset();
  if (zenMinDegVect.size() < 2) {
    log_fatal("Could not set zenith bands with less than two boundaries.\n");
  }
  hZenDegBands_.SetBins(zenMinDegVect.size()-1, &zenMinDegVect[0]);

  // TH1D::SetBins makes checks for whether bins are in order, etc., and
  // will revert to default if not.
  // We check that histogram is successful and spans full zenith range:
  if (hZenDegBands_.GetNbinsX() == int(zenMinDegVect.size()-1) && 
      hZenDegBands_.GetBinLowEdge(1) == 0. &&
      hZenDegBands_.GetBinLowEdge(hZenDegBands_.GetNbinsX()+1) == 180.)
  {
    CreateBands();
  } else {
    log_fatal("Could not set zenith bands correctly in SetZenithBandsDeg.\n"
	      "Check whether full range (0,180) of zenith is covered and\n"
	      "boundaries are in ascending order.\n");
  }
}

void ZenithEnergyProb::SetZenithBandsRad(const vector<double>& zenMinRadVect) {
  vector<double> tempVect;
  for (unsigned int i=0; i<zenMinRadVect.size(); ++i) {
    tempVect.push_back( zenMinRadVect[i] * TMath::RadToDeg() );
  }
  SetZenithBandsDeg(tempVect);
}

// Create vector of SimpleEnergyProbs, pass along the settings
void ZenithEnergyProb::CreateBands() {
  eProbVect_.clear();
  eProbVect_.resize( GetNZenithBands()+1 );
  // 0th bin corresponds to underflow bin of histogram, and is not used

  // These fns propagate the main settings to the individual eProbs
  SetConstrainSignal(optConstrainSignal_);
  SetLoadModeNew(optLoadModeNew_);
  if (rangeIsSet_) {
    SetEnergyGammaRangeAndBackFill(
	      GetEnergyBins(), GetEnergyMin(), GetEnergyMax(),
	      GetGammaBins(), GetGammaMin(), GetGammaMax(),
	      GetEnergyNBinStartBackFill() );
  }
}


void ZenithEnergyProb::SetConstrainSignal(bool opt) {
  SimpleEnergyProb::SetConstrainSignal(opt);
  for (int i=1; i<= GetNZenithBands(); ++i) {
    eProbVect_[i].SetConstrainSignal(opt);
  }
}

void ZenithEnergyProb::SetLoadModeNew(bool opt) {
  SimpleEnergyProb::SetLoadModeNew(opt);
  for (int i=1; i<= GetNZenithBands(); ++i) {
    eProbVect_[i].SetLoadModeNew(opt);
  }
}

void ZenithEnergyProb::SetEnergyGammaRangeAndBackFill(
	    int nBinsEnergy, double energyMin, double energyMax,
	    int nBinsGamma, double gammaMin, double gammaMax,
	    int energyNBinStartBackFill) {
  SimpleEnergyProb::SetEnergyGammaRangeAndBackFill(
	      nBinsEnergy, energyMin, energyMax,
	      nBinsGamma, gammaMin, gammaMax,
	      energyNBinStartBackFill);
  for (int i=1; i<= GetNZenithBands(); ++i) {
    eProbVect_[i].SetEnergyGammaRangeAndBackFill(
	              nBinsEnergy, energyMin, energyMax,
	              nBinsGamma, gammaMin, gammaMax,
	              energyNBinStartBackFill);
  }
}


void ZenithEnergyProb::SetTableGamma(TH2D* Aeff, vector< vector<TH1D*> > logEproxy){
  for (int zBand=1; zBand <= GetNZenithBands(); ++zBand) {
    double zenMinDeg = hZenDegBands_.GetBinLowEdge(zBand) - sourceZenWidthDeg_;
    double zenMaxDeg = hZenDegBands_.GetBinLowEdge(zBand+1)+sourceZenWidthDeg_;

    //Printf("zBand = %d, zenMinDeg = %f, zenMaxDeg = %f", zBand, zenMinDeg, zenMaxDeg);

    double decMinDeg = zenMinDeg - 90.0;
    double decMaxDeg = zenMaxDeg - 90.0;

    int hemi=0;
    if(decMaxDeg<-10) hemi=0;
    else if(decMaxDeg<10) hemi=1;
    else hemi=2;

    eProbVect_[zBand].SetTableGamma(Aeff, logEproxy[hemi], decMinDeg, decMaxDeg);
  }
}

void ZenithEnergyProb::SetTableGamma(TTree *srcTree,TCut cut,TString enString)
{
  if (recoZenRadName_ == "") {
    log_fatal("recoZenRadName is empty. Must be set in ZenithEnergyProb\n");
  }

  for (int zBand=1; zBand <= GetNZenithBands(); ++zBand) {
    double zenMinDeg= hZenDegBands_.GetBinLowEdge(zBand) - sourceZenWidthDeg_;
    double zenMaxDeg= hZenDegBands_.GetBinLowEdge(zBand+1)+sourceZenWidthDeg_;

    TString zenString = recoZenRadName_ + " >= ";
    zenString += ( zenMinDeg * TMath::DegToRad() );
    zenString += " && " + recoZenRadName_ + " < ";
    zenString += ( zenMaxDeg * TMath::DegToRad() );
    
    TCut zenCut(zenString);
    assert( VerifyExpression(srcTree, zenCut.GetTitle()) );

    eProbVect_[zBand].SetTableGamma(srcTree, zenCut*cut, enString);
  }

  // For now, we deliberately do *not* set the corresponding Bkg Table 
  // in the ZenithEnergyProb object itself.  (So any accidental reference
  // to it will be detected since it is unfilled)
}


void ZenithEnergyProb::SetTableBkg(const vector<I3Event>& eventVect)
{
  // Make separate EventPtrList for each zenith band.
  // Range is 1 through NBands
  vector<EventPtrList> selectedVectArray;
  selectedVectArray.resize( GetNZenithBands() +1 ); 

  // Separate events into different vectors
  for (unsigned int j=0; j<eventVect.size(); ++j) {
    int zBand = GetZenDegBand(eventVect[j]);
    selectedVectArray[zBand].AddEvent(&(eventVect[j]));
    // don't copy whole event, just point to it
  }

  // Fill each band
  for (int zBand=1; zBand <= GetNZenithBands() ; ++zBand) {
    eProbVect_[zBand].SetTableBkg( selectedVectArray[zBand] );
  }
  // For now, we deliberately do *not* set the corresponding Bkg Table 
  // in the ZenithEnergyProb object itself.  (So any accidental reference
  // to it will be detected since it is unfilled)
}  //
//
// TWO-VERSIONS RIGHT NOW
//

void ZenithEnergyProb::SetTableBkg(const EventPtrList& evList)
{
  // Make separate EventPtrList for each zenith band.
  // Range is 1 through NBands
  vector<EventPtrList> selectedVectArray;
  selectedVectArray.resize( GetNZenithBands() +1 ); 

  // Separate events into different vectors
  for (int j=0; j<evList.GetSize(); ++j) {
    int zBand = GetZenDegBand(*(evList.GetEvent(j)));
    selectedVectArray[zBand].AddEvent(evList.GetEvent(j));
    // don't copy whole event, just point to it
  }

  // Fill each band
  for (int zBand=1; zBand <= GetNZenithBands() ; ++zBand) {
    eProbVect_[zBand].SetTableBkg( selectedVectArray[zBand] );
  }
  // For now, we deliberately do *not* set the corresponding Bkg Table 
  // in the ZenithEnergyProb object itself.  (So any accidental reference
  // to it will be detected since it is unfilled)
}




void ZenithEnergyProb::SelectZenithBand(int band) {
  selectZenithBand_ = ValidateZenithBand(band);
}


//
// GET FUNCTIONS()   CONST
//

int ZenithEnergyProb::GetZenDegBand(double zenDeg) const {
  int band = hZenDegBands_.GetXaxis()->FindBin(zenDeg);
  if (band==0 || band == GetNZenithBands()+1) {
    band = -1;
    log_error("GetZenDegBand(%lg) went out of bounds. Setting band= -1\n",
	      zenDeg);
  }
  return band;
}

int ZenithEnergyProb::GetZenDegBand(const Event& event) const {
  const I3Event& i3event = dynamic_cast<const I3Event&>(event);
  double zenDeg = i3event.GetParams().recoZenithDeg;
  return GetZenDegBand(zenDeg);
}

int ZenithEnergyProb::ValidateZenithBand(int band) const {
  if (band>=1 && band<=GetNZenithBands() ) {
    return band;
  } else {
    log_error("ValidateZenithBand(%d) out of range.  Set to 1.\n",band);
    return 1;
  }
}


double ZenithEnergyProb::GetEnergyProbGamma
  (const Event& event, double gamma) const
{
  int band = GetZenDegBand(event);
  return eProbVect_[band].GetEnergyProbGamma(event, gamma);
}

double ZenithEnergyProb::GetEnergyProbBkg(const Event& event) const
{
  int band = GetZenDegBand(event);
  return eProbVect_[band].GetEnergyProbBkg(event);
}

double ZenithEnergyProb::GetEnergyMaxRatio(const Event& event) const
{
  int band = GetZenDegBand(event);
  return eProbVect_[band].GetEnergyMaxRatio(event);
}

