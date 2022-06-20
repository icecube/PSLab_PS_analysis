#include "llh/public/SimpleEnergyProb.h"

#include "TAxis.h"
#include "TCut.h"
#include "TEventList.h"
#include "TH1.h" // for HistFillIn function

#include "rootExt/public/FunctionsRoot.h" // for VerifyExpresion()
#include "rootExt/public/log_report.h"  // for log_warns, etc.
#include "rootExt/public/TreeReader.h"
#include "rootExt/public/MakeNewName.h"
#include "fluxus/public/FluxFunction.h"



// CONSTRUCTOR

SimpleEnergyProb::SimpleEnergyProb() :
  optConstrainSignal_(true),
  optLoadModeNew_(true),
  rangeIsSet_(false),
  energyNBinStartBackFill_(0)
{ }


//   PRIVATE FUNCTIONS   //


int SimpleEnergyProb::GetEnergyBin(double energy) const {
  if (rangeIsSet_) {
    int eBin = hProbBkg_.GetXaxis()->FindBin(energy);
    // if energy is in underflow or overflow region, move it back w/in range
    if (eBin==0) { eBin++; }
    if (eBin==GetEnergyBins()+1) { eBin--; }
    return eBin;
  }
  log_error("Error: no return value, energy range is not defined.\n");
  return -1;
}


//   PUBLIC FUNCTIONS   //


//
// SET FUNCTIONS
//

void SimpleEnergyProb::SetEnergyGammaRangeAndBackFill(
	 int nBinsEnergy, double energyMin, double energyMax,
	 int nBinsGamma, double gammaMin, double gammaMax,
	 int energyNBinStartBackFill) {
  // Reset empties any previous contents.  GetSum() now equals zero.
  hProbBkg_.Reset();
  hProbBkg_.SetBins(nBinsEnergy, energyMin, energyMax);

  hProbGamma_.Reset();
  hProbGamma_.SetBins(nBinsEnergy, energyMin, energyMax,
		      nBinsGamma, gammaMin, gammaMax);

  hProbGammaMax_.Reset();
  hProbGammaMax_.SetBins(nBinsEnergy, energyMin, energyMax);

  energyNBinStartBackFill_ = energyNBinStartBackFill;
  rangeIsSet_ = true;
}

TH1D* SimpleEnergyProb::CalculateEnergyPDF(TH2D* effectiveArea, double gamma, double decMinDeg, double decMaxDeg){
  TH2D Aeff;
  effectiveArea->Copy(Aeff);
  Int_t binymin = Aeff.GetYaxis()->FindBin( TMath::Sin( TMath::DegToRad()*(decMinDeg+0.1) ) );
  Int_t binymax = Aeff.GetYaxis()->FindBin( TMath::Sin( TMath::DegToRad()*(decMaxDeg-0.1) ) );
  int Nxbins = Aeff.GetNbinsX();
  PowerLawFlux pflux(1,-gamma);

  double newSinDecBin = Aeff.GetYaxis()->GetBinLowEdge(binymax+1)-Aeff.GetYaxis()->GetBinLowEdge(binymin);
  double effAreaBin, oldSinDecBin, energy, BinWidth;
  for(int binx=1; binx<Nxbins+1; binx++){
    for(int biny=binymin; biny<binymax+1; biny++){
      effAreaBin = Aeff.GetBinContent(binx, biny);
      oldSinDecBin = Aeff.GetYaxis()->GetBinLowEdge(biny+1)-Aeff.GetYaxis()->GetBinLowEdge(biny);
      energy = pow(10, Aeff.GetXaxis()->GetBinCenter(binx) );
      BinWidth = pow(10, Aeff.GetXaxis()->GetBinLowEdge(binx+1)) - pow(10, Aeff.GetXaxis()->GetBinLowEdge(binx));
      Aeff.SetBinContent(binx, biny, effAreaBin * oldSinDecBin/newSinDecBin * BinWidth * pflux.GetFlux(energy));
    }
  }

  return Aeff.ProjectionX("EffectiveArea", binymin, binymax);
}

void SimpleEnergyProb::SetTableGamma(TH2D *effectiveArea, vector<TH1D*> logEproxy, double decMinDeg, double decMaxDeg) {

  if (!rangeIsSet_) {
    log_error("Error: Range must be set before setting hProbGamma.\n");
    return;
  }

  hProbGamma_.Reset();        // remove contents
  hProbGammaMax_.Reset();

  /*TH1D* Aeff = CalculateEffectiveArea(effectiveArea, decMinDeg, decMaxDeg);
  assert( Aeff->GetNbinsX()==GetEnergyBins() );

  for(int eBin=1; eBin<=GetEnergyBins(); eBin++){
    double effAreaValue = Aeff->GetBinContent(eBin);
    double energy = pow(10, Aeff->GetXaxis()->GetBinCenter(eBin) );
    for (int gBin=1; gBin<=GetGammaBins(); gBin++) {
      double gamma = hProbGamma_.GetYaxis()->GetBinCenter(gBin);
      //Printf("energy = %.1e, gamma=%.1f", energy, gamma);
      PowerLawFlux pflux(1,-gamma);
      hProbGamma_.SetBinContent(eBin, gBin, effAreaValue*pflux.GetFlux(energy));
    }
  } 
  */
  for (int gBin=1; gBin<=GetGammaBins(); gBin++) {
    double gamma = hProbGamma_.GetYaxis()->GetBinCenter(gBin);
    PowerLawFlux pflux(1,-gamma);
    TH1D* Aeff = CalculateEnergyPDF(effectiveArea, gamma, decMinDeg, decMaxDeg);
    double norm = Aeff->Integral();
    assert(norm>0);
    for(int eBin=1; eBin<=GetEnergyBins(); eBin++){
      double energy   = pow(10, Aeff->GetXaxis()->GetBinCenter(eBin) );
      double PDFvalue = Aeff->GetBinContent(eBin);
      int energybin = int(Aeff->GetBinCenter(eBin)/0.5)-4;
      for(int i=0; i<int(1e5*PDFvalue/norm); ++i){
        double logMuEnergy = logEproxy[energybin]->GetRandom();
        hProbGamma_.Fill(logMuEnergy, gamma);
      }
    }
  }

  // Normalize each gamma pdf separately
  assert( hProbGamma_.GetSum() > 0 );
  for (int gBin=1; gBin<=GetGammaBins(); gBin++) {
    double sum = 0.;
    // To get complete sum, here we include the underflow (0) and 
    // overflow (n+1) bins, because they can be non-empty if the user
    // set optConstrainSignal = false.  
    // Otherwise if opt = true, they will be zero and have no effect
    for (int eBin=1; eBin<=GetEnergyBins(); eBin++) sum += hProbGamma_.GetBinContent(eBin, gBin);
    for (int eBin=1; eBin<=GetEnergyBins(); eBin++) {
      double content = hProbGamma_.GetBinContent(eBin, gBin);
      hProbGamma_.SetBinContent(eBin, gBin, content / sum);
    }
  }

  // Find maximum for each energy
  for (int eBin=1; eBin <= GetEnergyBins(); ++eBin) {
    double max = 0.;
    for (int gBin=1; gBin<= GetGammaBins(); ++gBin) {
      double value = hProbGamma_.GetBinContent(eBin,gBin);
      if (value > max) { max = value; }
    }
    hProbGammaMax_.SetBinContent(eBin,max);
  }
}

void SimpleEnergyProb::SetTableGamma(TTree *srcTree, TCut cut,
				       TString enString) {
  if (!rangeIsSet_) {
    log_error("Error: Range must be set before setting hProbGamma.\n");
    return;
  }

  hProbGamma_.Reset();        // remove contents
  hProbGammaMax_.Reset();     // remove contents

 if (optLoadModeNew_) {
  assert( VerifyExpression(srcTree, cut.GetTitle()) );
  assert( VerifyExpression(srcTree, "mcPrimary_Energy_GeV") );
  assert( VerifyExpression(srcTree, "mcOneWeight") );
  assert( VerifyExpression(srcTree, enString) );

  TString listName = MakeNewName("srcEvList");
  TEventList eList(listName);
  srcTree->Draw(">> "+listName,cut);

  TreeReader tr_mcEnergy(srcTree, "mcPrimary_Energy_GeV");
  TreeReader tr_mcOneWeight(srcTree, "mcOneWeight");
  TreeReader tr_energy(srcTree, enString);

  // SafeMode==False is much faster, but see TreeReader header file
  // for example of when this can go wrong.  If in doubt, set to True.
  tr_mcEnergy.SetSafeMode(false);
  tr_mcOneWeight.SetSafeMode(false);
  tr_energy.SetSafeMode(false);

  for (int i=0; i<eList.GetN(); ++i) {

    double mcPrimary_Energy_GeV = tr_mcEnergy.GetValue( eList.GetEntry(i) );
    double mcOneWeight = tr_mcOneWeight.GetValue( eList.GetEntry(i) );

    double energyToFill;
    {
      double energy = tr_energy.GetValue( eList.GetEntry(i) );
      if (optConstrainSignal_) {
	// this corrects for out-of-bound energy values by using GetEnergyBin()
	int eBin = GetEnergyBin(energy);
	energyToFill = hProbBkg_.GetBinCenter(eBin);
      } else {
	energyToFill = energy;
      }
    }

    for (int gBin=1; gBin<=GetGammaBins(); gBin++) {
      double gamma = hProbGamma_.GetYaxis()->GetBinCenter(gBin);
      double weight = pow(mcPrimary_Energy_GeV,-fabs(gamma)) * mcOneWeight;
      hProbGamma_.Fill( energyToFill, gamma, weight );
    }
  }

  // Normalize each gamma pdf separately
  if (hProbGamma_.GetSum() > 0) {   // if == 0,  should we worry? 
    for (int gBin=1; gBin<=GetGammaBins(); gBin++) {
      double sum = 0.;
      // To get complete sum, here we include the underflow (0) and 
      // overflow (n+1) bins, because they can be non-empty if the user
      // set optConstrainSignal = false.  
      // Otherwise if opt = true, they will be zero and have no effect
      for (int eBin=0; eBin<=GetEnergyBins()+1; eBin++) {
	sum += hProbGamma_.GetBinContent(eBin, gBin);
      }
      for (int eBin=1; eBin<=GetEnergyBins(); eBin++) {
	double content = hProbGamma_.GetBinContent(eBin, gBin);
	hProbGamma_.SetBinContent(eBin, gBin, content / sum);
      }
    }
  }
 } 
 else {
  TString hName = MakeNewName("h1Energy");
  TH1D h1Energy(hName,"",GetEnergyBins(), GetEnergyMin(), GetEnergyMax() );

  for (int gBin=1; gBin<=GetGammaBins(); gBin++) {
    double gamma = hProbGamma_.GetYaxis()->GetBinCenter(gBin);

    TString weightString = "pow(mcPrimary_Energy_GeV,-";
    weightString += fabs(gamma);
    weightString += ")*mcOneWeight";
    TCut weight(weightString);

    assert( VerifyExpression(srcTree, weight.GetTitle()) );
    assert( VerifyExpression(srcTree, cut.GetTitle()) );

    srcTree->Draw(enString+">>"+hName,cut*weight, "goff");

    h1Energy.Scale(1./ h1Energy.GetSum() );

    // now move under/overflow content inside main region, if requested
    // (note that scale to 1/sum  correctly handled these bins too)
    if (optConstrainSignal_) {
      h1Energy.AddBinContent(1, h1Energy.GetBinContent(0));
      h1Energy.SetBinContent(0, 0);
      int lastBin = GetEnergyBins();
      h1Energy.AddBinContent(lastBin, h1Energy.GetBinContent(lastBin+1));
      h1Energy.SetBinContent(lastBin+1, 0);
    }
    
    for (int eBin=1; eBin<=GetEnergyBins(); eBin++) {
      hProbGamma_.SetBinContent(eBin, gBin, h1Energy.GetBinContent(eBin) );
    }
  }
 }

  // Find maximum for each energy
  for (int eBin=1; eBin <= GetEnergyBins(); ++eBin) {
    double max = 0.;
    for (int gBin=1; gBin<= GetGammaBins(); ++gBin) {
      double value = hProbGamma_.GetBinContent(eBin,gBin);
      if (value > max) { max = value; }
    }
    hProbGammaMax_.SetBinContent(eBin,max);
  }
}


void SimpleEnergyProb::SetTableBkg(const vector<I3Event>& eventVect) {
  if (!rangeIsSet_) {
    log_warn("Bkg histogram must be defined before can be reset.\n");
    return;
  }

  hProbBkg_.Reset(); // remove contents

  for (unsigned int i=0; i<eventVect.size(); ++i) {
    // Here is the connection to I3Event format
    double energy = eventVect[i].GetParams().energyValue;
    int eBin = GetEnergyBin(energy);  // this corrects for out-of-bound values
    hProbBkg_.AddBinContent(eBin);
  }
  hProbBkg_.SetEntries(eventVect.size());
  // ... so that when GetEntries() checks, the histogram will not look empty

  hProbBkg_.Scale( 1./hProbBkg_.GetSum() );

  HistFillIn(&hProbBkg_, energyNBinStartBackFill_);
} //
// 
//  TWO-VERSIONS RIGHT NOW
//
void SimpleEnergyProb::SetTableBkg(const EventPtrList& evList) {
  if (!rangeIsSet_) {
    log_warn("Bkg histogram must be defined before can be reset.\n");
    return;
  }

  hProbBkg_.Reset(); // remove contents

  for (int i=0; i<evList.GetSize(); ++i) {
    // Here is the connection to I3Event format
    const I3Event* ev = dynamic_cast<const I3Event*>(evList.GetEvent(i));
    double energy = ev->GetParams().energyValue;
    int eBin = GetEnergyBin(energy);  // this corrects for out-of-bound values
    hProbBkg_.AddBinContent(eBin);
  }
  hProbBkg_.SetEntries(evList.GetSize());
  // ... so that when GetEntries() checks, the histogram will not look empty

  hProbBkg_.Scale( 1./hProbBkg_.GetSum() );

  HistFillIn(&hProbBkg_, energyNBinStartBackFill_);
}


//
// GET FUNCTIONS()  CONST
//

int SimpleEnergyProb::GetEnergyBins() const {
  if (rangeIsSet_) { return hProbBkg_.GetNbinsX(); }
  log_warn("Range not set.\n");
  return 0;
}

double SimpleEnergyProb::GetEnergyMin() const {
  if (rangeIsSet_) { return hProbBkg_.GetXaxis()->GetXmin(); }
  log_warn("Range not set.\n");
  return 0;
}

double SimpleEnergyProb::GetEnergyMax() const {
  if (rangeIsSet_) { return hProbBkg_.GetXaxis()->GetXmax(); }
  log_warn("Range not set.\n");
  return 0;
}

int SimpleEnergyProb::GetGammaBins() const {
  if (rangeIsSet_) { return hProbGamma_.GetNbinsY(); }
  log_warn("Range not set.\n");
  return 0;
}

double SimpleEnergyProb::GetGammaMin() const {
  if (rangeIsSet_) { return hProbGamma_.GetYaxis()->GetXmin(); }
  log_warn("Range not set.\n");
  return 0;
}

double SimpleEnergyProb::GetGammaMax() const {
  if (rangeIsSet_) { return hProbGamma_.GetYaxis()->GetXmax(); }
  log_warn("Range not set.\n");
  return 0;
}


double SimpleEnergyProb::GetEnergyProbGamma(const Event& event, double gamma) const 
{
  if ( !rangeIsSet_ || !hProbGamma_.GetEntries() ) {
    log_error("Error: ProbGamma table not defined or filled.\n");
    return 0.;
  }

  int gBin = hProbGamma_.GetYaxis()->FindBin(gamma);
  if (gBin==0 || gBin==hProbGamma_.GetNbinsY()+1) {
    printf("gamma = %0.2f\n",gamma);
    log_fatal("Error: gamma out of range\n");
    return 0.;
  }

  double energy = (dynamic_cast<const I3Event&>(event)).GetParams().energyValue;
  int eBin = GetEnergyBin(energy);

  // nominal value of prob, before we interpolate
  double valueAtThisBin = hProbGamma_.GetBinContent(eBin, gBin);

  double gCenterThisBin = hProbGamma_.GetYaxis()->GetBinCenter(gBin);

  // if we are on the edge of the histogram, interpolate to zero to kill prob.
  if ( gBin==1 && gamma<gCenterThisBin ) {
    double lowEdge = hProbGamma_.GetYaxis()->GetBinLowEdge(gBin);
    return valueAtThisBin * (gamma-lowEdge)/(gCenterThisBin-lowEdge);
  }
  if ( gBin == GetGammaBins() && gamma>gCenterThisBin ) {
    double hiEdge = hProbGamma_.GetYaxis()->GetBinLowEdge(gBin+1);
    return valueAtThisBin * (hiEdge-gamma)/(hiEdge-gCenterThisBin);
  }

  int dir = 1;
  if (gamma<gCenterThisBin) { dir = -1; }
  
  double valueAtNextBin = hProbGamma_.GetBinContent(eBin, gBin+dir);
  double gCenterNextBin = hProbGamma_.GetYaxis()->GetBinCenter(gBin+dir);

  double fractionOffset = 
    fabs((gamma-gCenterThisBin)/(gCenterNextBin-gCenterThisBin));

  return valueAtThisBin + (valueAtNextBin-valueAtThisBin)*fractionOffset;
}


double SimpleEnergyProb::GetEnergyProbBkg(const Event& event) const {
  double energy = (dynamic_cast<const I3Event&>(event)).GetParams().energyValue;
  if ( hProbBkg_.GetEntries() > 0 ) {
    return hProbBkg_.GetBinContent( GetEnergyBin(energy) );
  }
  log_fatal("Error: no return value, ProbBkg table not defined.\n");
  return 0.;
}


double SimpleEnergyProb::GetEnergyMaxRatio(const Event& event) const {
  double bkgValue = GetEnergyProbBkg(event);
  if (bkgValue <= 0) {
    if ( hProbBkg_.GetEntries() == 0) {
      log_fatal("Ratio undefined because ProbBkg empty.\n");
    } else {
      log_fatal("Ratio undefined, prob bkg value is zero for this energy.\n");
    }
  }

  double energy=(dynamic_cast<const I3Event&>(event)).GetParams().energyValue;
  double gammaValue = hProbGammaMax_.GetBinContent( GetEnergyBin(energy) );
  if (gammaValue <= 0) {
    if ( hProbGammaMax_.GetEntries() <= 0) {
      log_fatal("Ratio undefined because ProbGammaMax empty.\n");
    }
  // THIS WARNING IS NOW HANDLED BY LlhEnergy, which can turn it on or off 
  //else {
  //  log_error("Error: gamma table is zero for this energy.  Max Ratio=0.\n");
  //}
  }

  return gammaValue / bkgValue;
}


TH1D* SimpleEnergyProb::GetHistProbGamma(double gamma) const {
  TString px_name = hProbGamma_.GetName();
  px_name += "_";
  px_name += gamma;
  int bin = hProbGamma_.GetYaxis()->FindBin(gamma);
  return hProbGamma_.ProjectionX(px_name,bin,bin);
}



//  HELPER FUNCTION FOR BACK-FILLING HISTOGRAMS

// default for nBinStart is specified in header file
void HistFillIn(TH1* h, int nBinStart) {
  if (h->GetNbinsY()>1 || h->GetNbinsZ()>1) {
    cout << "Error. HistFillIn only works with 1-D histogram.\n";
    return;
  }
  int nBins = h->GetNbinsX();

  int n=nBinStart;
  while (n<=nBins) {
    if (h->GetBinContent(n)>0) {
      ++n;
      continue;  // this bin has content, go to next one
    }
    // found empty bin... see if any later bins have content
    int nsearch = n+1;
    while (nsearch <= nBins) {
      if (h->GetBinContent(nsearch)>0) {
	for (int i=n; i<=nsearch; ++i) {
	  h->SetBinContent(i, h->GetBinContent(nsearch)/(nsearch-n+1));
	}
	break;  // while loop
      }
      ++nsearch;
    }
    n = nsearch+1;
  }

}
