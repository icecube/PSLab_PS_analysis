#include "llhSN/public/SimpleProb.h"

#include "TAxis.h"
#include "TCut.h"
#include "TEventList.h"
#include "TH1.h" // for HistFillIn function

#include "rootExt/public/FunctionsRoot.h" // for VerifyExpresion()
#include "rootExt/public/log_report.h"  // for log_warns, etc.
#include "rootExt/public/TreeReader.h"
#include "rootExt/public/MakeNewName.h"



//   PRIVATE FUNCTIONS   //


int SimpleProb::GetBin(double energy) const {
  if (rangeIsSet_) {
    int eBin = hProbBkg_.GetXaxis()->FindBin(energy);
    // if energy is in underflow or overflow region, move it back w/in range
    if (eBin==0) { eBin++; }
    if (eBin==GetBins()+1) { eBin--; }
    return eBin;
  }
  log_error("Error: no return value, energy range is not defined.\n");
  return -1;
}


//   PUBLIC FUNCTIONS   //


//
// SET FUNCTIONS
//

void SimpleProb::SetRangeAndBackFill(int nBins, double sMin, double sMax,
	 int nBinStartBackFill) {
  // Reset empties any previous contents.  GetSum() now equals zero.
  hProbBkg_.Reset();
  hProbBkg_.SetBins(nBins, sMin, sMax);

  nBinStartBackFill_ = nBinStartBackFill;
  rangeIsSet_ = true;
}


void SimpleProb::SetTableBkg(const vector<SNEvent> eventVect) {

  hProbBkg_.Reset(); // remove contents

  for (unsigned int i=0; i<eventVect.size(); ++i) {
    // Here is the connection to I3Event format
    double s = eventVect[i].GetStrength();
    int eBin = GetBin(s);  // this corrects for out-of-bound values
    hProbBkg_.AddBinContent(eBin);
  }
  hProbBkg_.SetEntries(eventVect.size());
  // ... so that when GetEntries() checks, the histogram will not look empty

  hProbBkg_.Scale( 1./hProbBkg_.GetSum() );

  HistFillIn(&hProbBkg_, nBinStartBackFill_);
}

//
// 
//  TWO-VERSIONS RIGHT NOW
//
void SimpleProb::SetTableBkg(const EventPtrList& evList) {
  if (!rangeIsSet_) {
    log_warn("Bkg histogram must be defined before can be reset.\n");
    return;
  }

  hProbBkg_.Reset(); // remove contents

  for (int i=0; i<evList.GetSize(); ++i) {
    // Here is the connection to I3Event format
    const SNEvent* ev = dynamic_cast<const SNEvent*>(evList.GetEvent(i));
    double s = ev->GetStrength();
    int eBin = GetBin(s);  // this corrects for out-of-bound values
    hProbBkg_.AddBinContent(eBin);
  }
  hProbBkg_.SetEntries(evList.GetSize());
  // ... so that when GetEntries() checks, the histogram will not look empty

  hProbBkg_.Scale( 1./hProbBkg_.GetSum() );

  HistFillIn(&hProbBkg_, nBinStartBackFill_);
}


//
// GET FUNCTIONS()  CONST
//

int SimpleProb::GetBins() const {
  if (rangeIsSet_) { return hProbBkg_.GetNbinsX(); }
  log_warn("Range not set.\n");
  return 0;
}

double SimpleProb::GetMin() const {
  if (rangeIsSet_) { return hProbBkg_.GetXaxis()->GetXmin(); }
  log_warn("Range not set.\n");
  return 0;
}

double SimpleProb::GetMax() const {
  if (rangeIsSet_) { return hProbBkg_.GetXaxis()->GetXmax(); }
  log_warn("Range not set.\n");
  return 0;
}


double SimpleProb::GetProbBkg(const SNEvent * event) const {
  //double strength = (dynamic_cast<const SNEvent>(event)).GetStrength();
  double strength = event->GetStrength();
  if ( hProbBkg_.GetEntries() > 0 ) {
    return hProbBkg_.GetBinContent( GetBin(strength) );
  }
  log_fatal("Error: no return value, ProbBkg table not defined.\n");
  return 0.;
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
