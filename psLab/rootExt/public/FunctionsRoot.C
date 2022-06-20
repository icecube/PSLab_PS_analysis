#include "rootExt/public/FunctionsRoot.h"

#include <iostream>


bool VerifyExpression(TTree *tree, const char* expression) {
  TTreeFormula formula("","1", tree);  // initialize with expression="1", (OK)
  if (formula.Compile(expression)) {   // now test expresssion
    return false;
  }
  return true;
}


// Return the value of a particular variable (or expression) evaluated
// on the first entry in a TTree.  For example:
//
// "radius_cm = GetValueFromTree(tree, "InjectionRadius*100.");

double GetValueFromTree(TTree *tree, const char *varexp) {
  return GetValueFromTree(tree, varexp, 0);
}

// Or, evaluate on the i-th entry of the TTree.

double GetValueFromTree(TTree *tree, const char *varexp, int entry) {
  TTreeFormula formula(MakeNewName("tempFormula"), varexp, tree);
  tree->LoadTree(entry);
  formula.UpdateFormulaLeaves();
  // check if this entry is empty
  if (formula.GetNdata() == 0) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return formula.EvalInstance();
}

// specify substituteValue to return if tree entry is empty (rare)

double GetValueFromTree(TTree *tree, const char *varexp, int entry, 
			double substituteValue) {
  TTreeFormula formula(MakeNewName("tempFormula"), varexp, tree);
  tree->LoadTree(entry);
  formula.UpdateFormulaLeaves();
  // check if this entry is empty
  if (formula.GetNdata() == 0) {
    return substituteValue;
  }
  return formula.EvalInstance();
}





class HistIntegral {
  TH1* hIntegral;
  bool optNormalize;
  bool optContained;

public:

  HistIntegral(bool normalize = true, bool contained = false) {
    hIntegral = NULL;
    optNormalize = normalize;
    optContained = contained;
  }

  ~HistIntegral() {
    if (hIntegral) { delete hIntegral;}
  }

  void SetNormalize(bool normalize) { optNormalize = normalize; }

  void SetContained(bool contained) { optContained = contained; }

  TH1* Get() {
    if (hIntegral) { return hIntegral; }
    return NULL;
  }


  TH1* Ascending(const TH1* hMain) 
  {
    if (hIntegral) { delete hIntegral; }

    if (!hMain) {
      cout << "Error: input hist. for HistIntegral::Ascending  not found.\n";
      return NULL;
    }


    // Start by giving same binning of hMain to hIntegral

    hIntegral = dynamic_cast<TH1*> 
      ( hMain->Clone(MakeNewName(hMain->GetName())+"Integral") );
    hIntegral->Reset();
    // Don't inherit rebinning, because we want to be able to fill
    // underflow and overflow bins correctly, without histogram expanding
    hIntegral->SetCanExtend(TH1::kAllAxes);

    Int_t nBins = hMain->GetNbinsX();

    Int_t firstBin, lastBin;
    if (optContained) {
      firstBin=1;
      lastBin=nBins;
    } else {
      firstBin=0;        // the underflow bin
      lastBin=nBins+1;   // the overflow bin
    }

    Double_t sum=0.;
    for (Int_t i=firstBin; i<=lastBin; i++) {
      sum += hMain->GetBinContent(i);
      hIntegral->SetBinContent(i,sum);
    }
    // OLD:
    //    h_new->SetEntries(sum); // somewhat meaningless ....

    if (optNormalize) {
      hIntegral->Scale(1./sum);
    }

    if (optContained && hMain->GetBinContent(0) != 0.)
      cout << "Warning: there was underflow content, which not summed.\n";
    if (optContained && hMain->GetBinContent(nBins+1) != 0.)
      cout << "Warning:  there was overflow content, which not summed.\n";

    return hIntegral;
  }

  /*
  TH1* Descending() {
    if (hIntegral) { delete hIntegral; }
    // start by making hIntegral a bin-for-bin copy of hMain
    hIntegral = dynamic_cast<TH1*> 
      ( hMain->Clone(MakeNewName(hMain->GetName())+"Integral") );

    DescendingCumulate(hMain, hIntegral);
    return hIntegral;
  }

  TH1* DescendingNorm() {
    Descending();
    hIntegral->Scale(1./hIntegral->GetEntries());
    return hIntegral;
  }
  */

};



// Older version



TH1* AscendingCumulate(TH1 *h, const char* options) // default options=""
{
  int NBinMin = 0; // underflow bin
  int NBinMax = h->GetNbinsX()+1; // overflow bin
  TH1 *hnew = dynamic_cast<TH1*> (h->Clone());
  // Don't inherit rebinning, because we want to be able to fill
  // underflow and overflow bins correctly, without histogram expanding
  hnew->SetCanExtend(TH1::kAllAxes);

  double sum=0.;
  for (Int_t i=NBinMin; i<=NBinMax; i++) {
    sum += h->GetBinContent(i);
    hnew->SetBinContent(i,sum);
  }
  hnew->SetEntries(sum); // in case anyone asks...
  if (strstr(options,"Scale")) {
    hnew->Scale(1./sum);
  }
  return hnew;
}


TH1* DescendingCumulate(TH1 *h, const char* options) // default options=""
{
  int NBinMin = 0; // underflow bin
  int NBinMax = h->GetNbinsX()+1; // overflow bin
  TH1 *hnew = dynamic_cast<TH1*> (h->Clone());
  // Don't inherit rebinning, because we want to be able to fill
  // underflow and overflow bins correctly, without histogram expanding
  hnew->SetCanExtend(TH1::kAllAxes);

  double sum=0.;
  for (Int_t i=NBinMax; i>=NBinMin; i--) {
    sum += h->GetBinContent(i);
    hnew->SetBinContent(i,sum);
  }
  hnew->SetEntries(sum); // in case anyone asks...
  if (strstr(options,"Scale")) {
    hnew->Scale(1./sum);
  }
  return hnew;
}



// Oldest Version



// Make Cumulative Histograms

void AscendingCumulate(TH1 *h_old, TH1 *h_new)
{
  Double_t sum=0.;
  Int_t nBins = h_old->GetNbinsX();

  for (Int_t i=1; i<= nBins; i++) {
    sum += h_old->GetBinContent(i);
    h_new->SetBinContent(i,sum);
  }
  h_new->SetEntries(sum); // somewhat meaningless ....

  if (h_old->GetBinContent(0) != 0. || h_old->GetBinContent(nBins+1) != 0.) {
    printf(
     "Warning: cumulating histogram with non-zero underflow/overflow bin\n");
  }
}


void DescendingCumulate(TH1 *h_old, TH1 *h_new)
{
  Double_t sum=0.;
  Int_t nBins = h_old->GetNbinsX();

  for (Int_t i=nBins ; i>= 1; i--) {
    sum += h_old->GetBinContent(i);
    h_new->SetBinContent(i,sum);
  }
  h_new->SetEntries(sum); // somewhat meaningless ....

  if (h_old->GetBinContent(0) != 0. || h_old->GetBinContent(nBins+1) != 0.) {
    printf(
     "Warning: cumulating histogram with non-zero underflow/overflow bin\n");
  }
}



// TO DO:  Actually, these problems do not always mean you can't calculate
// the median... routine could be made smarter to figure out when....

Double_t histMedian(TH1 *h)
{
  // check for potential problems
  if (h->GetBinContent(0) != 0.) {
    cout << "ERROR: median not valid, due to contents in underflow bin\n";
    return 0.;
  }
  if (h->GetBinContent(h->GetNbinsX()+1) != 0.) {
    cout << "ERROR: median not valid, due to contents in overflow bin\n";
    return 0.;
  }

  Double_t half =0.5;
  Double_t median_value;
  h->GetQuantiles(1,&median_value,&half);
  return median_value;
}
