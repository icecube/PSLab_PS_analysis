#include "rootExt/public/TreeConverter.h"

#include <iostream>


void TreeConverter::SetTree(TTree *tree) {

  if (tree_) {
    cout << "Error. Tree has already been set in TreeConverter.\n";
    return;
  }

  tree_ = tree;

  nCuts_ = cutStringVect_.size();

  cutFormulaVect_.clear();
  for (int i=0; i<nCuts_; ++i) {
    TString name = TString("cut_formula_");
    name += i;
    cutFormulaVect_.push_back(
      new TTreeFormula(MakeNewName(name), cutStringVect_[i], tree_) );
  }

  nVars_ = varNameVect_.size();

  varFormulaMap_.clear();
  varResultMap_.clear();
  for (int i=0; i<nVars_; ++i) {
    // need MakeNewName
    TString name = varNameVect_[i];
    varFormulaMap_[name] = 
      new TTreeFormula(MakeNewName(name), varExpVect_[i], tree_);
    varResultMap_[name] = 0.;
  } 


  nTotalEntries_ = tree_->GetEntries();
  if (monitor_) {
    cout << "Selecting events from " << nTotalEntries_ << " entries\n";
  }

  tnumber_ = -1;   // starting value for keeping count of chained files
  entry_ = -1;
}



bool TreeConverter::ReadIthEntry(int ie) {

  resultPtr_ = NULL;
  // this will only get set if and when an event passes all cuts

  map<TString,TTreeFormula*>::iterator it; 

  tree_->LoadTree(ie);

  // UPDATE FORMULA LEAVES
  //
  // If we chained files and have moved from one file to the next,
  // then we need to update:
  if (tnumber_ != tree_->GetTreeNumber() ) // tree changed, update leaves 
  { 
    tnumber_ = tree_->GetTreeNumber();
    
    for (int i=0; i<nCuts_; ++i) {
      cutFormulaVect_[i]->UpdateFormulaLeaves();
    }

    for (it = varFormulaMap_.begin(); it != varFormulaMap_.end(); it++) {
      (it->second)->UpdateFormulaLeaves();
    }
  }

    
  // EVALUATE CUTS, AND STORE RESULTS IF PASS 

  bool pass=true;
  int i=0;
  while (pass && i<nCuts_) {
    pass = cutFormulaVect_[i]->EvalInstance();
    i++;
  }

  if (pass) {
    for (it = varFormulaMap_.begin(); it != varFormulaMap_.end(); it++) {
      TString varname = it->first;
      varResultMap_[varname] = (it->second)->EvalInstance();
    }

    resultPtr_ = &varResultMap_;
  }

  return pass;
}



bool TreeConverter::GetNextEntry() {

  resultPtr_ = NULL;
  // this will only get set when an event passes all cuts

  bool pass = false;

  while (!pass && entry_ < (nTotalEntries_-1) ) {
    ++entry_;

    // Give progress report...
    if (monitor_ && entry_%(nTotalEntries_/10) == 0) {
      cout << entry_/(nTotalEntries_/100) << "% " << flush;
    }

    pass = ReadIthEntry(entry_);
    // resultPtr_ gets set if entry passed cuts
  }

  return pass;

}
