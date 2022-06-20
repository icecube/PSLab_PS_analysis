#ifndef ROOTEXT_FUNCTIONSROOT_H_

//
// Retrieve information from TTree entries directly, 
// using same syntax as formulae in TTree->Draw() commands.
//
//
// EXAMPLE USAGE:    (assume you have a TTree *myTree  pointer already)
//
//  double sum = 0;
//  int nEvents = myTree->GetEntries();
//
//  for (int i=0; i<nEvents; ++i) {
//    zenError = GetValueFromTree(myTree, "recoZenith - MCPrimary_Zenith", i);
//    sum += zenError;
//  }
//  avg_zen_error = sum / nEvents;
//


// <limits> is needed for NaN functionality in GetValueFromTree
#include <limits>

#include "TTree.h"
#include "TTreeFormula.h"

double GetValueFromTree(TTree *tree, const char *varexp, int entry) {
  TTreeFormula formula("GetValueFromTree_tempFormula", varexp, tree);
  tree->LoadTree(entry);
  formula.UpdateFormulaLeaves();
  // check if this entry is empty
  if (formula.GetNdata() == 0) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return formula.EvalInstance();
}


#endif
