#ifndef ROOTEXT_TREEREADER_H_
#define ROOTEXT_TREEREADER_H_

#include "TTree.h"
#include "TTreeFormula.h"

// TreeReader offers improved speed over the GetValueFromTree function,
// because the formula is created only once.  Example Usage:
// 
// TreeReader tr(myTree, "log10(energy)");
// for (int i=0; i<myTree->GetEntries(); ++i) {
//   sum = tr.GetValue(i);
// }
// meanLogEnergy = sum / myTree->GetEntries(); 
//
//
// It can be made even faster by turning SafeMode off, read on:
//
// SafeMode is true by default.
// The Formula leaves are updated every call of GetValue.
//
// If SafeMode is turned off, GetValue still checks to see if the treeNumber
// (i.e. chained file) has changed since the last call, and updates the leaves
// if so.
// The danger is if the treeNumber is the same as during the last call, but
// meanwhile another piece of code has changed the chained file.
// Example:
//
// TreeReader tr1(myChainPtr, "energy");
// TreeReader tr2(myChainPtr, "energy");
// tr1.SetSafeMode(false);
// tr2.SetSafeMode(false);
// tr1.GetValue(0);
// tr2.GetValue(0);
// tr1.GetValue(9999);  // different file in the chain!
// tr2.GetValue(0);     // tr2 doesn't "see" that the chain has been
//                      // pointed somewhere else in the meantime, 
//                      // because this entry is the same as last time.
//                      // So the formula leaves are not updated.
//
// The behavior in this case ranges from retrieving the wrong value, to 
// an execution error.
//
// Turning SafeMode off can speed up GetValue by an order of magnitude,
// but make sure you know what you are doing.  It should be okay in 
// situations where the chain is accessed *only* by TreeReaders, acting
// *always* in concert on the same entries, e.g.:
//
// for (int i=0; i<10000; ++i) {
//   double x = tr1.GetValue(i);
//   double y = tr2.GetValue(i);
//   ...
// }
//


class TreeReader {
 public: 
  TreeReader(TTree *tree);
  TreeReader(TTree *tree, const char *varexp);

  virtual ~TreeReader() { if (formula_) { delete formula_; } }

  void SetSafeMode(bool mode) { safeMode_ = mode; }
  bool GetSafeMode() const { return safeMode_; }

  void SetVarexp(const char *varexp);

  double GetValue(int entry);

 private:
  TTree *tree_;
  int treeNumber_;
  TTreeFormula *formula_;
  bool safeMode_;
};

#endif // ROOTEXT_TREEREADER_H_
