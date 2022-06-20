#include "rootExt/public/TreeReader.h"

#include <limits>

#include "rootExt/public/MakeNewName.h"



TreeReader::TreeReader(TTree *tree) : 
  tree_(tree),
  treeNumber_(-1),
  formula_(NULL),
  safeMode_(true)
{ }

TreeReader::TreeReader(TTree *tree, const char *varexp) : 
  tree_(tree),
  treeNumber_(-1),
  formula_(NULL),
  safeMode_(true)
{
  SetVarexp(varexp);
}


void TreeReader::SetVarexp(const char *varexp) 
{
  if (formula_) { delete formula_; }
  formula_ = new TTreeFormula(MakeNewName("TreeReaderFormula"),varexp, tree_);
  treeNumber_ = -1;
}


double TreeReader::GetValue(int entry)
{
  tree_->LoadTree(entry);

  if ( safeMode_ ||
       tree_->GetTreeNumber() != treeNumber_ ) {
       // If we chained files and have moved to a different file
       // since previous call, then we need to update formula leaves:
    treeNumber_ = tree_->GetTreeNumber();
    formula_->UpdateFormulaLeaves();
  }

  // check if this entry is empty
  if (formula_->GetNdata() == 0) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  return formula_->EvalInstance();
}
