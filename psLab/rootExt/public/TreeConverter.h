#ifndef ROOTEXT_TREECONVERTER_H_
#define ROOTEXT_TREECONVERTER_H_

#include <vector>
#include <map>

#include "TString.h"
#include "TTree.h"
#include "TTreeFormula.h"

#include "rootExt/public/MakeNewName.h"


class TreeConverter {
  int nCuts_;
  vector<TString> cutStringVect_;
  vector<TTreeFormula*> cutFormulaVect_;

  int nVars_;
  vector<TString> varNameVect_;
  vector<TString> varExpVect_;
  map<TString,TTreeFormula*> varFormulaMap_;
  map<TString,double> varResultMap_;
  map<TString,double>* resultPtr_;

  TTree *tree_;
  int tnumber_;
  int nTotalEntries_;
  int entry_;

  bool monitor_;

 public:

  TreeConverter() {
    nCuts_ = 0;
    nVars_ = 0;
    resultPtr_ = NULL;
    tree_ = NULL;
    nTotalEntries_ = 0;
    tnumber_ = -1;
    monitor_ = true;
  }

  ~TreeConverter() {
    for (unsigned int i=0; i<cutFormulaVect_.size(); ++i) {
      delete cutFormulaVect_[i];
    }
    map<TString,TTreeFormula*>::iterator it;
    for (it = varFormulaMap_.begin(); it != varFormulaMap_.end(); it++) {
      delete (it->second);
    }
  }

  void ResetEntry() {entry_ = -1;}

  void AddCut(const char* cut) { 
    cutStringVect_.push_back(cut); 
  }

  void AddVar(const char* name, const char* exp) { 
    varNameVect_.push_back(name);
    varExpVect_.push_back(exp);
  }

  void SetMonitor(bool monitor) {
    monitor_ = monitor;
  }

  void SetTree(TTree *tree);

  bool ReadIthEntry(int entry);

  bool GetNextEntry();

  const map<TString,double>* GetResultMap() {return resultPtr_;}
  
};

#endif // ROOTEXT_TREECONVERTER_H_

