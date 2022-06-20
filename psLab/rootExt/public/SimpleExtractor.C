
// 2008/3/1 : - Change TNtuple to TNtupleD
//              i.e., store doubles instead of floats
//            - Added monitoring... set to zero to turn off


// SimpleExtractor takes an input tree or ntuple (usually a chain of files),
// and a list of expressions,
// then writes a single output file with just these expressions
// (also applying any cuts specified)

// Example Usage:

/*

  SimpleExtractor s;


  s.AddVar("Nchan");

  s.AddVar("mcZenDeg", "MaxPrimary.dir_.zenith_*TMath::RadToDeg()");

  // Less commonly used //
  s.AddVar("reco3Zen", "reco_3_zenith", -1.);


  s.SetCut("paraboloidSigmaDeg<2.5 && paraboloidNDir>=9");
 

  TTree *inputTree;  // point this to the tree to be extracted
  char *outputFileName = "extracted.root";
  char *outputTreeName = "ntuple";  // analogous to typical name "tree" 
  char *outputTreeTitle = "My Extracted tree";

  s.MakeTree(inputTree, outputFileName, outputTreeName, outputTreeTitle);

*/

// Explanation of AddVar usage:
//
// 1-Parameter Example: name
//   s.AddVar("Nchan");
// output variable name is same as input tree variable name (or alias name)
//
// 2-Parameter Example: output name, input expression
//   s.AddVar("mcZenDeg", "MaxPrimary.dir_.zenith_*TMath::RadToDeg()");
// input variable name can be an expression, function call, etc.
//
// 3-Parameter Example:  output name, input expression, substitution value
//   s.AddVar("reco3Zen", "reco_3_zenith", -1.);
// If the entry has an empty value for the expression (rare), 
// specify what value should be stored in the output.
//
// Note that by default (if not specified as third parameter),
// an empty value in a tree entry will be stored as NaN in the output tree







// <limits> needed for using NaN's in case of empty entries
#include <limits>

#include <iostream>

#include "TCut.h"
#include "TFile.h"
#include "TNtupleD.h"
#include "TTree.h"
#include "TTreeFormula.h"

#include "rootExt/public/CountMonitor.h"
#include "rootExt/public/FunctionsRoot.h"
#include "rootExt/public/MakeNewName.h"



class SimpleExtractor {

 private:

  vector<TString> varNameVect_;  // the variable names
  vector<TString> varExpVect_;   // the expressions evaluated on each entry

  // values to substitute if the tree entry is empty (rarely used)
  vector<Double_t> subValueVect_;

  TCut cut_;

  double monitorPercent_;

 public:

  SimpleExtractor() { 
    cut_ = "1"; 
    monitorPercent_ = 10.;
  }

  // If varName == varExp, just use one parameter function call:
  void AddVar(const char* name) {  AddVar(name,name);  }

  void AddVar(const char* varName, const char* varExp) {
    Double_t subValue = numeric_limits<double>::quiet_NaN();
    AddVar(varName, varExp, subValue);
  }

  void AddVar(const char* varName, const char* varExp, Double_t subValue) {
    varNameVect_.push_back(varName);
    varExpVect_.push_back(varExp);
    subValueVect_.push_back(subValue);
  }


  void SetCut(TCut cut) { cut_ = cut; }

  // This is just a courtesy for the user.
  // Specify frequency of monitoring report, in percentage, 
  // i.e. maximum is 100.  (not 1.)
  // if ==0, no monitoring status is given
  void SetMonitorPercent(double percent) { monitorPercent_ = percent; }

  void MakeTree(TTree *inputTree, 
		char *outputFileName,
		char *outputTreeName, 
		char *OutputTreeTitle);
};




void SimpleExtractor::MakeTree (TTree *inputTree, 
				char *outputFileName,
				char *outputTreeName, 
				char *outputTreeTitle)
{

  // ROOT experts recommend this before passing TChain to TTreeFormula
  // (not having this line was causing problems on macs)
  inputTree->LoadTree(0);


  TFile f(outputFileName,"RECREATE");

  const int nVars = varNameVect_.size();

  TString varlist;
  TTreeFormula *formulaArray[nVars];
  bool emptyWarning[nVars];  // keep track of empty occurences
  bool emptyCutWarning = false;

  // In this loop:
  // 1) Make "varlist" of variable names, to initialize output ntuple
  // 2) Set up corresponding formulas to evaluate on input tree
  // 3) initialize warnings to false (= not yet activated)
  for (int i=0; i<nVars; ++i) {

    // varlist is one long string e.g. "var1:var2:var3"
    varlist += varNameVect_[i];
    if (i+1<nVars) {
      varlist += ":";
    }

    formulaArray[i] = 
      new TTreeFormula(MakeNewName("formula"),varExpVect_[i], inputTree);

    emptyWarning[i] = false;
  } 

  TNtupleD *ntuple = new TNtupleD(outputTreeName,outputTreeTitle,varlist);


  TTreeFormula *formulaCut = 
    new TTreeFormula(MakeNewName("cut"),cut_.GetTitle(), inputTree);



  Double_t values[nVars];  // result of expressions stored here

  int totalEntries = inputTree->GetEntries();

  CountMonitor countMon(monitorPercent_, totalEntries);
  if (monitorPercent_) { cout << "Monitoring Progress...: " << flush; }


  // Main Loop:

  Int_t tnumber = -1;   // starting value for keeping count of chained files

  for (int entry=0; entry < totalEntries ; ++entry) {

    inputTree->LoadTree(entry);

    // if we chained files and have moved to the next one, we need to update:
    if (tnumber != inputTree->GetTreeNumber() ) // tree changed, update leaves 
    { 
      tnumber = inputTree->GetTreeNumber();
      for (int i=0; i<nVars; ++i) {
	formulaArray[i]->UpdateFormulaLeaves();
      }
      formulaCut->UpdateFormulaLeaves();
    }


    // first check that the cut even evaluates to a result
    if ( formulaCut->GetNdata() == 0 ) {
      if ( !emptyCutWarning ) {
	cout << "\nWARNING: cut result is Empty for at least one entry.";
	cout << " (equivalent to failing cut.)\n";
	emptyCutWarning = true;
      }
    }
    else {
      // Okay, cut *can* be evaluated, so now actually evaluate it

      if (formulaCut->EvalInstance() ) {

	// Passed cut, so calculate each variable for this tree entry
	for (int i=0; i<nVars; ++i) {

	  // check that this variable is not empty
	  if ( formulaArray[i]->GetNdata() > 0 ) {
	    values[i] = formulaArray[i]->EvalInstance();
	  }
	  else { 
	    // variable for this entry is empty!  Handle accordingly: 
	    values[i] = subValueVect_[i];

	    // Let user know that an empty expression was encountered
	    if (!emptyWarning[i]) {
	      cout << "WARNING: Substituting ";
	      cout << subValueVect_[i] << " for empty expression \"";
	      cout << formulaArray[i]->GetTitle() << "\"\n";
	      emptyWarning[i] = true;
	    }
	  }
	}

	ntuple->Fill(values);
      } // closes if-cut-passes block

    } // close if-cut-not-empty block

    // if monitorPercent_ == 0, this never prints any status report
    countMon.UpdateCount();

  } // closes loop over all entries


  // Clean up

  if (monitorPercent_) { cout << endl; }

  for (int i=0; i<nVars; ++i) {
    delete formulaArray[i];
  } 
  delete formulaCut;

  f.Write();
  f.Close();

}
