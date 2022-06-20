// To Do:
// Template the Get Histograms functions
// Delete old histogram upon re-plotting


#include <iostream>

#include "TCut.h"
#include "TH1.h"
#include "TPad.h"
#include "TROOT.h"
#include "TString.h"
#include "TTree.h"

#include "rootExt/public/MakeNewName.h"

#include "rootDev/public/SetPad.h"


static const int NMAX = 100;

// recall: static global variables have scope only for the remainder 
// of the file they are in


static const TString replaceString = "$$";

// wherever this occurs in varexp, it will be replaced with a unique 
// name for the object being drawn.  the name will be accessible via
// GetName(int i);
//
// I believe "$$" would never occur in a normal varexp, but perhaps
// if there is a clash in the future, the replaceString could be changed

static const TString sumString = "$SUM$";

// Example:   myMH.SetVarexp(5,"$SUM$ 0,1,2");  
//  Entry 5 will be a sum of entries 0, 1, 2




class MH {

private:
  bool active_[NMAX];
  TTree *tree_[NMAX];
  TString varexp_[NMAX];
  TCut selection_[NMAX];
  TCut cut2_[NMAX];
  TString option_[NMAX];
  int color_[NMAX];
  int linewidth_[NMAX];
  int linestyle_[NMAX];
  TString label_[NMAX];

  TString nameOfDrawnObject_[NMAX];

  bool verbose_;
  bool optLogy_;
  bool optYmin_;
  bool optYmaxAuto_;
  bool optYmaxSpecific_;
  double ymin_;
  double ymax_;

  double norm_;
public:

  MH() {
    for (int i=0; i<NMAX; ++i) {
      active_[i] = false;
      tree_[i] = NULL;
      varexp_[i] = "";
      selection_[i] = "";
      cut2_[i] = "";
      option_[i] = "";
      color_[i] = i+1;
      linewidth_[i] = 2;
      linestyle_[i] = 1;
      label_[i] = "";
    }
    verbose_ = true;
    optLogy_ = false;
    optYmin_ = false;
    optYmaxAuto_ = true;
    optYmaxSpecific_ = false;
    ymin_ = 0.;
    ymax_ = 0.;  // if ymax<=ymin, pad is rescaled to contain all hists
    norm_ = 0.;
  }

  void SetVerbose(bool verbose) {verbose_ = verbose;}

  void SetLogy(bool optLogy) {
    optLogy_ = optLogy; 
  }
  void SetYMin(double ymin) {
    optYmin_ = true;
    ymin_ = ymin;
  }
  void SetYMaxAuto(bool optYmax) {
    optYmaxAuto_ = optYmax;
    optYmaxSpecific_ = false;
  }
  void SetYMax(double ymax) {
    optYmaxAuto_ = false;
    optYmaxSpecific_ = true;
    ymax_ = ymax;
  }
  void SetNorm(double norm=1) {
    norm_ = norm;
  }
  void SetNormOff() {
    norm_ = 0.;
  }


  bool Valid(int i) {
    if (i>=0 && i<NMAX) { return true; }
    printf("i=%d is out of range [0,%d]\n",i,NMAX);
    return false;
  }

  void Activate(int i) { if (Valid(i)) { active_[i] = true; } }
  void DeActivate(int i) { if (Valid(i)) { active_[i] = false; } }
  bool IsActive(int i) {
    if (Valid(i)) { return active_[i]; }
    return false;
  }


  void SetEntry(int i, 
		TTree *tree, 
		TString varexp, 
		TCut selection = "", 
		TCut cut2 = "",
		TString option = "") {
    if (Valid(i)) {
      tree_[i] = tree;
      varexp_[i] = varexp;
      selection_[i] = selection;
      cut2_[i] = cut2;
      option_[i] = option;
      // by default, setting an entry with this method will activate it
      Activate(i);
    }
  }

  // SelectActive("2,4,6"), // only 2,4, and 6 will be active

  int SelectActive(TString activeString) {
    for (int i=0; i<NMAX; ++i) {
      DeActivate(i);
    }
    TObjArray* activeArray = activeString.Tokenize(",");
    if (!activeArray) return 0;

    for (int i=0; i < activeArray->GetEntries(); ++i) {
      TString s = ((*activeArray)[i])->GetName();
      this->Activate(s.Atoi());
    }
    return activeArray->GetEntries();
  }


  void SetTree(int i, TTree *tree) { if (Valid(i)) tree_[i] = tree; }
  void SetTree(TTree *tree) { SetTree(0, tree); }

  void SetVarexp(int i, TString varexp) {if (Valid(i)) varexp_[i] = varexp; }
  void SetVarexp(TString varexp) { SetVarexp(0, varexp); }

  void SetSelection(int i, TCut selection) 
  { if (Valid(i)) selection_[i] = selection; }
  void SetSelection(TCut selection) { SetSelection(0, selection); }

  void SetCut2(int i, TCut cut2) {if (Valid(i)) cut2_[i] = cut2; }
  void SetCut2(TCut cut2) { SetCut2(0, cut2); }

  void SetOption(int i, TString option) { if (Valid(i)) option_[i] = option; }
  void SetOption(TString option) { SetOption(0, option); }


  void SetColor(int i, int color) { if (Valid(i)) color_[i] = color; }

  void SetLineWidth(int i, int width) { if (Valid(i)) linewidth_[i] = width; }

  void SetLineStyle(int i, int style) { if (Valid(i)) linestyle_[i] = style; }


  void Draw(TString title);
  void Draw() { Draw(""); }
  // recall that for title you can do: "Main Title;X-axis label;Y-axis label"


  TString GetName(int i) {
    if (Valid(i)) { return nameOfDrawnObject_[i]; }
    return "";
  }

  TH1F* GetTH1F(int i) {
    if (Valid(i)) {
      return dynamic_cast<TH1F*>(gROOT->FindObject(GetName(i)));
    }
    return NULL;
  }



};


void MH::Draw(TString title) {

  // Clear names first, so even if Draw fails, there are no accidental
  // references to objects drawn earlier
  for (int i=0; i<NMAX; ++i) {
    nameOfDrawnObject_[i] = "";
  }

  bool firstAlreadyDrawn = false;

  
  for (int i=0; i<NMAX; ++i) {

    if (active_[i]) {

      // Start with common settings

      TTree *thisTree = tree_[0];
      TString thisVarexp = varexp_[0];
      TCut thisSelection = selection_[0];
      TCut thisCut2 = cut2_[0];
      TString thisOption = option_[0];


      // Change to other settings, if any are specified

      if (i>0) {
	if (tree_[i] != NULL) thisTree = tree_[i];
	if ( strcmp(varexp_[i],"") ) thisVarexp = varexp_[i];
	if ( strcmp(selection_[i].GetTitle(),"") ) thisSelection = selection_[i];
	if ( strcmp(cut2_[i].GetTitle(),"") ) thisCut2 = cut2_[i];
	if ( strcmp(option_[i],"") ) thisOption = option_[i];
      }


      // replace with unique name for object, if requested

      if (thisVarexp.Contains(replaceString)) {
	nameOfDrawnObject_[i] = MakeNewName("MH");
	thisVarexp.ReplaceAll(replaceString,nameOfDrawnObject_[i]);
      }



      if (thisVarexp.Contains(sumString)) {
	TH1 *hSum = NULL;
	thisVarexp.ReplaceAll(sumString,"");
	const TObjArray* objList = thisVarexp.Tokenize(",");
	if (!objList) {
	  printf("ERROR: $SUM$ was requested but not understood.\n");
	  return;
	}
	else {
	  for (int j = 0; j < objList->GetEntries(); ++j) {
	    TString s = ((*objList)[j])->GetName();
	    int hEntry = s.Atoi(); // s is a string such as "1" or "4"
	    if (hEntry>=i) {
	      printf("ERROR: $SUM$ called for histogram entry %d ",hEntry);
	      printf("not yet plotted.\n");
	      return;
	    }
	    TString hName = GetName(hEntry);
	    TH1 *hOrig = dynamic_cast<TH1*>(gROOT->FindObject(hName));
	    if (!hOrig) {
	      printf("ERROR: $SUM$ called for a histogram entry %d ",hEntry);
	      printf("not found.\n");
	      return;
	    }
	    if (j==0) {
	      nameOfDrawnObject_[i] = MakeNewName("MH");
	      hSum = dynamic_cast<TH1*>(hOrig->Clone(nameOfDrawnObject_[i]));
	      hSum->SetLineColor(color_[i]);
	      hSum->SetLineWidth(linewidth_[i]);
	      hSum->SetLineStyle(linestyle_[i]);
	    } else {
	      hSum->Add(hOrig);
	    }
	  }
	}
	hSum->Draw(thisOption+"same");
	gPad->Update();
	continue;
      }



      if (firstAlreadyDrawn) {
	thisOption = thisOption+"same";
      } else {
	firstAlreadyDrawn = true;
      }


      thisTree->SetLineColor(color_[i]);
      thisTree->SetMarkerColor(color_[i]); // in case Graph is drawn
      thisTree->SetLineWidth(linewidth_[i]);
      thisTree->SetLineStyle(linestyle_[i]);

      TCut finalSelection = thisSelection * thisCut2;

      if (verbose_) {
	cout << i << ": ";
	cout << "\"" << thisVarexp << "\", ";
	cout << "\"" << finalSelection << "\", ";
	cout << "\"" << thisOption << "\"\n";
      }

      thisTree->Draw(thisVarexp, finalSelection, thisOption);
      if (norm_ != 0) {
	GetTH1F(i)->Scale( norm_ / GetTH1F(i)->GetSum() );
      }
      gPad->Modified();
    }
  }

  // make nicer

  if (firstAlreadyDrawn) {   // just make sure something is there
    if (title != "") {
      GetFirstHistInPad()->SetTitle(title);
    }
    if (optLogy_) { 
      gPad->SetLogy();
    }
    if (optYmin_) {
      SetPadYMin(ymin_); 
    }
    if (optYmaxAuto_) {
	SetPadYMax();
    }
    if (optYmaxSpecific_) {
      SetPadYMax(ymax_);
    }
    gPad->Modified();
  }

}

