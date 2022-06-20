#include <vector>
#include <iostream>

#include "TCut.h"
#include "TEventList.h"
#include "TTree.h"
#include "TString.h"

#include "rootExt/public/MakeNewName.h"


class EventController {

public:
  vector<TTree*> treeVect;
  vector<TEventList*> evListVect;


  EventController() { }

  void ClearLists() {
    for (unsigned int i=0; i<treeVect.size(); ++i) {
      if (evListVect[i]) {
	delete evListVect[i];
	evListVect[i] = NULL;
      }
      treeVect[i]->SetEventList(NULL);
    }
  }

  ~EventController() {
    ClearLists();
  }


  void AddTree(TTree *tree) {
    treeVect.push_back(tree);
    evListVect.push_back(NULL);
  }

  
  void SetLists(TCut cut, bool verbose=true) {
    ClearLists();
    for (unsigned int i=0; i<treeVect.size(); ++i) {
      TString listName = MakeNewName("EventControlList");
      TEventList *listPtr = new TEventList(listName);
      treeVect[i]->Draw(">> "+listName, cut);
      evListVect[i] = listPtr;
      treeVect[i]->SetEventList(evListVect[i]);
      if (verbose) {
	cout << "Entries selected for tree " << i;
	cout << "(" << treeVect[i]->GetTitle() << ") = ";
	cout << evListVect[i]->GetN() << endl;
      }
    }
  }

};


