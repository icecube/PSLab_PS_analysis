
// This is useful for testing possible speed gains with safe mode on or off,
// and comparing with GetValueFromTree.
// 
// Of course, you must COMPILE THIS FUNCTION
// if you want to see any useful difference at all!


#include <iostream>

#include "TStopwatch.h"

#include "TreeReader.h"
#include "FunctionsRoot.h"



void TreeReaderTimeTest(TTree *tree, const char* varexp, 
			bool safeMode=true, bool getValueTest=false) {
  double sum = 0.;
  cout << "TreeReader\n";
  TStopwatch ts1;
  TreeReader tr(tree, varexp);
  tr.SetSafeMode(safeMode);
  for (int i=0; i<tree->GetEntries(); ++i) {
    sum += tr.GetValue(i);
  }
  ts1.Print();
  cout << "Mean: " << sum / tree->GetEntries() << "\n\n";

  if (getValueTest) {
    cout << "GetValueFromTree\n";
    TStopwatch ts2;
    sum = 0.;
    for (int i=0; i<tree->GetEntries(); ++i) {
      sum += GetValueFromTree(tree,varexp,i);
    }
    ts2.Print();
    cout << "Mean: " << sum / tree->GetEntries() << "\n\n";
  }
}
