#ifndef ROOTEXT_MAKENEWNAME_H_
#define ROOTEXT_MAKENEWNAME_H_

#include "TDirectory.h"
#include "TString.h"

// Take input name, check if it is already in use by another object or not.
// - If not, return the name.
// - If it is, make a new name by incrementing the trailing number until
// a clash is not found.
//
// Note that TStrings automatically cast to and from char* and string,
// so an expression like:
//
//   TH1D(MakeNewName("myHist"),"This is my histogram",10,0,1);
//
// should work smoothly where char* is expected.


TString MakeNewName(TString name) {
  int n = 0;
  TString newName = name;
  while ( gDirectory->FindObject(newName) ) {
    ++n;
    newName = name + "_";
    newName += n;
  }
  return newName;
}


#endif // ROOTEXT_MAKENEWNAME_H_
