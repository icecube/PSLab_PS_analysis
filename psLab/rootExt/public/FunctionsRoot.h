#ifndef ROOTEXT_FUNCTIONSROOT_H_
#define ROOTEXT_FUNCTIONSROOT_H_


// <limits> is needed for NaN functionality in GetValueFromTree
#include <limits>

#include "TH1.h"
#include "TString.h"
#include "TTree.h"
#include "TTreeFormula.h"

#include "rootExt/public/MakeNewName.h"


// return false if expression cannot be evaluated correctly in this tree
bool VerifyExpression(TTree *tree, const char* expression);

double GetValueFromTree(TTree *tree, const char *varexp);
double GetValueFromTree(TTree *tree, const char *varexp, int entry);
double GetValueFromTree(TTree *tree, const char *varexp, int entry,
                        double substituteValue);


TH1* AscendingCumulate(TH1 *h, const char* options="");
TH1* DescendingCumulate(TH1 *h, const char* options="");

void AscendingCumulate(TH1 *h_old, TH1 *h_new);
void DescendingCumulate(TH1 *h_old, TH1 *h_new);

Double_t histMedian(TH1 *h);

#endif // ROOTEXT_FUNCTIONSROOT_H_

