#ifndef ROOTDEV_SETPAD_H_
#define ROOTDEV_SETPAD_H_


#include <iostream>

#include "TH1.h"
#include "TList.h"
#include "TObject.h"
#include "TPad.h"


TH1* GetFirstHistInPad() {
  TList *list = gPad->GetListOfPrimitives();
  TIter iter(list);

  TObject *obj;
  TH1 *h;

  while ( (obj = iter.Next()) ) {
    h = dynamic_cast<TH1*> (obj);
    if (h) {
      return h;
    }
  }

  cout << "ERROR: No histogram found in pad.\n";
  return NULL;
}





void SetPadYMin(double yMin) {
  TH1* h = GetFirstHistInPad();
  if (h) 
    h->SetMinimum(yMin);
}


void SetPadYMax(double yMax) {
  TH1* h = GetFirstHistInPad();
  if (h) 
    h->SetMaximum(yMax);
}


void SetPadYMax() {
  TList *list = gPad->GetListOfPrimitives();
  TIter iter(list);

  TObject *obj;
  TH1 *h;

  int nhist = 0;
  double yMax = 0.;
  double yMaxCurrent = 0.;
  double yMinCurrent = 0.;

  while ( (obj = iter.Next()) ) {
    h = dynamic_cast<TH1*> (obj);
    if (h) {
      if (nhist==0) {
	// assume first is current
	yMaxCurrent = h->GetMaximum();  // this is y-range of plot 
	yMinCurrent = h->GetMinimum();
	yMax = yMaxCurrent;
      }
      else 
      if (yMax < h->GetMaximum()) {
	yMax = h->GetMaximum();  // this is max y content of other hisograms
      }
      ++nhist;
    }
  }

  if (nhist==0) {
    cout << "Did not find maximum to set yMax... no histogram found in pad.\n";
    return;
  }

  double yMaxNew;

  if (gPad->GetLogy()) {
    if (yMinCurrent > 0.) {
      double log10_yDiff = log10(yMax)-log10(yMinCurrent);
      yMaxNew = yMinCurrent*pow(10., log10_yDiff*1.05);
    }
    else { // have to get ymin from pad itself 
      double log10_yDiff = log10(yMax) - gPad->GetUymin();
      yMaxNew = pow(10., gPad->GetUymin()+log10_yDiff*1.05);
    }
  } else {
    double yDiff = yMax - yMinCurrent;
    yMaxNew = yMinCurrent+yDiff*1.05;
  }


  if (yMaxNew > yMaxCurrent) {
    SetPadYMax(yMaxNew);
  }

  return;
}




void SetPadRebin(Int_t ngroup = 2, bool resetYMax = true) {
  TList *list = gPad->GetListOfPrimitives();
  TIter iter(list);

  TObject *obj;
  TH1 *h;

  while ( (obj = iter.Next()) ) {
    h = dynamic_cast<TH1*> (obj);
    if (h) {
      h->Rebin(ngroup);
    }
  }

  if (resetYMax) {
    SetPadYMax();
  }
}




void SetPadTitle(const char *title) { // N.B. const needed to cast TStrings
  TH1* h = GetFirstHistInPad();
  if (h) 
    h->SetTitle(title);
}



#endif // ROOTDEV_SETPAD_H_
