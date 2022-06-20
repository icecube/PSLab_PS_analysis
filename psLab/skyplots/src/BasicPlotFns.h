#ifndef BASICPLOTFNS_H_
#define BASICPLOTFNS_H_


#include "TH1.h"
#include "TList.h"
#include "TPaletteAxis.h"

void log_warn(char* c) { printf("%s",c); }



TPaletteAxis* GetHistogramPalette(TH1* h) {
  TPaletteAxis* paxis = dynamic_cast<TPaletteAxis*>
    ( h->GetListOfFunctions()->FindObject("palette") );
  if (paxis == NULL) {
    log_warn("TPaletteAxis* is null.  You need to draw the histogram or\n"
	     "  update the pad first.");
  }
  return paxis;
}


#endif // BASICPLOTFNS_H_
