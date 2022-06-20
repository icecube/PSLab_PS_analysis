#ifndef ROOTEXT_ROOT_DRAW_TOOLS_H_
#define ROOTEXT_ROOT_DRAW_TOOLS_H_

#include "TLine.h"

// Create vertical line ranging top to bottom (ok for lin or log Pad)
// Example Usage: vline(5.)->Draw();

TLine* vline(Double_t x, Color_t lcolor=1, Width_t lwidth=1, Style_t lstyle=1);

// Create horizontal line ranging top to bottom (ok for lin or log Pad)
// Example Usage: vline(5.)->Draw();

TLine* hline(Double_t y, Color_t lcolor=1, Width_t lwidth=1, Style_t lstyle=1);


#endif // ROOTEXT_ROOT_DRAW_TOOLS_H_

