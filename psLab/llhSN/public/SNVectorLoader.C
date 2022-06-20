// VectorLoader.C
//
// Root can compile vectors without much trouble, but the CINT interpreter 
// has trouble.  Solution:  .L VectorLoader.C
//
// Just add here the functionality you require on the command line 
// or in a macro.

// NOTE: currently I have not been able to get 
// the iterator operations == or != to work, which makes the iterator 
// pretty useless..


#include <vector>

#include "llhSN/public/SNEvent.h"
//#include "$LAB_MAIN_DIR/llhSN/public/SNEvent.h"

#pragma link C++ class vector<SNEvent>;

