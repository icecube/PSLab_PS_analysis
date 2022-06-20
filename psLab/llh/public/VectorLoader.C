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

#include "llh/public/classes.h"
#include "llh/public/CoordEquatorialDeg.h"
#include "llh/public/Time.h"
#include "llh/public/I3Event.h"
#include "llh/public/LlhFunctionsBase.h"


#pragma link C++ class vector<EquatorialDeg>;
#pragma link C++ class vector<I3Event>;
#pragma link C++ class vector<AnalysisSet*>;
#pragma link C++ class vector<MinuitParDef>;
#pragma link C++ class vector<Time>;

#pragma link C++ class vector< vector<double> >;
