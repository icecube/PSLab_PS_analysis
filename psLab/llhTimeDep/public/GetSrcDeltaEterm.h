#ifndef LLHTIMEDEP_GETSRCDELTAETERM_H_
#define LLHTIMEDEP_GETSRCDELTAETERM_H_

#include <vector>

#include "llh/public/I3Event.h"

double GetSrcDeltaEterm(vector<I3Event>* evVectPtr, 
			double pcmin=0.05, double pcmax=0.95);


#endif // LLHTIMEDEP_GETSRCDELTAETERM_H_
