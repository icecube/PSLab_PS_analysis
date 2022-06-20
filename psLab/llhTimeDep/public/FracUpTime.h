#ifndef LLHTIMEDEP_FRACUPTIME_H_
#define LLHTIMEDEP_FRACUPTIME_H_

#include <string>

#include "TH1D.h"


// This takes in a histogram of the uptimes and limits you're interested in,
// and returns the fraction of uptime for things like fluence scaling.

// First off, uptime between some min to max times

double FracUpTime(TH1D* a, double xmin, double xmax);


// This does the same thing but is specifically for a block function.
// farther down it does the same for a Gaussian.

double BlockFracUpTime(TH1D* a, string fname, double thresh);

double FracUpTimeGaus(TH1D* a, double mean, double sigma);


#endif // LLHTIMEDEP_FRACUPTIME_H_
