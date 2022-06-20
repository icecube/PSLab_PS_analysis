#ifndef LLH_NEWLLHBOXSTACK_H_
#define LLH_NEWLLHBOXSTACK_H_

#include "TMinuit.h"
#include "TStopwatch.h" // This may be temporary, for optimization purposes

#include "llh/public/classes.h"
#include "llh/public/MinuitAnalysisFn.h"
#include "llh/public/MultiAnalysisSet.h"

#include "llhTimeDep/public/NewLlhBoxTime.h"
#include "llh/public/CoordEquatorialDeg.h"

// Forward Declarations (when feasible, more efficient than including headers)
//class EnergyProb;

class NewLlhBoxStack : public NewLlhBoxTime {

 public:
  // New for stacking/extended sources:
  vector<vector<double> > srcWeightsTable_;
  vector<EquatorialDeg> srcCoords_;
  vector<double> srcSigmas_;
  vector<vector<double> > spaceRatioVects_;
  vector<vector<double> > timeRatioVects_;
  vector<double> tminSrc_;
  vector<double> tmaxSrc_;
  double tminPer_;
  double tmaxPer_;

  void OptimizeEventSelection();

  NewLlhBoxStack();
  virtual ~NewLlhBoxStack() { }

  void PrepareAnalysis();

  // These are exact copies from NewLlhEnergy, but are required because fcn_ is of
  //  a different type (overloaded).  
  virtual double EvalFCN(const vector<double>& parVect) const;

  double EvaluateLlh(double nSrc, double gamma);

  virtual double EvaluateLlh(double* parArray) {
    return EvaluateLlh(parArray[0], parArray[1]);
  }

  // New for stacking:
  /* @brief Copies external vector of source EquatorialDeg coordinates to 
   / interal member
  */
  void SetSourceCoords(vector<EquatorialDeg>& sourceCoords);

  /* @brief Copies external vector of source sigmas to interal member
  */
  void SetSourceSigmas(vector<double>& sourceSigmas);

  void SetSourceTminTmax(vector<double> minTime, vector<double> maxTime);
  void SetPeriodTminTmax(double tmin, double tmax);

  /* @brief A routine to set the Source Weights Table, the weight of each source
   /  depending on source location and gamma.  This table is first calculated
   /  by something like macro_StackingWeights.C -JD
  */
  void SetStackedWeightTable(vector<vector<double> > srcWeightsArray);

  vector< vector<double> > GetSpaceRatios() { return spaceRatioVects_;}

};

NewLlhBoxStack *tdBoxStack;
#endif // LLH_NEWLLHBOXSTACK_H_
