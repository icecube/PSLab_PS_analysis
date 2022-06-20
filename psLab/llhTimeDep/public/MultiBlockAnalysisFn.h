#ifndef LLH_MULTIBLOCKANALYSISFN_H_
#define LLH_MULTIBLOCKANALYSISFN_H_

//#include <exception>

#include "rootExt/public/generalfunctions.h"
#include "llh/public/LlhFunctionsBase.h"
#include "llh/public/MinuitAnalysisFn.h"
#include "llh/public/MultiAnalysisSet.h"
#include "llhTimeDep/public/BlockLevel.h"
#include "llhTimeDep/public/NewLlhBlockTime.h"

class FluxBase;
//class MultiBlockAnalysisFCN;


class MultiBlockAnalysisFn : public AnalysisFn {
 protected:
  TMinuit* minuit_;
  //MultiBlockAnalysisFCN* fcn_;

  vector<AnalysisFn*> analysisFnVect_;
  const ParTranslator* parTrans_;
  vector<MinuitParDef> parDefVect_;

  double logLambdaBest_;
  void StoreLogLambdaBest();

  int minuitOut_;

 public:

  int nPar;
  double laglimit_;
  string blocksFile_;
  int Ndof;

  MultiBlockAnalysisFn();
  virtual ~MultiBlockAnalysisFn() {
    if (minuit_) delete minuit_;
  }

  TMinuit* Minuit() { return minuit_; }

  const ParTranslator* GetParTranslator() { return parTrans_; }
  vector<AnalysisFn*> GetAnalysisFn()     { return analysisFnVect_;}

  virtual void AddAnalysisFn(NewLlhBlockTime* llh) {

    for (int j=0;j <= (int) analysisFnVect_.size(); j++) {
        if (llh->ndof_!=Ndof && Ndof!=-1){
            cout << "ERROR: the newly added llh has a different number of ndof !!! stopping." << endl;
            exit(0);
        }
    }
    Ndof=llh->ndof_;
    cout << "for calculating p values will use Ndof= " << Ndof << endl;
    analysisFnVect_.push_back((AnalysisFn*) llh);
  }

  virtual void SetAnalysisSet(AnalysisSet*) {
    log_error("use AddAnalysis  instead of  SetAnalysisSet(aSet)\n");
  }

  virtual void SetSearchCoord(const Coord& coord) {
    srcCoord_ = &coord;
    for (int i=0; i<int(analysisFnVect_.size()); ++i) {
      analysisFnVect_[i]->SetSearchCoord(coord);
    }
  }

  virtual void SetParDefs(vector<MinuitParDef>& parDefVect) {
    parDefVect_ = parDefVect;
  }

  virtual void SetParTranslator(const ParTranslator* pt) { parTrans_ = pt; }


  virtual void PrepareAnalysis() {
    for (int i=0; i<int(analysisFnVect_.size()); ++i) {
      analysisFnVect_[i]->PrepareAnalysis();
    }
  }

  virtual void MaximizeLlh();

  virtual double EvalFCN(const vector<double>& parVect) const;

  virtual double EvaluateLlh(double *parValueArray);

  virtual double GetPar(int i) const { double par=0., err=0.; return minuit_->GetParameter(i, par, err); }
  virtual double Get_logLambdaBest() const { return logLambdaBest_; }
  virtual double GetTestStatistic() const { return Get_logLambdaBest(); }
  virtual double GetEstProb() const {
    double chiSq = 2. * Get_logLambdaBest();
    double p_temp, p;
    chisq_prob(chiSq, Ndof, &p_temp, &p);
    return p / 2.;  // one-sided chi-sq prob
  }

  vector<I3Event> GetAllEvents();
  vector<I3Event> GetAllSrcEvents();
  //  void GetFlareGuess(double & Guess_nsrc, double & Guess_gamma, double & Guess_mean, double & Guess_rms);

  double SearchForLag(double laglimit);
  void SearchBlockSpace(string blocksFile, double lag, double & initV, double & maxT);

};



/*class MultiBlockAnalysisFCN : public ROOT::Minuit2::FCNBase {
 private:
  MultiBlockAnalysisFn* ptr;
 public:
  // Pure Virtual Fn inherited from ROOT::Minuit2::FCNBase
  // This is related to how "errors" are defined, see Minuit2 for documentation
  virtual double Up() const { return 0.5; }

  // Pure Virtual Fn inherited from ROOT::Minuit2.:FCNBase
  // This is what gets minimized: you have to define your likelihood
  virtual double operator() (const vector<double>& par) const;

  // Here's where the connection is made to FCN can access the data
  virtual void Point(AnalysisFn* fn) {
    ptr = dynamic_cast<MultiBlockAnalysisFn*>(fn);
  }
};*/

MultiBlockAnalysisFn *tdMLlhBlock;
#endif // LLH_MULTIBLOCKANALYSISFN_H_

