#ifndef LLH_MULTIBOXANALYSISFNSTACK_H_
#define LLH_MULTIBOXANALYSISFNSTACK_H_

#include "llh/public/MinuitAnalysisFn.h"
#include "llhTimeDep/public/NewLlhBoxTime.h"
#include "llhTimeDep/public/NewLlhBoxTimeStack.h"
#include "llh/public/LlhFunctionsBase.h"
#include "rootExt/public/generalfunctions.h"

//#include "llh/public/I3Event.h"

#include "TMinuit.h"

class FluxBase;

class MultiBoxAnalysisFnStack : public AnalysisFn {
    public: 

      vector<NewLlhBoxStack*> analysisFnVect_; // this needs to be specifically
                                              // for periodic analyses, not an <I3Analysis>
      const NewLlhBoxTime_ParTranslator* parTrans_;
      vector<MinuitParDef> parDefVect_;
      double tmin_;
      double tmax_;

      double nEventsTot_;

      double nsrcGuess_;
      double gammaGuess_;
      double boxMinGuess_;
      double boxMaxGuess_;

      double gammaMin_;
      double gammaMax_;
      double gammaFixed_;

      double boxMinFixed_;
      double boxMaxFixed_;
      double maxClusterLength;

      double srcFracMax_;

      double logLambdaBest_;
      void StoreLogLambdaBest(TMinuit *minuit);

      double nSrcBest_;
      double gammaBest_;
      double boxMinBest_;
      double boxMaxBest_;
      double nsSet_[5];

      int monitorLevel_;
      int  minuitOut_;   //store minimization result (0=fit converged, <0=fit not converged)                                    
      bool histoForProb_;
      bool JimsTerm_;
      bool optUseEnergy_;
      bool optStoreRatios_;
                                                     
      vector<const EnergyProb*> eProbVect_;
      vector<double> spaceRatioVect_;
      vector<double> bkgSpaceProbVect_;
      vector<double> lcBkgProbVect_;
      vector<double> spaceWeightVect_;
      vector<double> enWeightVect_;
      vector<double> timeWeightVect_;
      vector<double> raVect_;
      vector<double> decVect_;
      vector<double> angErrVect_;
      vector<double> eneVect_;
      vector<double> timeVect_;
      vector<int> eventID_;
 
      bool optParAuto_[6]; //flags to keep the llh params free (true) to vary or fixed (false) during minimization 

      TH1 * pvalHisto_;   // distribution to estimate p-values
      TH1 * nullTestStat_; // a few parameters for using a custom teststat

      double seedWtMin;
      int nPar;

      MultiBoxAnalysisFnStack();
      virtual ~MultiBoxAnalysisFnStack() {
      }

      void ClearAnalysisFnStack(){
        analysisFnVect_.clear();
      }

      virtual void AddAnalysisFn(AnalysisFn* llh) {
          NewLlhBoxStack* llh1 = dynamic_cast<NewLlhBoxStack*>(llh);
          dynamic_cast<I3Analysis*>(llh1->GetAnalysisSet())->SetIsStacking(true);
          analysisFnVect_.push_back(llh1);
      }
        
      virtual void SetAnalysisSet(AnalysisSet*) {log_error("use AddAnalysis  instead of  SetAnalysisSet(aSet)\n");}
 
      virtual void GetFlareGuess(bool UseE, double & Guess_nsrc, double & Guess_gamma);

      int GetMinimizationResult() {return minuitOut_;}

      virtual void SetParTranslator(const NewLlhBoxTime_ParTranslator* pt) { parTrans_ = pt; }
	  
      void SetOptStoreRatios(bool flag) {optStoreRatios_ = flag;} 
    
      void SetGammaFixed(double g) {gammaFixed_ = g; optParAuto_[1] = false;}
      void SetBoxMinFixed(double flaremin) {boxMinFixed_ = flaremin; optParAuto_[2] = false;}
      void SetBoxMaxFixed(double flaremax) {boxMaxFixed_ = flaremax; optParAuto_[3] = false;}


      virtual void PrepareAnalysis() {
        nEventsTot_ = 0.;
        for (int i=0; i<int(analysisFnVect_.size()); ++i) {
            analysisFnVect_[i]->PrepareAnalysis();
            nEventsTot_ += analysisFnVect_[i]->Get_nEvents();
            cout<<"sample "<<i<<" has "<<analysisFnVect_[i]->Get_nEvents()<<" events, " << analysisFnVect_[i]->selectedList_.GetSize() << " are selected by PrepareAnalysis."<<endl;
        }
      }
                
      virtual void MaximizeLlh(); 

      vector<I3Event> GetAllEvents();
      virtual double EvaluateLlh(double *parValueArray);
          double EvaluateLlh( double a, double b) {    
          double dd[] = {a, b }; 
          double Llh = EvaluateLlh(dd);
          return Llh;
      }
      double Get_nSrcBest() const {return nSrcBest_;}
      double Get_gammaBest() const {return gammaBest_;}
      double Get_BoxMinBest() const {return boxMinBest_;}
      double Get_BoxMaxBest() const {return boxMaxBest_;}

      double GetNsrcGuess(){ return nsrcGuess_;}
      double GetBoxMinGuess() { return boxMinGuess_; }
      double GetBoxMaxGuess(){ return boxMaxGuess_;}

      double GetNsSet(int iset) {return nsSet_[iset];}
      double Get_nEvents() const {return nEventsTot_;}

      virtual double EvalFCN(const vector<double>& parVect) const;
        
      double GetProbFromHisto(double teststat) const;
      void SetNullTestStat(TH1D * inputhisto);

      bool GetOptStoreRatios() { return optStoreRatios_; }
  
      virtual double GetPar(int i) const { 
         return i;
      }
        
      virtual double Get_logLambdaBest() const { return logLambdaBest_; }
      virtual double GetTestStatistic() const { return Get_logLambdaBest(); }
      virtual double GetEstProb() const {
            if (histoForProb_)  return GetProbFromHisto( Get_logLambdaBest() ); 
            double chiSq = 2. * Get_logLambdaBest();
            if(chiSq<0)  chiSq=0.; 
            double p_temp, p;
            int nDoF = analysisFnVect_[0]->ndof_;
            chisq_prob(chiSq, nDoF, &p_temp, &p);
            return p / 2.;  // one-sided chi-sq prob
        }
  
      void SetTimeBounds(double tmin, double tmax) {
            tmin_ = tmin;
            tmax_ = tmax;
            for (int i=0; i<int(analysisFnVect_.size()); ++i) {
              analysisFnVect_[i]->SetTimeBounds(tmin,tmax);
            }
         }

      
      vector<double> GetSpatialWeights(){return spaceWeightVect_;}
      vector<double> GetEnergyWeights() {return enWeightVect_;}
      vector<double> GetTimeWeights()   {return timeWeightVect_;}
      vector<double> GetraVect()        {return raVect_;}
      vector<double> GetdecVect()       {return decVect_;}
      vector<double> GetAngErrVect()    {return angErrVect_;}
      vector<double> GettimeVect()      {return timeVect_;}
      vector<double> GetEneVect()       {return eneVect_;}
      vector<int> GetEventID()          {return eventID_;}

      virtual void SetSearchCoord(const Coord& coord) {
          srcCoord_ = &coord;
          for (int i=0; i<int(analysisFnVect_.size()); ++i) {
              analysisFnVect_[i]->SetSearchCoord(coord);
          }
      }
      
      void SetSeedWtMin(double value) {
          seedWtMin = value;
          for (int i=0; i<int(analysisFnVect_.size()); ++i) {
                analysisFnVect_[i]->SetSeedWtMin(value);
          }
      }
};

#endif // LLH_MULTIBOXANALYSISFNSTACK_H_
