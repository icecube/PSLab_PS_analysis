#ifndef LLH_MULTIGAUSSANALYSISFN_H_
#define LLH_MULTIGAUSSANALYSISFN_H_

#include "llh/public/MinuitAnalysisFn.h"
#include "llhTimeDep/public/NewLlhGausTime.h"
#include "llh/public/LlhFunctionsBase.h"
#include "rootExt/public/generalfunctions.h"

//#include "llh/public/I3Event.h"

#include "TMinuit.h"

class FluxBase;
class MultiGaussAnalysisFCN;

class MultiGaussAnalysisFn : public AnalysisFn {
    protected:
      TMinuit* minuit4_;
      TMinuit* minuit6_;
  
      EquatorialDeg* srcCoordGuess_;
      double srcCoordUncRA_;
      double srcCoordUncDEC_;
      double stepSizeSrcCoord_;

      MultiGaussAnalysisFCN* fcn_; 

      vector<NewLlhGausTime*> analysisFnVect_; // this needs to be specifically
                                              // for periodic analyses, not an <I3Analysis>
      const NewLlhGausTime_ParTranslator* parTrans_;
      vector<MinuitParDef> parDefVect_;
      double tmin_;
      double tmax_;

      double nEventsTot_;

      int nFlareGuess_;
      double TSthr_;
      vector<double> nsrcGuess_;
      vector<double> gammaGuess_;
      vector<double> meanGuess_;
      vector<double> sigmaGuess_;
      double sigmamin_;
      double sigmamax_;
      double sigmaFixed_;
      double meanFixed_;
      bool   isSetLowLimSigma_;
      bool   isSetUpLimSigma_;

      double gammaMin_;
      double gammaMax_;
      double gammaFixed_;

      double srcFracMax_;
      double nSigmaTrunc_;

      double logLambdaBest_;
      void StoreLogLambdaBest(TMinuit *minuit);

      vector<pair<double, double> > nSrcBest_;
      vector<pair<double, double> > gammaBest_;
      vector<pair<double, double> > meanBest_;
      vector<pair<double, double> > sigmaBest_;
      double nsSet_[5];

      int  minuitOut_;   //store minimization result (0=fit converged, <0=fit not converged)                                    

      bool histoForProb_;
      bool fitSrcCoord_;
      bool doInitGauss_;   //flag to use the class method to initialize parameters                                  
      bool optUsePrior_;   //flag to choose between standard/bayesian likelihood                
	                     //with a prior term modelled with a 2D (a)symmetric Gaussian PDF    
      bool JimsTerm_;
      bool optStoreRatios_;
      const Coord* priorCoord_; //coord  of the event used to give a prior term                                                 
      double sigmaX_;   //X extension of the area associated to prior term            
      double sigmaY_;   //Y extension of the area associated to prior term              
      double theta_;    //rotation angle of area associated to prior term                                   
	                  //(in degrees, from X axis counterclockwise)                                                
      vector<const EnergyProb*> eProbVect_;
      vector<double> spaceRatioVect_;
      vector<double> bkgSpaceProbVect_;
      vector<double> lcBkgProbVect_;
      vector<double> sigmaSpaceVect_;
      vector<double> spaceWeightVect_;
      vector<double> enWeightVect_;
      vector<double> timeWeightVect_;
      vector<double> raVect_;
      vector<double> decVect_;
      vector<double> angErrVect_;
      vector<double> eneVect_;
      vector<double> timeVect_;
      vector<int> eventID_;
      vector<int> runID_;
 
      bool optParAuto_[6]; //flags to keep the llh params free (true) to vary or fixed (false) during minimization 
      bool isSingleFlare_;

      TH1 * pvalHisto_;   // distribution to estimate p-values
      TH1 * nullTestStat_; // a few parameters for using a custom teststat

    public:
      double seedWtMin;
      int nPar;

      MultiGaussAnalysisFn();
      virtual ~MultiGaussAnalysisFn() {
      }

      void ClearAnalysisFn(){
        analysisFnVect_.clear();
      }

      virtual void AddAnalysisFn(AnalysisFn* llh) {
          NewLlhGausTime* llh1 = dynamic_cast<NewLlhGausTime*>(llh);
          analysisFnVect_.push_back(llh1);
      }
        
      virtual void SetAnalysisSet(AnalysisSet*) {log_error("use AddAnalysis  instead of  SetAnalysisSet(aSet)\n");}

      virtual void SetSearchCoord(const Coord& coord) {
          srcCoord_ = &coord;
          for (int i=0; i<int(analysisFnVect_.size()); ++i) {
              analysisFnVect_[i]->SetSearchCoord(coord);
          }
      }
      double calcSpatialPrior(double testRA, double testDEC);
      virtual void GetFlareGuessGauss(vector<double> & Guess_nsrc, vector<double> & Guess_gamma, vector<double> & Guess_mean, vector<double> & Guess_rms, double & sigmamin);
      bool GetJimsTerm() {return JimsTerm_;}
      int GetMinimizationResult() {return minuitOut_;}

      virtual void SetParDefs(vector<MinuitParDef>& parDefVect) { parDefVect_ = parDefVect;}
      virtual void AddParDef(MinuitParDef parDef){
          parDefVect_.push_back(parDef);
      }
      virtual void SetParTranslator(const NewLlhGausTime_ParTranslator* pt) { parTrans_ = pt; }

      void SetUseFitSrc(bool flag) {//to activate the fit of the src position in the llh                    
          fitSrcCoord_ = flag;
          optParAuto_[4] = flag;
          optParAuto_[5] = flag;
          for (int i=0; i<int(analysisFnVect_.size()); ++i) {
            analysisFnVect_[i]->SetUseFitSrc(flag);
          }
        }
        void SetFitSrc(EquatorialDeg& srcCoord, double uncRA, double uncDEC,
		       double stepSize) {
          //to activate the fit of the src position in the llh                                              
          //and set their init values, ranges and minuit step                                                
          fitSrcCoord_   = true;
          optParAuto_[4] = true;
          optParAuto_[5] = true;
          srcCoordGuess_ = &srcCoord;
          srcCoordUncRA_    = uncRA;
          srcCoordUncDEC_   = uncDEC;
          stepSizeSrcCoord_ = stepSize;
          for (int i=0; i<int(analysisFnVect_.size()); ++i) {
            analysisFnVect_[i]->SetFitSrc(srcCoord,uncRA,uncDEC,stepSize);
          }
        }
        void SetParamGuess(vector<double> nsGuess, vector<double> gammaGuess,
              vector<double> timeGuess, vector<double> sigmaGuess,
              bool optParAuto) {
          //in case the initialization of the (4) llh parameters is done manually                                     
          doInitGauss_= false;
          nsrcGuess_  = nsGuess;
          gammaGuess_ = gammaGuess;
          meanGuess_  = timeGuess;
          size_t nel = sigmaGuess.size();
          for(size_t i=0; i<nel; i++) sigmaGuess_[i] = pow(10,sigmaGuess[i]);
	      //optParAuto_[i] = true  : param. i free to vary in the fit minimization                          
	      //optParAuto_[i] = false : param. i kept fixed in the fit minimization                            
	      optParAuto_[0] = optParAuto;
	      optParAuto_[1] = optParAuto;
	      optParAuto_[2] = optParAuto;
	      optParAuto_[3] = optParAuto;
	    }
	  void SetParamGuess(vector<double> nsGuess, vector<double> gammaGuess,
               vector<double> timeGuess, vector<double> sigmaGuess,
               double srcRAGuess, double srcDECGuess,
               bool optParAuto) {
	      //in case the initialization of the (6) llh parameters is done manually                  
              doInitGauss_= false;
              nsrcGuess_  = nsGuess;
              gammaGuess_ = gammaGuess;
	      meanGuess_  = timeGuess;
	      size_t nel = sigmaGuess.size();
              for(size_t i=0; i<nel; i++) sigmaGuess_[i] = pow(10,sigmaGuess[i]);

	      optParAuto_[0] = optParAuto;
	      optParAuto_[1] = optParAuto;
	      optParAuto_[2] = optParAuto;
	      optParAuto_[3] = optParAuto;
	      SetUseFitSrc(true);
	      srcCoordGuess_->SetCoords(srcRAGuess,srcDECGuess);
	  }
	  void SetOptParAuto(bool flag) {
	      optParAuto_[0] = flag;
	      optParAuto_[1] = flag;
	      optParAuto_[2] = flag;
	      optParAuto_[3] = flag;
              for (int i=0; i<int(analysisFnVect_.size()); ++i) {
                analysisFnVect_[i]->SetOptParAuto(flag);
              }
	  }
          void SetOptParAuto(int i, bool flag) {
              if(i>=0 && i<4){
                  optParAuto_[i] = flag; 
           
                  for (int j=0; j<int(analysisFnVect_.size()); ++j) { 
                    analysisFnVect_[j]->SetOptParAuto(i, flag); 
                  }
              }
              else Printf("Cannot fix parameter %d: choose an index between 0 and 3", i);
          }
	  void SetSpatialPrior(const Coord& priorCoord, double sigmaX, double sigmaY, double theta/*in degree*/) {
	      optUsePrior_ = true;
	      priorCoord_  = &priorCoord;
	      sigmaX_ = sigmaX;
	      sigmaY_ = sigmaY;
	      theta_  = theta;
              for (int i=0; i<int(analysisFnVect_.size()); ++i) {
	        analysisFnVect_[i]->SetSpatialPrior(priorCoord,sigmaX,sigmaY,theta);
              }
	  }
          void SetOptStoreRatios(bool flag) {optStoreRatios_ = flag;}
	  void UnsetParamGuess() {
	      doInitGauss_= true;
              for (int i=0; i<int(analysisFnVect_.size()); ++i) {
                analysisFnVect_[i]->UnsetParamGuess();
              }
	    }

          void SetNFlareGuess(int nflares){nFlareGuess_ = nflares;}
    
          void SetGammaFixed(double g) {gammaFixed_ = g;}
          void SetMeanFixed(double m) {meanFixed_ = m;}
          void SetSigmaFixed(double s) {sigmaFixed_ = s;}

	  bool IsFitSrc() {return fitSrcCoord_;}

          virtual void PrepareAnalysis() {
            nEventsTot_ = 0.;
            for (int i=0; i<int(analysisFnVect_.size()); ++i) {
                analysisFnVect_[i]->PrepareAnalysis();
                nEventsTot_ += analysisFnVect_[i]->Get_nEvents();
                cout<<"sample "<<i<<" has "<<analysisFnVect_[i]->Get_nEvents()<<" events, " << analysisFnVect_[i]->selectedList_.GetSize() << " are selected by PrepareAnalysis."<<endl;
            }
          }
                
        virtual void MaximizeLlh();
        
        void ClearBestFitParams() {
            nSrcBest_.clear();
            gammaBest_.clear();
            meanBest_.clear();
            sigmaBest_.clear();
        }

        vector<I3Event> GetAllEvents();
        virtual double EvaluateLlh(double *parValueArray);
        double EvaluateLlh( double a, double b, double c, double d) {    
            double dd[] = {a, b, c, log10(d) }; //use log10(sigma) in minimizer
            double Llh = EvaluateLlh(dd);
            return Llh;
        }
        vector<pair<double,double> > Get_nSrcBest() const {return nSrcBest_;}
        vector<pair<double,double> > Get_gammaBest() const {return gammaBest_;}
        vector<pair<double,double> > Get_meanBest() const {return meanBest_;}
        vector<pair<double,double> > Get_sigmaBest() const {return sigmaBest_;}

        vector<double> GetNsrcGuess(){ return nsrcGuess_;}
        vector<double> GetMeanGuess() { return meanGuess_; }
        vector<double> GetSigmaGuess(){ return sigmaGuess_;}

        double GetNsSet(int iset) {return nsSet_[iset];}
        double Get_nEvents() const {return nEventsTot_;}
        int GetNFlareGuess() {return nFlareGuess_;}
        double GetTmin() {return tmin_;}
        double GetTmax() {return tmax_;}

        virtual double EvalFCN(const vector<double>& parVect) const;
        
        double GetProbFromHisto(double teststat) const;
        void SetNullTestStat(TH1D * inputhisto);

	bool GetUsePrior() const { return optUsePrior_; }
        bool GetOptStoreRatios() { return optStoreRatios_; }
  
        virtual double GetPar(int i) const { 
          double par=0., err=0.;
          minuit4_->GetParameter(i, par, err);
          return par;
        }
        
        virtual double Get_logLambdaBest() const { return logLambdaBest_; }
        virtual double GetTestStatistic() const { return Get_logLambdaBest(); }
        virtual double GetSigmaMin() { return sigmamin_;}
        virtual double GetEstProb() const {
            if (histoForProb_)  return GetProbFromHisto( Get_logLambdaBest() ); 
            double chiSq = 2. * Get_logLambdaBest();
            if(chiSq<0)  chiSq=0.; 
            double p_temp, p;
            int nDoF = analysisFnVect_[0]->ndof_;
            chisq_prob(chiSq, nDoF, &p_temp, &p);
            return p / 2.;  // one-sided chi-sq prob
        }
 
        const NewLlhGausTime_ParTranslator* GetParTranslator() { return parTrans_; }
        vector<NewLlhGausTime*> GetAnalysisFn()     { return analysisFnVect_;}
 
        void SetTimeBounds(double tmin, double tmax) {
            tmin_ = tmin;
            tmax_ = tmax;
            for (int i=0; i<int(analysisFnVect_.size()); ++i) {
              analysisFnVect_[i]->SetTimeBounds(tmin,tmax);
            }
         }

        void SetLowLimitSigma(double limit) {
          sigmamin_ = limit;
          isSetLowLimSigma_ = true;
          for (int i=0; i<int(analysisFnVect_.size()); ++i) {
            analysisFnVect_[i]->SetLowLimitSigma(limit);
          }
        }
        void SetUpLimitSigma(double limit) {
          sigmamax_ = limit;
          isSetUpLimSigma_ = true;
          for (int i=0; i<int(analysisFnVect_.size()); ++i) {
            analysisFnVect_[i]->SetUpLimitSigma(limit);
          }
        }

        void SetIsSingleFlare(bool b){ isSingleFlare_ = b; }

        vector<double> GetSpatialWeights(){return spaceWeightVect_;}
        vector<double> GetEnergyWeights() {return enWeightVect_;}
        vector<double> GetTimeWeights()   {return timeWeightVect_;}
        vector<double> GetraVect()        {return raVect_;}
        vector<double> GetdecVect()       {return decVect_;}
        vector<double> GetAngErrVect()    {return angErrVect_;}
        vector<double> GettimeVect()      {return timeVect_;}
        vector<double> GetEneVect()       {return eneVect_;}
        vector<int> GetEventID()          {return eventID_;}
        vector<int> GetRunID()          {return runID_;}

        void SetNSigmaTrunc(double nSig){
          nSigmaTrunc_=nSig;
          for (int i=0; i<int(analysisFnVect_.size()); ++i) {
            analysisFnVect_[i]->SetNSigmaTrunc(nSig);
          }
        }
        void SetTSthr(double TSthr) {TSthr_= TSthr;}
        void SetSeedWtMin(double value) {
            seedWtMin = value;
	    for (int i=0; i<int(analysisFnVect_.size()); ++i) {
                analysisFnVect_[i]->SetSeedWtMin(value);
            }
        }
};

#endif // LLH_MultiGAUSSANALYSISFN_H_
