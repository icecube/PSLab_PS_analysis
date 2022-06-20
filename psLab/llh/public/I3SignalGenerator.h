#ifndef LLH_I3SIGNALGENERATOR_H_
#define LLH_I3SIGNALGENERATOR_H_

#include <vector>

#include "rootExt/public/log_report.h"
#include "TH1D.h"
#include "TH2D.h"

//#include "fluxus/public/FluxFunction.h"
// Forward Declarations (when feasible, more efficient than including headers)
class FluxBase;

#include "llh/public/I3Event.h"
#include "llh/public/TimePdf.h"


//
// ABSTRACT BASE CLASS
//



class I3SignalGenerator : public SourceModule {
 public:
  virtual ~I3SignalGenerator() { }

  virtual I3SignalGenerator* Clone() const = 0;
  // This allows us to copy a derived class object w/o knowing what it is.
  // THIS MUST BE DEFINED SEPARATELY for each derived class.
  // For simple class (not requiring a deep copy) this should suffice:
  // { return new DerivedClass(*this); }

  // "SET" FUNCTIONS

  virtual void SetLivetime(double livetime) = 0;

  virtual void SetRecoLogEproxyDistribution(vector< vector<TH1D*> >& histo) = 0;
  virtual void SetRecoAngErrDistribution(vector< vector<TH1D*> >& histo) = 0;
  virtual void SetPSFDistribution(vector< vector<TH1D*> >& histo) = 0;
  virtual void SetEffectiveAreaDistribution(TH2D* histo) = 0;
  virtual void SaveDetectorResponseHisto(char* outfile) = 0;

  // "GET" FUNCTIONS

  virtual double GetLivetime() const = 0;

  virtual double GetMeanSrcNev() const = 0;
  // mean number of signal events generated for flux model "installed"

  virtual double GetMeanSrcNevForFluxModel(const FluxBase& fluxModel) const = 0;
  // mean number of signal events... for any other model specified

  // how much to scale model flux up or down to match Nev
  virtual double GetFluxScaleForNev(double Nev) const {
    double meanNev = GetMeanSrcNev();
    if (meanNev>0) { return Nev / meanNev; }
    return 0.;
  }

  virtual double GetFluenceNormalization() = 0;

  virtual void SetTimePdf(TimePdf * timePdf) {
    if ( timePdf->GetLivetime() ) { }
    cout << "setTimePdf not implemented!\n";
  }
  virtual TimePdf * GetTimePdf() const = 0;
  virtual void ResetTimePdf() = 0;
  virtual vector<TimePdf*> GetTimePdfVect() const = 0;

  // GENERATOR

  virtual I3Event GenerateEvent() = 0;
};



//
// MULTI-SIGNAL GENERATOR: 
//
// Combine different I3SignalGenerator objects into one signal.
//
//
// You can add as many I3SignalGenerator's of any type as you like.
// When you GenerateEvent() with Multi, it will randomly choose between 
// the different sources, based on their relative weighting.

// To begin with, their weighting is based on their event normalizations: 
// GetMeanSrcNev().  This depends only on how the source itself was 
// configured.  There are two ways this can be altered in Multi:

// 1) change the enhanceFactor for the source in Multi, which will 
// multiply the GetMeanSrcNev normalization.

// 2) the more obscure way is to change the livetime of one of the sources.  
// Each source is added with the livetime it started with, but you can 
// change it later with AccessSignalPtr(i)->SetLivetime(livetime).  (You
// can in fact get complete access to the source with this method.)

// Note that the function I3MultiSignalGenerator::SetLivetime(livetime) will 
// reset all signals with the same global livetime.


class I3MultiSignalGenerator : public I3SignalGenerator {
 private:
  double livetime_;
  vector<I3SignalGenerator*> signalPtrVect_;
  vector<double> enhanceFactorVect_;
  TimePdf * timePdf_; // with blocks and such needing to be smart about 
                      // normalizations, it sounds like a good idea to keep
                      // an extra copy here.
  int nSources_; //Number of sources, = size of signalPtrVect_

 public:
  I3MultiSignalGenerator() : livetime_(0.), timePdf_(NULL), nSources_(0) { }

  virtual ~I3MultiSignalGenerator() {
    for (unsigned int i = 0; i< signalPtrVect_.size(); ++i) {
      delete signalPtrVect_[i];
    }
  }

  // Need a deep copy, since new copies must be made for signalPtr's
  virtual I3SignalGenerator* Clone() const;

  // "SET" FUNCTIONS

  // This sets all livetimes simultaneously;
  // otherwise they can be independent of each other (e.g. a burst?)
  virtual void SetLivetime(double livetime);
  
  void SetTimePdf(TimePdf * tPdf);

  void AddSignal(const I3SignalGenerator& signal, 
		    double enhanceFactor = 1.);

  void SetEnhanceFactor(unsigned int i, double factor) {
    if (i<enhanceFactorVect_.size()) { enhanceFactorVect_[i] = factor; }
    else { log_error("%d enhanceFactor not assigned.\n",i); }
  }

  virtual void SetRecoLogEproxyDistribution(vector< vector<TH1D*> >& histo){
    for(unsigned int i = 0; i< signalPtrVect_.size(); ++i) {
      signalPtrVect_[i]->SetRecoLogEproxyDistribution(histo);
    }
  }

  virtual void SetRecoAngErrDistribution(vector< vector<TH1D*> >& histo){
    for(unsigned int i = 0; i< signalPtrVect_.size(); ++i) {
      signalPtrVect_[i]->SetRecoAngErrDistribution(histo);
    }
  }

  virtual void SetPSFDistribution(vector< vector<TH1D*> >& histo){
    for(unsigned int i = 0; i< signalPtrVect_.size(); ++i) {
      signalPtrVect_[i]->SetPSFDistribution(histo);
    }
  }

  virtual void SetEffectiveAreaDistribution(TH2D* histo){
    for(unsigned int i = 0; i< signalPtrVect_.size(); ++i) {
      signalPtrVect_[i]->SetEffectiveAreaDistribution(histo);
    }
  }

  virtual void SaveDetectorResponseHisto(char* outfile){
    signalPtrVect_[0]->SaveDetectorResponseHisto(outfile);
  }

  // "GET" FUNCTIONS

  virtual double GetLivetime() const {return livetime_;}
  
  TimePdf * GetTimePdf() const {return timePdf_;}
  virtual void ResetTimePdf() {};
  vector<TimePdf*> GetTimePdfVect() const {return vector<TimePdf*>();}

  virtual double GetMeanSrcNev() const;
  double GetMeanSrcNev(int isig);
  // adds up enhanceFactorVect_[i] * signalVect_[i].GetMeanSrcNev()

  virtual double GetMeanSrcNevForFluxModel(const FluxBase& fluxModel) const;
    // same as above, for any model
  
  virtual double GetMeanSrcNevBase() const;
  // adds up signalVect_[i].GetMeanSrcNev()

  virtual double GetMeanSrcNevForFluxModelBase(const FluxBase& fluxModel) const;
    // same as above, for any model

  virtual double GetFluenceNormalization();
  double GetFluenceNormalization(int isig);

  I3SignalGenerator* AccessSignalPtr(unsigned int i) {
    if (i<signalPtrVect_.size()) { return signalPtrVect_[i]; }
    else { log_warn("%d signalPtr not assigned.\n",i); return NULL; }
  }

  double GetEnhanceFactor(unsigned int i) {
    if (i<enhanceFactorVect_.size()) { return enhanceFactorVect_[i]; }
    else { log_warn("%d enhanceFactor not assigned.\n",i); return 0.; }
  }

  int GetNSources() { return nSources_; }

  // GENERATOR
  I3Event GenerateEvent();
};
  

// 
// SIMPLE GENERATOR FOR A POINT SOURCE
//


class I3PointGenerator : public I3SignalGenerator {
 private:
  double livetime_;
  EquatorialDeg sourceCoord_;
  bool isTimeDep_;
  bool injectEnergyRange_;
  double Emin_, Emax_;
  double angErrSrc_;
  bool extendedSrc_;

  int timeAzBins_; // to make sure we're injecting signal
                   // with the right local coordinate structure
    
  TimePdf * sourceTimePdf_;
  vector<TimePdf*> sourceTimePdfVect_;
   
  double fluenceSum_;

  vector< vector<TH1D*> > h_recoLogEproxy_;
  vector< vector<TH1D*> > h_recoAngErr_;
  vector< vector<TH1D*> > h_PSF_;
  TH2D* h_Aeff_;
  const FluxBase *fluxModel_;
  double nsFluenceRatio_;

 public:

  I3PointGenerator();
  I3PointGenerator(const FluxBase& fluxModel,
		   const EquatorialDeg& sourceCoord, 
		   double livetime,
                   TH2D* Aeff);
  I3PointGenerator(const FluxBase& fluxModel,
                   const EquatorialDeg& sourceCoord,
                   double livetime,
                   TimePdf * tPdf);
  I3PointGenerator(const FluxBase& fluxModel,
		   const EquatorialDeg& sourceCoord, 
		   double livetime,
		   TimePdf * tPdf,
                   TH2D* Aeff);
  I3PointGenerator(const FluxBase& fluxModel,
                   const EquatorialDeg& sourceCoord,
                   double livetime,
                   TimePdf * tPdf,
                   double emin,
                   double emax,
                   TH2D* Aeff); 
  I3PointGenerator(const FluxBase& fluxModel,
                   const EquatorialDeg& sourceCoord,
                   double livetime,
                   vector<TimePdf*> tPdfVect,
                   TH2D* Aeff); 
  ~I3PointGenerator() {}

  I3SignalGenerator* Clone() const { return new I3PointGenerator(*this); }
  // This allows us to copy a derived class object w/o knowing what it is.
  // THIS MUST BE DEFINED SEPARATELY for each derived class, thusly:
  // { return new DerivedClass(*this); }

  // "SET" FUNCTIONS

  // If you set events, you have to set weights too...
  void SetLivetime(double livetime) { livetime_ = livetime; }
 
  void SetTimePdf(TimePdf * tPdf) {
    sourceTimePdf_ = tPdf->Clone();
    isTimeDep_ = true;
  }

  virtual void ResetTimePdf() { delete sourceTimePdf_; sourceTimePdf_ = NULL;}

  void SetTimePdfVect(vector<TimePdf*> tPdfVect) {
    sourceTimePdfVect_.clear();
    for(unsigned int i=0; i<tPdfVect.size(); ++i){
      sourceTimePdfVect_.push_back(tPdfVect[i]->Clone());
      sourceTimePdfVect_[i]->fillHisto(1e3);
    }
  }

  void SetTimePdfVect(vector<TimePdf*> tPdfVect, float tmin, float tmax) {
    sourceTimePdfVect_.clear();
    for(unsigned int i=0; i<tPdfVect.size(); ++i){
      sourceTimePdfVect_.push_back(tPdfVect[i]->Clone());
      sourceTimePdfVect_[i]->CheckTimeBounds(tmin, tmax);
    }
  }
 
  void SetTimeAzBins(int nbins) { timeAzBins_ = nbins; }
 
  void SetSourceCoord(const EquatorialDeg& srcCoords){ sourceCoord_ = srcCoords; }

  void SetEnergyRange(double emin, double emax){
    injectEnergyRange_=true;
    Emin_ = emin;
    Emax_ = emax;
  }
 
  void SetSourceSigmas(double sigma){ angErrSrc_ = sigma; extendedSrc_=true;}

  void SetFluxModel(const FluxBase& flux){ fluxModel_ = &flux; }

  virtual void SetRecoLogEproxyDistribution(vector< vector<TH1D*> >& histo){
    h_recoLogEproxy_.clear();
    h_recoLogEproxy_ = histo;
  }

  virtual void SetRecoAngErrDistribution(vector< vector<TH1D*> >& histo){
    h_recoAngErr_.clear();
    h_recoAngErr_ = histo;
  }

  virtual void SetPSFDistribution(vector< vector<TH1D*> >& histo){
    h_PSF_.clear();
    h_PSF_ = histo;
  }

  virtual void SetEffectiveAreaDistribution(TH2D* histo){
    h_Aeff_ = new TH2D(*histo);
    SetNsFluenceRatio( GetFluenceNormalization() );
  }
 
  virtual void SaveDetectorResponseHisto(char* outfile);

  void SetNsFluenceRatio(double ratio){ nsFluenceRatio_ = ratio; }

  // "GET" FUNCTIONS
  
  virtual double GetLivetime() const {return livetime_;}
  virtual double GetMeanSrcNev() const {
    if(!isTimeDep_) { return nsFluenceRatio_ * livetime_; }
    else { return nsFluenceRatio_; }
  }

  virtual double GetMeanSrcNevForFluxModel(const FluxBase& fluxModel) const;
  TimePdf * GetTimePdf() const {return sourceTimePdf_;}
  vector<TimePdf*> GetTimePdfVect() const {return sourceTimePdfVect_;}
  double GetFluenceSum() {return fluenceSum_;}
  
  EquatorialDeg GetSourceCoord() const {return sourceCoord_;}

  vector< vector<TH1D*> > GetPSFDistribution(){        return h_PSF_;           }
  vector< vector<TH1D*> > GetRecoAngErrDistribution(){ return h_recoAngErr_;    }
  vector< vector<TH1D*> > GetRecoEnergyDistribution(){ return h_recoLogEproxy_; }
  TH2D* GetEffectiveAreaDistribution(){               return h_Aeff_;          }
  virtual double GetFluenceNormalization();

  // GENERATOR

  I3Event GenerateEvent();
  bool CheckTimeInGRL(double mjd);
  TH1D* GenerateEnuPDF();
};


#endif // LLH_I3SIGNALGENERATOR_H_
