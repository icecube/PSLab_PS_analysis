#ifndef MACRO_ARK_H_
#define MACRO_ARK_H_

#include "llh/public/classes.h"
#include "llh/public/CoordEquatorialDeg.h"
#include "llh/public/SimpleEnergyProb.h"
#include "llh/public/DecBkgProb.h"
#include "fluxus/public/FluxFunction.h"
#include "llh/public/I3Analysis.h"
#include "llh/public/MultiAnalysisSet.h"
#include "llh/public/EventLoaderExt.h"
#include "llhTimeDep/public/LocalCoordBkgProb.h"
#include "TFile.h"

class Ark {
    public:
    double livetime;
    double tmin;
    double tmax;
    AnalysisSet* psData;
    EquatorialDeg mySrcLocation;

    Ark(): livetime(0) { }
    virtual ~Ark() {
        if (psData) { delete psData; };
    }  // take care of psData in non-abstract derived classes

    virtual void SetPointSource(EquatorialDeg srcCoord, const FluxBase& flux){};

    virtual void SetPointSource(EquatorialDeg srcCoord, const FluxBase& flux, TimePdf * tPdf, bool SetNormFromGRL=false) = 0;

    virtual void SetPointSource(EquatorialDeg srcCoord, const FluxBase& flux, vector<TimePdf*> tPdf, bool SetNormFromGRL=false){};

    virtual void SetPointSource(EquatorialDeg srcCoord, const FluxBase& flux, TimePdf* tPdfVect, double emin, double emax, bool SetNormFromGRL=false){};

    // just some user-friendly alternatives

    void SetPointSource(double raDeg, double decDeg, double fluxConstant, double index) {
        PowerLawFlux flux(fluxConstant, index);
        printf("FluxConstant: %lg   PowerLaw Index: %lg   (GeV^-1 cm^-2 s^-1)\n", flux.GetSpectralIndex(), flux.GetFluxConstant() );
        SetPointSource(EquatorialDeg(raDeg, decDeg), flux);
    }

    void SetPointSource(double raDeg, double decDeg, const char* formula) {
        FormulaFlux flux(formula);
        printf("Formula: %s\n  where x is GeV, and flux is GeV^-1 cm^-2 s^-1\n", flux.GetFormula() );
        SetPointSource(EquatorialDeg(raDeg, decDeg), flux);
    }
  
    void SetPointSource(double raDeg, double decDeg, const char* formula, TimePdf * tPdf, bool SetNormFromGRL=false) {
        FormulaFlux flux(formula);
        printf("Formula: %s\n  where x is GeV, and flux is GeV^-1 cm^-2 s^-1\n", flux.GetFormula() );
        SetPointSource(EquatorialDeg(raDeg, decDeg), flux, tPdf, SetNormFromGRL);
    }

    void SetSource(SourceModule *src) {
        (dynamic_cast<I3Analysis*> (psData))->SetSource(*src);
    }
};

class MultiArk : public Ark {
    public:
    Ark* arkPtrArray[10];  // 10 is good for now!
    int n;

    MultiArk() {
        psData = new MultiAnalysisSet();
        for (int i=0; i<10; ++i) { 
        arkPtrArray[i] = NULL;
        }
        n = 0;
    }
    ~MultiArk() { 
        for (int i=0; i<10; ++i) { 
            if (arkPtrArray[i] != NULL ) {break;}
            delete arkPtrArray[i];
        }
    }
  
    void AddArk(Ark& ark) {
        arkPtrArray[n] = &ark;
        (dynamic_cast<MultiAnalysisSet*> (psData))->AddAnalysisSet(ark.psData);
        ++n;
    }

    void SetPointSource(EquatorialDeg srcCoord, const FluxBase& flux) {
        mySrcLocation = srcCoord;
        int i=0;
        while (arkPtrArray[i]) {
            arkPtrArray[i]->SetPointSource(srcCoord, flux);
            ++i;
        }
    }

    void SetPointSource(EquatorialDeg srcCoord, const FluxBase& flux, TimePdf * tPdf, bool SetNormFromGRL=false) {
        mySrcLocation = srcCoord;
        int i=0;
        while (arkPtrArray[i]) {
            arkPtrArray[i]->SetPointSource(srcCoord, flux, tPdf, SetNormFromGRL);
            ++i;
        }
    }

    void SetPointSource(EquatorialDeg srcCoord, const FluxBase& flux, vector<TimePdf*> tPdfVect, bool SetNormFromGRL=false) {
        mySrcLocation = srcCoord;
        int i=0;
        while (arkPtrArray[i]) {
            arkPtrArray[i]->SetPointSource(srcCoord, flux, tPdfVect, SetNormFromGRL);
            ++i;
        }
    }

    void SetPointSource(EquatorialDeg srcCoord, const FluxBase& flux, TimePdf* tPdf, double emin, double emax, bool SetNormFromGRL=false) {
        mySrcLocation = srcCoord;
        int i=0;
        while (arkPtrArray[i]) {
            arkPtrArray[i]->SetPointSource(srcCoord, flux, tPdf, emin, emax, SetNormFromGRL);
            ++i;
        }
    }

    void SetPointSource(EquatorialDeg srcCoord, const FluxBase& flux, vector<TimePdf*> tPdfVect, int iArk, bool SetNormFromGRL=false) {
        mySrcLocation = srcCoord;
        arkPtrArray[iArk]->SetPointSource(srcCoord, flux, tPdfVect, SetNormFromGRL);
    }

    void SetSource(SourceModule *src) {
        int i=0;
        while (arkPtrArray[i]) {
            arkPtrArray[i]->SetSource(src);
            ++i;
        }
    }  
};

class I3Ark : public Ark {
    public:
    EventLoaderExt evLoader;
    vector<I3Event> baseEvents;
    SimpleEnergyProb* eProb;
    DecBkgProb decBkgProb;
    LocalCoordBkgProb lcBkgProb;
    vector<double> startMissRuns;
    vector<double> stopMissRuns;
    vector< vector<TH1D*> > h_recoLogEproxy_;
    vector< vector<TH1D*> > h_PSF_;
    vector< vector<TH1D*> > h_recoAngErr_;
    TH2D* h_Aeff_;

    double sourceZenWidthDeg;
    // This is the +/- zenith range in degrees to select MC events from,
    // and then rotate to desired source location.
    // (Smaller range is more accurate, but trade-off is lower statistics
    // for signal simulation.  You have to pick something appropriate for the 
    // statistics you have.)

    I3Ark() : Ark() ,eProb(NULL) ,  sourceZenWidthDeg(0)   { }

    ~I3Ark() {  }


    void SetPointSource(EquatorialDeg srcCoord, const FluxBase& flux, TH2D* effectiveArea) {
        mySrcLocation = srcCoord;
        cout << "mySrcLocation set to:  "<< mySrcLocation.GetRa()<<" r.a., ";
        cout << mySrcLocation.GetDec() << " dec.\n";
    
        I3PointGenerator i3point(flux, mySrcLocation, livetime, effectiveArea);
        (dynamic_cast<I3Analysis*> (psData))->SetSource(i3point);
    }

    virtual void SetPointSource(EquatorialDeg srcCoord, const FluxBase& flux, TimePdf * tPdf, bool SetNormFromGRL=false) {
        mySrcLocation = srcCoord;
        cout << "mySrcLocation set to:  "<< mySrcLocation.GetRa()<<" r.a., ";
        cout << mySrcLocation.GetDec() << " dec.\n";

        TimePdf * tPdf1 = tPdf->Clone();
        if(SetNormFromGRL) tPdf1->SetNormFromGRL(startMissRuns, stopMissRuns);
        tPdf1->CheckTimeBounds(tmin,tmax);
        Printf("TIMES: min = %f, max = %f", tmin, tmax);
        I3PointGenerator i3point(flux, mySrcLocation, livetime, tPdf1, h_Aeff_);
        i3point.SetRecoLogEproxyDistribution(h_recoLogEproxy_);
        i3point.SetRecoAngErrDistribution(h_recoAngErr_);
        i3point.SetPSFDistribution(h_PSF_);
        i3point.SetEffectiveAreaDistribution(h_Aeff_);
        //Printf("i3point, effective area X bins = %d", h_Aeff_->GetNbinsX());
        (dynamic_cast<I3Analysis*> (psData))->SetSource(i3point);
    }

    void SetPointSource(EquatorialDeg srcCoord, const FluxBase& flux, TimePdf * tPdf, TH2D* effectiveArea, bool SetNormFromGRL=false) {
        mySrcLocation = srcCoord;
        cout << "mySrcLocation set to:  "<< mySrcLocation.GetRa()<<" r.a., ";
        cout << mySrcLocation.GetDec() << " dec.\n";
        
        TimePdf * tPdf1 = tPdf->Clone();
        if(SetNormFromGRL) tPdf1->SetNormFromGRL(startMissRuns, stopMissRuns);
        tPdf1->CheckTimeBounds(tmin,tmax);
        Printf("TIMES: min = %f, max = %f", tmin, tmax);
        I3PointGenerator i3point(flux, mySrcLocation, livetime, tPdf1, effectiveArea);
        (dynamic_cast<I3Analysis*> (psData))->SetSource(i3point);
    }

    void SetPointSource(EquatorialDeg srcCoord, const FluxBase& flux, vector<TimePdf*> tPdfVect, TH2D* effectiveArea, bool SetNormFromGRL=false) {
        mySrcLocation = srcCoord;
        cout << "mySrcLocation set to:  "<< mySrcLocation.GetRa()<<" r.a., ";
        cout << mySrcLocation.GetDec() << " dec.\n";

        vector<TimePdf*> tPdfVect1;
        Printf("Loading time pdf with %d source flares", int(tPdfVect.size()));
        for(unsigned int iFl=0; iFl<tPdfVect.size(); iFl++){
          tPdfVect1.push_back(tPdfVect[iFl]->Clone());
          if(SetNormFromGRL) tPdfVect1[iFl]->SetNormFromGRL(startMissRuns, stopMissRuns);
          tPdfVect1[iFl]->CheckTimeBounds(tmin,tmax);
        }
        I3PointGenerator i3point(flux, mySrcLocation, livetime, tPdfVect1, effectiveArea);
        i3point.SetRecoLogEproxyDistribution(h_recoLogEproxy_);
        i3point.SetRecoAngErrDistribution(h_recoAngErr_);
        i3point.SetPSFDistribution(h_PSF_);
        i3point.SetEffectiveAreaDistribution(h_Aeff_);
        (dynamic_cast<I3Analysis*> (psData))->SetSource(i3point);
        (dynamic_cast<I3Analysis*> (psData))->SetNFlares(tPdfVect.size());
    }

    void SetPointSource(EquatorialDeg srcCoord, const FluxBase& flux, TimePdf* tPdf, double emin, double emax, TH2D* effectiveArea, bool SetNormFromGRL=false) {
        mySrcLocation = srcCoord;
        cout << "mySrcLocation set to:  "<< mySrcLocation.GetRa()<<" r.a., ";
        cout << mySrcLocation.GetDec() << " dec.\n";

        TimePdf* tPdf1;
        tPdf1 = tPdf->Clone();
        if(SetNormFromGRL) tPdf1->SetNormFromGRL(startMissRuns, stopMissRuns);
        tPdf1->CheckTimeBounds(tmin,tmax);
        
        I3PointGenerator i3point(flux, mySrcLocation, livetime, tPdf1, effectiveArea);
        i3point.SetEnergyRange(emin, emax);
        (dynamic_cast<I3Analysis*> (psData))->SetSource(i3point);
    }

    void SetPSFDistribution(vector< vector<TH1D*> > histo){
        h_PSF_ = histo;
    }

    void SetRecoAngErrDistribution(vector< vector<TH1D*> > histo){
        h_recoAngErr_ = histo;
    }

    void SetRecoEnergyDistribution(vector< vector<TH1D*> > histo){
        h_recoLogEproxy_ = histo;
    }

    void SetEffectiveAreaDistribution(TH2D* histo){
        h_Aeff_ = histo;
    }

    void SaveDetectorResponseHisto(char* outfile){
        TFile *fileOutput = new TFile(outfile, "recreate");

        cout << "Writing Histograms  to: " << outfile << endl;

        for(int i=0; i<3; ++i){
            for(unsigned int j=0; j<h_recoLogEproxy_[i].size(); ++j){ h_recoLogEproxy_[i][j]->Write(); } }

        for(int i=0; i<3; ++i){
            for(unsigned int j=0; j<h_recoLogEproxy_[i].size(); ++j){ h_PSF_[i][j]->Write(); } }

        for(int i=0; i<3; ++i){
            for(unsigned int j=0; j<h_recoLogEproxy_[i].size(); ++j){ h_recoAngErr_[i][j]->Write(); } }
 
        h_Aeff_->Write();

        fileOutput->Close();
    }

    vector< vector<TH1D*> > GetPSFDistribution(){        return h_PSF_;           }
    vector< vector<TH1D*> > GetRecoAngErrDistribution(){ return h_recoAngErr_;    }
    vector< vector<TH1D*> > GetRecoEnergyDistribution(){ return h_recoLogEproxy_; }
    TH2D* GetEffectiveAreaDistribution(){               return h_Aeff_;          }
};

#endif // MACRO_ARK_H_
