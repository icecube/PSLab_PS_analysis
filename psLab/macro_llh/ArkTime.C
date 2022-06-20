#include "llh/public/classes.h"
#include "llh/public/CoordEquatorialDeg.h"
#include "llh/public/SimpleEnergyProb.h"
#include "llh/public/DecBkgProb.h"
#include "fluxus/public/FluxFunction.h"
#include "llh/public/I3Analysis.h"
#include "llh/public/MultiAnalysisSet.h"
#include "llh/public/EventLoader.h"
#include "llhTimeDep/public/LocalCoordBkgProb.h"

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

    virtual void SetPointSource(EquatorialDeg srcCoord, const FluxBase& flux, TimePdf * tPdf, bool SetNormFromGRL=false){};

    virtual void SetPointSource(EquatorialDeg srcCoord, const FluxBase& flux, vector<TimePdf*> tPdfVect, bool SetNormFromGRL=false){};

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
    EventLoader evLoader;
    vector<I3Event> baseEvents;
    SimpleEnergyProb* eProb;
    DecBkgProb decBkgProb;
    LocalCoordBkgProb lcBkgProb;
    vector<double> startMissRuns;
    vector<double> stopMissRuns;

    double sourceZenWidthDeg;
    // This is the +/- zenith range in degrees to select MC events from,
    // and then rotate to desired source location.
    // (Smaller range is more accurate, but trade-off is lower statistics
    // for signal simulation.  You have to pick something appropriate for the 
    // statistics you have.)

    I3Ark() : Ark() ,eProb(NULL) ,  sourceZenWidthDeg(0)   { }

    ~I3Ark() {  }


    void SetPointSource(EquatorialDeg srcCoord, const FluxBase& flux) {
        mySrcLocation = srcCoord;
        cout << "mySrcLocation set to:  "<< mySrcLocation.GetRa()<<" r.a., ";
        cout << mySrcLocation.GetDec() << " dec.\n";
    
        vector<I3Event> sourceEvents;
        evLoader.LoadSourceEvents(sourceEvents, mySrcLocation);
        I3PointGenerator i3point(sourceEvents, flux, mySrcLocation, livetime);
        (dynamic_cast<I3Analysis*> (psData))->SetSource(i3point);
    }

    void SetPointSource(EquatorialDeg srcCoord, const FluxBase& flux, TimePdf * tPdf, bool SetNormFromGRL=false) {
        mySrcLocation = srcCoord;
        cout << "mySrcLocation set to:  "<< mySrcLocation.GetRa()<<" r.a., ";
        cout << mySrcLocation.GetDec() << " dec.\n";
        
        vector<I3Event> sourceEvents;
        evLoader.LoadSourceEvents(sourceEvents, mySrcLocation);
        TimePdf * tPdf1 = tPdf->Clone();
        if(SetNormFromGRL) tPdf1->SetNormFromGRL(startMissRuns, stopMissRuns);
        tPdf1->CheckTimeBounds(tmin,tmax);
        Printf("TIMES: min = %f, max = %f", tmin, tmax);
        I3PointGenerator i3point(sourceEvents, flux, mySrcLocation, livetime, tPdf1);
        (dynamic_cast<I3Analysis*> (psData))->SetSource(i3point);
    }

    void SetPointSource(EquatorialDeg srcCoord, const FluxBase& flux, vector<TimePdf*> tPdfVect, bool SetNormFromGRL=false) {
        mySrcLocation = srcCoord;
        cout << "mySrcLocation set to:  "<< mySrcLocation.GetRa()<<" r.a., ";
        cout << mySrcLocation.GetDec() << " dec.\n";

        vector<I3Event> sourceEvents;
        evLoader.LoadSourceEvents(sourceEvents, mySrcLocation);
        vector<TimePdf*> tPdfVect1;
        Printf("Loading time pdf with %d source flares", int(tPdfVect.size()));
        for(unsigned int iFl=0; iFl<tPdfVect.size(); iFl++){
          tPdfVect1.push_back(tPdfVect[iFl]->Clone());
          if(SetNormFromGRL) tPdfVect1[iFl]->SetNormFromGRL(startMissRuns, stopMissRuns);
          tPdfVect1[iFl]->CheckTimeBounds(tmin,tmax);
        }
        I3PointGenerator i3point(sourceEvents, flux, mySrcLocation, livetime, tPdfVect1);
        (dynamic_cast<I3Analysis*> (psData))->SetSource(i3point);
        (dynamic_cast<I3Analysis*> (psData))->SetNFlares(tPdfVect.size());
    }

    void SetPointSource(EquatorialDeg srcCoord, const FluxBase& flux, TimePdf* tPdf, double emin, double emax, bool SetNormFromGRL=false) {
        mySrcLocation = srcCoord;
        cout << "mySrcLocation set to:  "<< mySrcLocation.GetRa()<<" r.a., ";
        cout << mySrcLocation.GetDec() << " dec.\n";

        vector<I3Event> sourceEvents;
        evLoader.LoadSourceEvents(sourceEvents, mySrcLocation);
        TimePdf* tPdf1;
        tPdf1 = tPdf->Clone();
        if(SetNormFromGRL) tPdf1->SetNormFromGRL(startMissRuns, stopMissRuns);
        tPdf1->CheckTimeBounds(tmin,tmax);
        
        I3PointGenerator i3point(sourceEvents, flux, mySrcLocation, livetime, tPdf1);
        i3point.SetEnergyRange(emin, emax);
        (dynamic_cast<I3Analysis*> (psData))->SetSource(i3point);
    }
};


