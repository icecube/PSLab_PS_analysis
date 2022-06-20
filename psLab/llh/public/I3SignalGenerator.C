#include "llh/public/I3SignalGenerator.h"
#include "llh/public/Kent.C"
#include "llh/public/CoordinateTransform.h"
#include "llh/public/CoordEquatorialDeg.h"

#include "rootExt/public/randomfunctions.h"
#include "fluxus/public/FluxFunction.h"
#include "TFile.h"

bool sortpair(const pair<double, double> &a, const pair<double,double> &b)
{
  return (a.first < b.first);
}

//
// I3MultiSignalGenerator
//



// Need a deep copy, since new copies must be made for signalPtr's
I3SignalGenerator* I3MultiSignalGenerator::Clone() const {
  I3MultiSignalGenerator* newMultiPtr = new I3MultiSignalGenerator();
  newMultiPtr->livetime_ = livetime_;
  for (unsigned int i = 0; i< signalPtrVect_.size(); ++i) {
    newMultiPtr->AddSignal(*(signalPtrVect_[i]), enhanceFactorVect_[i]);
  }
  return newMultiPtr;
}


void I3MultiSignalGenerator::SetLivetime(double livetime) { // THIS IS FRAUGHT
  livetime_ = livetime;
  for (unsigned int i = 0; i< signalPtrVect_.size(); ++i) { 
    signalPtrVect_[i]->SetLivetime(livetime_);
  }
}

void I3MultiSignalGenerator::SetTimePdf(TimePdf * tPdf) {
  timePdf_ = tPdf;
  for (unsigned int i = 0; i< signalPtrVect_.size(); ++i) {
    signalPtrVect_[i]->SetTimePdf(tPdf);
  }
  
  double factor;
  for (unsigned int i = 0; i< signalPtrVect_.size(); ++i) {  
    factor = signalPtrVect_[i]->GetTimePdf()->GetNorm();
    SetEnhanceFactor(i, factor);
  }
}

void I3MultiSignalGenerator::AddSignal(const I3SignalGenerator& signal, 
					  double enhanceFactor) 
{
  // This makes a private copy of the added signal; must delete in destructor
  I3SignalGenerator* newSignalPtr = signal.Clone();
  signalPtrVect_.push_back(newSignalPtr);
  enhanceFactorVect_.push_back(enhanceFactor);
  nSources_ = signalPtrVect_.size();
}

double I3MultiSignalGenerator::GetMeanSrcNevBase() const { // HMMM
  double nev = 0.;
  for (unsigned int i = 0; i< signalPtrVect_.size(); ++i) {
    nev += signalPtrVect_[i]->GetMeanSrcNev();
  }
  return nev;
}

double I3MultiSignalGenerator::GetMeanSrcNev() const {
  double nev = 0.;
  for (unsigned int i = 0; i< signalPtrVect_.size(); ++i) {
    if(signalPtrVect_[i]->GetTimePdf()) nev += enhanceFactorVect_[i] * signalPtrVect_[i]->GetMeanSrcNev() * signalPtrVect_[i]->GetTimePdf()->GetNorm();
    else nev += enhanceFactorVect_[i] * signalPtrVect_[i]->GetMeanSrcNev();
  }
  return nev; 
}

double I3MultiSignalGenerator::GetMeanSrcNev(int isig) {
  if(isig>=int(signalPtrVect_.size())) log_error("ERROR: accessing a signal injector outside signalPtrVect_ size\n"); 
  return enhanceFactorVect_[isig] * signalPtrVect_[isig]->GetMeanSrcNev();
}


double I3MultiSignalGenerator::
GetMeanSrcNevForFluxModelBase(const FluxBase& fluxModel) const { // HMMM
  double nev = 0.;
  for (unsigned int i = 0; i<signalPtrVect_.size(); ++i) {
    nev += signalPtrVect_[i]->GetMeanSrcNevForFluxModel(fluxModel);
  }
  return nev;
}

double I3MultiSignalGenerator::
GetMeanSrcNevForFluxModel(const FluxBase& fluxModel) const {
  double nev = 0.;
  for (unsigned int i = 0; i<signalPtrVect_.size(); ++i) {
    nev += enhanceFactorVect_[i] * 
      signalPtrVect_[i]->GetMeanSrcNevForFluxModel(fluxModel);
  }
  return nev;
}

double I3MultiSignalGenerator::GetFluenceNormalization(int isig){
  if(isig>=int(signalPtrVect_.size())) log_error("ERROR: accessing a signal injector outside signalPtrVect_ size\n"); 
  return enhanceFactorVect_[isig] * signalPtrVect_[isig]->GetFluenceNormalization();
}

double I3MultiSignalGenerator::GetFluenceNormalization(){
  double nev = 0.;
  for (unsigned int i = 0; i< signalPtrVect_.size(); ++i) {
    nev += enhanceFactorVect_[i] * signalPtrVect_[i]->GetFluenceNormalization();
  }
  return nev;
}



I3Event I3MultiSignalGenerator::GenerateEvent() {
  double meanNevSum = GetMeanSrcNev();
  double ranNum = random_uniform(0., meanNevSum);

  // start adding up the nev (as in GetMeanSrcNev), and stop when 
  // random number is less than sum

  double sum = 0.;
  unsigned int i;
  for (i = 0; i< signalPtrVect_.size(); ++i) {
    sum += enhanceFactorVect_[i] * signalPtrVect_[i]->GetMeanSrcNev() * signalPtrVect_[i]->GetTimePdf()->GetNorm();
    if (ranNum <= sum) {
      break;
    }
  }
  // 'i' is now the signal which we want to use:
  return signalPtrVect_[i]->GenerateEvent();
}




//
// I3PointGenerator
//



I3PointGenerator::I3PointGenerator() : 
  livetime_(0),
  isTimeDep_(false),
  injectEnergyRange_(false),
  angErrSrc_(0.0),
  timeAzBins_(1),
  sourceTimePdf_(NULL),
  nsFluenceRatio_(0) 
{
  sourceTimePdfVect_.push_back(NULL);
}


I3PointGenerator::I3PointGenerator(
  const FluxBase& fluxModel,
  const EquatorialDeg& sourceCoord, 
  double livetime,
  TH2D* Aeff) :
  livetime_(livetime), isTimeDep_(true), injectEnergyRange_(false), angErrSrc_(0.0), extendedSrc_(false), timeAzBins_(1), nsFluenceRatio_(0)
{
  sourceTimePdf_ = NULL;
  sourceTimePdfVect_.push_back(NULL);
  SetSourceCoord(sourceCoord);
  SetFluxModel(fluxModel);
  SetEffectiveAreaDistribution(Aeff);
}

I3PointGenerator::I3PointGenerator(
  const FluxBase& fluxModel,
  const EquatorialDeg& sourceCoord,
  double livetime,
  TimePdf * tPdf) :
  livetime_(livetime), isTimeDep_(true), injectEnergyRange_(false), angErrSrc_(0.0), extendedSrc_(false), timeAzBins_(1), nsFluenceRatio_(0)
{
  SetTimePdf(tPdf);
  sourceTimePdfVect_.push_back(NULL);
  SetSourceCoord(sourceCoord);
  SetFluxModel(fluxModel);
}

I3PointGenerator::I3PointGenerator(
  const FluxBase& fluxModel,
  const EquatorialDeg& sourceCoord,
  double livetime,
  TimePdf* tPdf,
  TH2D* Aeff) : 
  livetime_(livetime), isTimeDep_(true), injectEnergyRange_(false), angErrSrc_(0.0), extendedSrc_(false), timeAzBins_(1), nsFluenceRatio_(0)
{
  SetTimePdf(tPdf);
  sourceTimePdfVect_.push_back(NULL);
  SetSourceCoord(sourceCoord);
  SetFluxModel(fluxModel);
  SetEffectiveAreaDistribution(Aeff);
}


I3PointGenerator::I3PointGenerator(
  const FluxBase& fluxModel,
  const EquatorialDeg& sourceCoord, 
  double livetime,
  TimePdf *tPdf,
  double emin,
  double emax,
  TH2D* Aeff) :
  livetime_(livetime), isTimeDep_(true), injectEnergyRange_(false), Emin_(emin), Emax_(emax), angErrSrc_(0.0), extendedSrc_(false), timeAzBins_(1), nsFluenceRatio_(0)
{
  SetTimePdf(tPdf);
  sourceTimePdfVect_.push_back(NULL);
  SetSourceCoord(sourceCoord);
  SetFluxModel(fluxModel);
  SetEffectiveAreaDistribution(Aeff);
}

I3PointGenerator::I3PointGenerator(
  const FluxBase& fluxModel,
  const EquatorialDeg& sourceCoord,
  double livetime,
  vector<TimePdf*> tPdfVect,
  TH2D* Aeff) :
  livetime_(livetime), isTimeDep_(true), injectEnergyRange_(false), angErrSrc_(0.0), extendedSrc_(false), timeAzBins_(1), nsFluenceRatio_(0)
{
  sourceTimePdf_ = NULL;
  SetTimePdfVect(tPdfVect);
  SetSourceCoord(sourceCoord);
  SetFluxModel(fluxModel);
  SetEffectiveAreaDistribution(Aeff);
}

double I3PointGenerator::GetMeanSrcNevForFluxModel(const FluxBase& fluxModel) const {
  double decmin=sourceCoord_.GetDec()-0.5;
  if(decmin<-90) decmin=-90;
  double decmax=sourceCoord_.GetDec()+0.5;
  if(decmax>90) decmax=90;

  TH2D Aeff;
  dynamic_cast<I3PointGenerator*>(Clone())->GetEffectiveAreaDistribution()->Copy( Aeff );
  Int_t binymin = Aeff.GetYaxis()->FindBin( TMath::Sin( TMath::DegToRad()*(decmin+0.1) ) );
  Int_t binymax = Aeff.GetYaxis()->FindBin( TMath::Sin( TMath::DegToRad()*(decmax-0.1) ) );
  int Nxbins = Aeff.GetNbinsX();

  double newSinDecBin = Aeff.GetYaxis()->GetBinLowEdge(binymax+1)-Aeff.GetYaxis()->GetBinLowEdge(binymin);
  double effAreaBin, oldSinDecBin, energy, BinWidth;
  double nsFluenceRatio = 0;
  for(int binx=1; binx<Nxbins+1; binx++){
    for(int biny=binymin; biny<binymax+1; biny++){
      effAreaBin = Aeff.GetBinContent(binx, biny);
      oldSinDecBin = Aeff.GetYaxis()->GetBinLowEdge(biny+1)-Aeff.GetYaxis()->GetBinLowEdge(biny);
      energy = pow(10, Aeff.GetXaxis()->GetBinCenter(binx) );
      BinWidth = pow(10, Aeff.GetXaxis()->GetBinLowEdge(binx+1)) - pow(10, Aeff.GetXaxis()->GetBinLowEdge(binx));
      nsFluenceRatio += effAreaBin * oldSinDecBin/newSinDecBin * fluxModel.GetFlux(energy) * BinWidth;
    }
  }
  return nsFluenceRatio;
}

double I3PointGenerator::GetFluenceNormalization(){
    TH1D* AeffTimesFluence = GenerateEnuPDF();
    int Nxbins = AeffTimesFluence->GetNbinsX();
    double nsFluenceRatio=0;

    for(int binx=1; binx<Nxbins+1; binx++) nsFluenceRatio += AeffTimesFluence->GetBinContent(binx);

    return nsFluenceRatio; 
}

TH1D* I3PointGenerator::GenerateEnuPDF(){
    double decmin=sourceCoord_.GetDec()-0.5;
    if(decmin<-90) decmin=-90;
    double decmax=sourceCoord_.GetDec()+0.5;
    if(decmax>90) decmax=90;

    TH2D Aeff;
    GetEffectiveAreaDistribution()->Copy( Aeff );
    Int_t binymin = Aeff.GetYaxis()->FindBin( TMath::Sin( TMath::DegToRad()*(decmin+0.1) ) );
    Int_t binymax = Aeff.GetYaxis()->FindBin( TMath::Sin( TMath::DegToRad()*(decmax-0.1) ) );
    int Nxbins = Aeff.GetNbinsX();

    double newSinDecBin = Aeff.GetYaxis()->GetBinLowEdge(binymax+1)-Aeff.GetYaxis()->GetBinLowEdge(binymin);
    double effAreaBin, oldSinDecBin, energy, BinWidth;
    for(int binx=1; binx<Nxbins+1; binx++){
      for(int biny=binymin; biny<binymax+1; biny++){
        effAreaBin = Aeff.GetBinContent(binx, biny);
        oldSinDecBin = Aeff.GetYaxis()->GetBinLowEdge(biny+1)-Aeff.GetYaxis()->GetBinLowEdge(biny);
        energy = pow(10, Aeff.GetXaxis()->GetBinCenter(binx) );
        BinWidth = pow(10, Aeff.GetXaxis()->GetBinLowEdge(binx+1)) - pow(10, Aeff.GetXaxis()->GetBinLowEdge(binx));
        Aeff.SetBinContent(binx, biny, effAreaBin*oldSinDecBin/newSinDecBin*BinWidth*fluxModel_->GetFlux(energy));
      }
    }

    return Aeff.ProjectionX("logEnuPDF", binymin, binymax);
}

I3Event I3PointGenerator::GenerateEvent() {
  //Randomly sample relevant parameters from IceCube response distributions
  I3Event e;

  //true neutrino direction = source location
  I3MCParameters mcPars;
  mcPars.PS_FlatSpectrumRate = 1;
  mcPars.srcWeight = 1;

  I3EventParameters evtPars;
  evtPars.runID   = 0;
  evtPars.eventID = 0;

  //determine hemisphere
  int hemi; //hemi=0 -> south (dec<-10), hemi=1 -> horizon (-10<dec<10), hemi=2 -> north (dec>10)
  if(sourceCoord_.GetDec()<-10.) hemi=0;
  else if(sourceCoord_.GetDec()<10.) hemi=1;
  else hemi=2;

  //Printf("*** Generating signal event ***");
  //Printf("dec = %f, hemi = %d", sourceCoord_.GetDec(), hemi);

  //Generate time
  Time t;
  do{
      t = sourceTimePdf_->GenerateEventTime();
  } while(!CheckTimeInGRL(t.GetMJD()));

  e.SetTime( t );
  //Printf("time = %f MJD", t.GetMJD());

  //neutrino direction: ...
  double nuTrueRA = 0, nuTrueDec = 0;
  if(extendedSrc_){
    //...random from Kent distribution, for extended sources with angular size angErrSrc_...
    TH2D *hKent = (TH2D*)Kent(sourceCoord_.GetDec()+90., sourceCoord_.GetRa(), angErrSrc_);
    hKent->GetRandom2(nuTrueRA, nuTrueDec); // Random Phi and Theta in radians
    nuTrueRA  = nuTrueRA*TMath::RadToDeg();
    nuTrueDec = nuTrueDec*TMath::RadToDeg()-90;
  }
  else{
    //...or fixed at source location for point-like sources
    nuTrueRA  = sourceCoord_.GetRa();
    nuTrueDec = sourceCoord_.GetDec();
  }
  //Printf("neutrino direction: RA = %.3f deg, dec = %.3f deg", nuTrueRA, nuTrueDec);   

  double thetaNuRad = (90.-nuTrueDec)*TMath::DegToRad(); //theta = zenith, in spherical coordinates
  double phiNuRad   = (nuTrueRA-180.)*TMath::DegToRad(); //phi = azimuth, in spherical coordinates, [-180, 180)
  Spherical nuDir(thetaNuRad, phiNuRad);

  //Generate neutrino energy
  TH1D* nuEPDF = GenerateEnuPDF();
  double trueLogEnuGeV = nuEPDF->GetRandom();
  const int energybin = int(trueLogEnuGeV/0.5)-4;
  mcPars.mcEnergy = pow(10, trueLogEnuGeV);
  //Printf("Nu random energy: %f GeV", pow(10, trueLogEnuGeV));

  //Generate reco muon direction 
  double muNuAngDist = h_PSF_[hemi][energybin]->GetRandom();
  //Printf("mu-nu angular separation = %.2f", muNuAngDist);
  Spherical muRecoDir = nuDir.RandomDirectionAtRadius(muNuAngDist);
  Double_t muRecoRA  = muRecoDir.GetPhiDeg() + 180.;
  Double_t muRecoDec = 90. - muRecoDir.GetThetaDeg();
  e.SetCoord(EquatorialDeg(muRecoRA, muRecoDec)); 
  //Printf("Reco muon direction: RA = %.3f deg, dec = %.3f deg", muRecoRA, muRecoDec); 
  double muRecoZen, muRecoAzi;
  EqToLocal(muRecoRA, muRecoDec, t.GetMJD(), muRecoZen, muRecoAzi);
  //Printf("Reco muon direction: Az = %.3f deg, Zen = %.3f deg", muRecoAzi, muRecoZen);

  evtPars.parafitSigmaDeg = h_recoAngErr_[hemi][energybin]->GetRandom();
  evtPars.recoZenithDeg   = muRecoZen;
  evtPars.recoAzimuthDeg  = muRecoAzi;

  //Generate muon energy
  double recoLogEmuon = h_recoLogEproxy_[hemi][energybin]->GetRandom();
  //Printf("Reco log10GeV muon energy = %.2f", recoLogEmuon);
  evtPars.energyValue = recoLogEmuon;

  //Assign evtPars and mcPars to event
  e.SetParams(evtPars);
  e.SetMCParams(mcPars);

  return e; 
}

bool I3PointGenerator::CheckTimeInGRL(double mjd){
  TimePdf *tPdf=0;
  if(sourceTimePdf_) tPdf = sourceTimePdf_->Clone();
  else if(sourceTimePdfVect_[0]) tPdf = sourceTimePdfVect_[0]->Clone();
  else log_error("ERROR: no src time PDF assigned\n");

  vector<double> startMissRunVect = tPdf->GetStartMissRunVect();
  vector<double> stopMissRunVect = tPdf->GetStopMissRunVect();

  for(unsigned int i=0; i<startMissRunVect.size(); ++i){
    if(mjd > startMissRunVect[i] && mjd < stopMissRunVect[i]) return false;
  }
  return true;
}

void I3PointGenerator::SaveDetectorResponseHisto(char* outfile){
  TFile *fileOutput = new TFile(outfile, "recreate");

  cout << "Writing Histograms  to: " << outfile << endl;

  for(unsigned int i=0; i<h_recoLogEproxy_.size(); ++i){
      for(unsigned int j=0; j<h_recoLogEproxy_[i].size(); ++j){ h_recoLogEproxy_[i][j]->Write(); } }

  for(unsigned int i=0; i<h_PSF_.size(); ++i){
      for(unsigned int j=0; j<h_PSF_[i].size(); ++j){ h_PSF_[i][j]->Write(); } }

  for(unsigned int i=0; i<h_recoAngErr_.size(); ++i){
      for(unsigned int j=0; j<h_recoAngErr_[i].size(); ++j){ h_recoAngErr_[i][j]->Write(); } }
  
  h_Aeff_->Write();

  fileOutput->Close();
}
