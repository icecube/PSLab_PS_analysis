#include "macro_llh/ArkTimeExt.C"
#include "llhTimeDep/public/MultiGaussAnalysisFn.h"
#include "TTree.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"

class SimpleAnalysis_multiSet {

public:

  TTree tFitCoord;

  int  ns_;
  int  nset_;

  double srcDEC_;
  double srcRA_;
  Double_t injRA_;
  Double_t injDEC_;

  Double_t TS_;
  vector<double> gammaBest_;
  vector<double> nsBest_;
  vector<double> sigmaBest_;
  vector<double> timeBest_;
  vector<double> gammaErr_;
  vector<double> nsErr_;
  vector<double> sigmaErr_;
  vector<double> timeErr_;
  Double_t nsSet_[5];
  Double_t fluenceScale_;
  Double_t fluxScale_;
  Int_t nFlares_;

  bool branchesNotSet_;
  bool rndOnlyTimes_;
  bool isUnblind_;

  char* fileToMostSigTreeName_;
  char* mostSigTreeName_;

  SimpleAnalysis_multiSet() : ns_(0), branchesNotSet_(true), isUnblind_(false) { }
  ~SimpleAnalysis_multiSet() {}
  
  void SetDoUnblind() {
    isUnblind_ = true;
  }

  void SetRndOnlyTimes(bool b){
    rndOnlyTimes_ = b;
  }

  void Initialize(int ns, double srcRA, double srcDEC, int nset, double fluencescale, double fluxscale){
    ns_     = ns;
    srcRA_  = srcRA;
    srcDEC_ = srcDEC;
    injRA_  = srcRA_;
    injDEC_ = srcDEC_;
    nset_   = nset;
    fluenceScale_ = fluencescale;
    fluxScale_    = fluxscale;
  }
  void Execute(MultiArk& ark, MultiGaussAnalysisFn& llh);
  bool Write(char *filename, char* fileoption);

  void SetFileName(char*filename);
  void SetTreeName(char*treename);
};

void SimpleAnalysis_multiSet::Execute(MultiArk& ark, MultiGaussAnalysisFn& llh)
{
  if(!isUnblind_){
    MultiAnalysisSet* multiPsData = (dynamic_cast<MultiAnalysisSet*> (ark.psData));
    //multiPsData->SetRndOnlyTime(rndOnlyTimes_);
    //multiPsData->SetSelEvFileName(fileToMostSigTreeName_);
    //multiPsData->SetSelEvTreeName(mostSigTreeName_);
    multiPsData->GenerateDataSet_with_nSrcEvents(ns_);
  }
  
  std::vector<double> spaceWeightVect;
  std::vector<double> eneWeightVect;
  std::vector<double> timeWeightVect;
  std::vector<double> raVect;
  std::vector<double> decVect;
  std::vector<double> timeVect;
  std::vector<int>    eventID;
  std::vector<int>    runID;

  if(branchesNotSet_) {
    tFitCoord.SetName("tFitCoord");
    tFitCoord.Branch("gammaBest",&gammaBest_);
    tFitCoord.Branch("nsBest",&nsBest_);
    tFitCoord.Branch("timeBest",&timeBest_);
    tFitCoord.Branch("sigmaBest",&sigmaBest_);
    tFitCoord.Branch("gammaErr",&gammaErr_);
    tFitCoord.Branch("nsErr",&nsErr_);
    tFitCoord.Branch("timeErr",&timeErr_);
    tFitCoord.Branch("sigmaErr",&sigmaErr_);
    tFitCoord.Branch("nFlares",&nFlares_, "nFlares/I");
    tFitCoord.Branch("TS",&TS_,"TS/D");
    tFitCoord.Branch("spaceWeightVect",&spaceWeightVect );
    tFitCoord.Branch("eneWeightVect",&eneWeightVect);
    tFitCoord.Branch("timeWeightVect",&timeWeightVect);
    tFitCoord.Branch("raVect",&raVect);
    tFitCoord.Branch("decVect",&decVect);
    tFitCoord.Branch("timeVect",&timeVect);
    tFitCoord.Branch("eventID",&eventID);
    tFitCoord.Branch("runID",&runID);
    tFitCoord.Branch("srcRA",&srcRA_, "srcRA/D");
    tFitCoord.Branch("srcDEC",&srcDEC_, "srcDEC/D");
    tFitCoord.Branch("fluenceScale",&fluenceScale_, "fluenceScale/D");
    tFitCoord.Branch("fluxScale",&fluxScale_, "fluxScale/D");
    branchesNotSet_=false;
  }

  
  //maximizing the llh
  llh.MaximizeLlh();
  int ires=llh.GetMinimizationResult();
  if(ires==0){
    double result = -TMath::Log10(llh.GetEstProb());
    double resultLLH = 2*llh.GetTestStatistic();
    Printf("Fill the tree with:");
    nFlares_=llh.GetNFlareGuess();
    for(int ifl=0; ifl<llh.Get_nSrcBest().size(); ifl++){
      Printf("*** flare %d ***", ifl+1);
      nsBest_.push_back(llh.Get_nSrcBest().at(ifl).first);
      gammaBest_.push_back(llh.Get_gammaBest().at(ifl).first);
      timeBest_.push_back(llh.Get_meanBest().at(ifl).first);
      sigmaBest_.push_back(llh.Get_sigmaBest().at(ifl).first);
      nsErr_.push_back(llh.Get_nSrcBest().at(ifl).second);
      gammaErr_.push_back(llh.Get_gammaBest().at(ifl).second);
      timeErr_.push_back(llh.Get_meanBest().at(ifl).second);
      sigmaErr_.push_back(llh.Get_sigmaBest().at(ifl).second);
      Printf("\tns Best %f +- %f",nsBest_[ifl], nsErr_[ifl]);
      Printf("\tgammaBest %f +- %f",gammaBest_[ifl], gammaErr_[ifl]);
      Printf("\ttimeBest %f +- %f",timeBest_[ifl], timeErr_[ifl]);
      Printf("\tsigmaBest %f +- %f",sigmaBest_[ifl], sigmaErr_[ifl]);
    }
    TS_        = 2*llh.GetTestStatistic();
 
    Printf("\tTS %f",TS_);

    spaceWeightVect.clear();
    eneWeightVect.clear();
    timeWeightVect.clear();
    raVect.clear();
    decVect.clear();
    timeVect.clear(); 
    eventID.clear();
    runID.clear();

    spaceWeightVect = llh.GetSpatialWeights();
    eneWeightVect = llh.GetEnergyWeights();
    timeWeightVect = llh.GetTimeWeights();
    raVect = llh.GetraVect();
    decVect = llh.GetdecVect();
    timeVect = llh.GettimeVect();  
    eventID = llh.GetEventID();
    runID = llh.GetRunID();

    tFitCoord.Fill();

    nsBest_.clear();
    gammaBest_.clear();
    timeBest_.clear();
    sigmaBest_.clear();
    nsErr_.clear();
    gammaErr_.clear();
    timeErr_.clear();
    sigmaErr_.clear();
  }
  else{
    Printf("fit did not converge!");
    double result = 0;
    double resultLLH = 0;
    Printf("Fill the tree with:");
    nFlares_=llh.GetNFlareGuess();
    Printf("*** flare 1 ***");
    nsBest_.push_back(0);
    gammaBest_.push_back(0);
    timeBest_.push_back(0);
    sigmaBest_.push_back(0);
    nsErr_.push_back(0);
    gammaErr_.push_back(0);
    timeErr_.push_back(0);
    sigmaErr_.push_back(0);
    Printf("\tns Best %f +- %f",nsBest_[0], nsErr_[0]);
    Printf("\tgammaBest %f +- %f",gammaBest_[0], gammaErr_[0]);
    Printf("\ttimeBest %f +- %f",timeBest_[0], timeErr_[0]);
    Printf("\tsigmaBest %f +- %f",sigmaBest_[0], sigmaErr_[0]);
    //Printf("\tns/flux = %.2e", nsFluxRatio);
    
    TS_        = 0;

    //Printf("pvalue %f",pvalue_);
    Printf("\tTS %f",TS_);

    //filling the weights
    spaceWeightVect.clear();
    eneWeightVect.clear();
    timeWeightVect.clear();
    raVect.clear();
    decVect.clear();
    timeVect.clear();
    eventID.clear();
    runID.clear();

    spaceWeightVect = llh.GetSpatialWeights();
    eneWeightVect = llh.GetEnergyWeights();
    timeWeightVect = llh.GetTimeWeights();
    raVect = llh.GetraVect();
    decVect = llh.GetdecVect();
    timeVect = llh.GettimeVect();
    eventID = llh.GetEventID();
    runID = llh.GetRunID();

    tFitCoord.Fill();

    nsBest_.clear();
    gammaBest_.clear();
    timeBest_.clear();
    sigmaBest_.clear();
    nsErr_.clear();
    gammaErr_.clear();
    timeErr_.clear();
    sigmaErr_.clear();
  }
}


bool SimpleAnalysis_multiSet::Write(char*filename, char* fileoption) {
  TFile *fileOutput = new TFile(filename, fileoption);
  if (fileOutput->IsZombie()) {
    cout << "Try using 'recreate' option?\n";
    return false; // no file saved                                                                                                                                                                          
  }

  cout << "Writing Histograms  to: " << filename << endl;

  tFitCoord.Write();
  
  fileOutput->Close();
  return true;  
                                                                                                                                                                 
}

void SimpleAnalysis_multiSet::SetFileName(char*filename) { fileToMostSigTreeName_ = filename; }
 
void SimpleAnalysis_multiSet::SetTreeName(char*treename) { mostSigTreeName_ = treename; } 
