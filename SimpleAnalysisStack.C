#include "macro_llh/ArkTimeExt.C"
#include "llhTimeDep/public/NewLlhBoxTime.h"
#include "llhTimeDep/public/NewLlhBoxTimeStack.h"
#include "llhTimeDep/public/MultiBoxAnalysisFn.h"
#include "llhTimeDep/public/MultiBoxAnalysisFnStack.h"
#include "TTree.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"

class SimpleAnalysis_multiSet {

public:

  TTree tFitCoord;

  //TH1F *hfluxScale_;

  int  ns_;
  int  nset_;

  //Double_t pvalue_;
  Double_t TS_;
  Double_t gammaBest_;
  Double_t nsBest_;
  Double_t fluxScale_;
  Double_t fluenceScale_;
  Double_t nsSet_[5];

  bool branchesNotSet_;
  bool rndOnlyTimes_;
  bool isUnblind_;

  char* fileToMostSigTreeName_;
  char* mostSigTreeName_;

  SimpleAnalysis_multiSet() : ns_(0), branchesNotSet_(true), isUnblind_(false) { }
  ~SimpleAnalysis_multiSet() { 
    //delete hfluxScale_;
  }
  
  void SetDoUnblind() {
    isUnblind_ = true;
  }

  void SetRndOnlyTimes(bool b){
    rndOnlyTimes_ = b;
  }

  void Initialize(int ns, int nset, double fluencescale, double fluxscale){
    ns_     = ns; 
    nset_   = nset;
    fluenceScale_ = fluencescale;
    fluxScale_    = fluxscale;
    //hfluxScale_ = new TH1F("hfluxScale_","Flux scale",1,0,1);
   
  }
  void Execute(MultiArk& ark, MultiBoxAnalysisFnStack& llh);
  bool Write(char *filename, char* fileoption);

  void SetFileName(char*filename);
  void SetTreeName(char*treename);
};

void SimpleAnalysis_multiSet::Execute(MultiArk& ark, MultiBoxAnalysisFnStack& llh)
{
  //generate scramble
  if(!isUnblind_){
    MultiAnalysisSet* multiPsData = (dynamic_cast<MultiAnalysisSet*> (ark.psData));
    multiPsData->GenerateDataSetStacking_with_nSrcEvents(ns_);
  }

  if(branchesNotSet_) {
    tFitCoord.SetName("tFitCoord");
    tFitCoord.Branch("gammaBest",&gammaBest_, "gammaBest/D");
    tFitCoord.Branch("nsBest",&nsBest_, "nsBest/D");
    tFitCoord.Branch("fluxScale",&fluxScale_, "fluxScale/D");
    tFitCoord.Branch("fluenceScale",&fluenceScale_, "fluenceScale/D");
    tFitCoord.Branch("TS",&TS_,"TS/D");
    branchesNotSet_=false;
  }

  
  //maximizing the llh
  llh.MaximizeLlh();
  double resultLLH = 2*llh.GetTestStatistic();
  Printf("Fill the tree with:");
  nsBest_ = llh.Get_nSrcBest();
  gammaBest_ = llh.Get_gammaBest();
  Printf("\tns Best %f",nsBest_);
  Printf("\tgammaBest %f",gammaBest_);
  TS_        = 2*llh.GetTestStatistic();
 
  Printf("\tTS %f",TS_);

  tFitCoord.Fill(); 
}


bool SimpleAnalysis_multiSet::Write(char*filename, char* fileoption) {
  TFile *fileOutput = new TFile(filename, fileoption);
  if (fileOutput->IsZombie()) {
    cout << "Try using 'recreate' option?\n";
    return false; // no file saved                                                                                                                                                                          
  }

  cout << "Writing Histograms  to: " << filename << endl;

  //hfluxScale_->Write();
  tFitCoord.Write();
  
  fileOutput->Close();
  return true;  
                                                                                                                                                                 
}

void SimpleAnalysis_multiSet::SetFileName(char*filename) { fileToMostSigTreeName_ = filename; }
 
void SimpleAnalysis_multiSet::SetTreeName(char*treename) { mostSigTreeName_ = treename; } 
