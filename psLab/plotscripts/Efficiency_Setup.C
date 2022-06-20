
{
  // Example:

  gROOT->ProcessLine(".L $LAB_DATA_DIR/ic40_L2_sim/TreeLoader_IC40_L2_sim.C");
  TTree *treeRaw = LoadTree_IC40_L2_nugen_numu_10k_small();
  int n = treeRaw->GetEntries();
  cout << n << " events loaded in treeRaw.\n";

  gROOT->ProcessLine(".L $LAB_DATA_DIR/ic40ps/TreeLoader_IC40_CutA5_Full_final.C");
  TTree *treeCut = LoadTree_IC40_CutA5_nugen_numu_10k_small();
  int n = treeCut->GetEntries();
  cout << n << " events loaded in treeCut.\n";

  TCut weight = "pow(mcPrimary_Energy_GeV,-1)*mcOneWeight/1e5";
  //  TCut cutRaw = "(FilterMinBias_08==1)*FMBScale";
  TCut cutRaw = "(EHEFilter_08==1 || ICMuonFilter_08==1)";
  TCut cutCut = "1";

}
