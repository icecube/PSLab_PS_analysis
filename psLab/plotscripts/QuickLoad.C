
{
  gROOT->ProcessLine(".x loadlibs.C");

  // OLD
  gROOT->ProcessLine(".L /home/cfinley/data/ic40_6months/TreeLoader_IC40_CutA5_6months_final.C");
  TTree *numuIC40old = LoadTree_febIC40_CutA5_nugen1882_small();


  // NEW: numu, nutau
  gROOT->ProcessLine(".L /home/cfinley/data/ic40_6months/CutA5_v2/TreeLoader_IC40_CutA5_v2.C");
  TTree *numuIC40new = LoadTree_CutA5_v2_Nugen1882_01000_small();

  TTree *nutauIC40 = LoadTree_CutA5_v2_nugen_nutau_2183_05000_small();


  // NEW: numu + nutau combined
  gROOT->ProcessLine(".L /home/cfinley/data/ic40_6months/CutA5_v2/TreeLoader_IC40_CutA5_v2_Combine_numu_nutau.C");
  TTree *nuCombinedsmall = LoadTree_CutA5_v2_numu_1k_plus_nutau_1k_small_alternate();
  TTree *nuCombined = LoadTree_CutA5_v2_numu_1k_plus_nutau_5k_small();



  // Default pointer for *tree in effective area script, etc.
  TTree *tree = nuCombined;
}
