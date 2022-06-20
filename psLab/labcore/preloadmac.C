
{
  #ifndef LABCORE_PRELOADMAC_C_
  #define LABCORE_PRELOADMAC_C_

  cout << "preloadmac.C\n";

  gSystem->Load("libTreePlayer.so");

  // This fixes a problem when using TTreeFormula.  On my mac, this library
  // does not get automatically loaded, causing problems later...


  gSystem->SetSoExt("dl");  // for mac only?  otherwise use default "so" ?

  // this seems to fix a problem, on my mac anyway, where changes to the
  // code aren't reflected when the file is recompiled... as if a ghost
  // of the original file keeps getting used for the compiled .so library.

  gROOT->Macro("$LAB_CORE_DIR/preload.C");

  #endif // LABCORE_PRELOADMAC_C_
}
