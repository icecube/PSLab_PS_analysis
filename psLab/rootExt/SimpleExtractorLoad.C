
{
  
  // note that for macs, you may prefer to run ".x loadmac.C" first,
  // which seems to fix some uncommon bugs


  // SET DIRECTORIES

  TString currentDir = gSystem->pwd();

  TString rootExtDir = currentDir+"/rootExt";  // modify as necessary

  cout << "cd to rootExt directory: " << rootExtDir << endl;
  gSystem->cd(rootExtDir);

  gROOT->ProcessLine(".L MakeNewName.h+");
  gROOT->ProcessLine(".L TStringify.h+");
  gROOT->ProcessLine(".L FunctionsRoot.C+");
  gROOT->ProcessLine(".L FunctionsGeneral.C+");
  gROOT->ProcessLine(".L SimpleExtractor.C+");

  cout << "cd to original directory: " << currentDir << endl;
  gSystem->cd(currentDir);

}
