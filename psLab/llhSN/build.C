
{

  // THIS PROJECT

  TString projectSrcDir = "llhSN/public";

  // LIST OF LIBRARIES

  vector<string> buildList;

  buildList.push_back("SNEvent.h");
  buildList.push_back("SimpleProb.C");
  buildList.push_back("SNAnalysis.C");
  buildList.push_back("NewLlhSN.C");

  // For using vectors in CINT macros:
  // (Put any new classes which you need from this project in here
   buildList.push_back("SNVectorLoader.C"); 
  // why does this work when I place vector<SNEvent> in llh but not llhSN?

  LOADSUCCESS = BasicBuild(projectSrcDir, buildList);
}
