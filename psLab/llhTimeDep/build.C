
{

  // THIS PROJECT

  projectSrcDir = "llhTimeDep/public";

  // LIST OF LIBRARIES

  buildList.clear();

  buildList.push_back("TimePdfCollection.h");
  buildList.push_back("BlockLevel.h");
  buildList.push_back("FracUpTime.C");
  buildList.push_back("LocalCoordBkgProb.C");
  buildList.push_back("NewLlhBlockTime.C");
  buildList.push_back("NewLlhGausTime.C");
  buildList.push_back("NewLlhPeriodicTime.C");
  buildList.push_back("GetSrcDeltaEterm.C");
  buildList.push_back("MultiBlockAnalysisFn.C");
  buildList.push_back("MultiGaussAnalysisFn.C");
  buildList.push_back("NewLlhBoxTime.C");
  buildList.push_back("MultiPeriodicAnalysisFn.C");
  buildList.push_back("NewLlhBoxTimeStack.C");
  buildList.push_back("MultiBoxAnalysisFn.C");
  buildList.push_back("MultiBoxAnalysisFnStack.C"); 

  LOADSUCCESS = BasicBuild(projectSrcDir, buildList);
}
