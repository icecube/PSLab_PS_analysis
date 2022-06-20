
{

  // THIS PROJECT

  projectSrcDir = "llh/public";

  // LIST OF LIBRARIES

  buildList.clear();
 
  buildList.push_back("CoordinateTransform.C");
  buildList.push_back("CoordClasses.C");
  buildList.push_back("CoordEquatorialDeg.h");
  buildList.push_back("classes.h");
  //buildList.push_back("Time.h");
  buildList.push_back("EventTimeModule.C");
  buildList.push_back("EventTimeModuleDiscrete.C");
  //buildList.push_back("TimePdf.h"); 
  //buildList.push_back("I3Event.h");
  buildList.push_back("EventLoader.C");
  buildList.push_back("I3SignalGenerator.C");
  buildList.push_back("DeclinationDensityMap.C");
  buildList.push_back("EnergyProb.C"); // abstract base class only
  buildList.push_back("SimpleEnergyProb.C");
  buildList.push_back("ZenithEnergyProb.C");
  buildList.push_back("BkgSpaceProb.h");
  buildList.push_back("DecBkgProb.C");
  buildList.push_back("I3Analysis.C");

  buildList.push_back("MinuitWrapperClass.h");
  buildList.push_back("LlhFunctionsBase.C");
  buildList.push_back("LlhFunctions.C");
  buildList.push_back("LlhEnergy.C");

  buildList.push_back("MinuitAnalysisFn.C");
  buildList.push_back("NewLlhEnergy.C");
  buildList.push_back("MultiAnalysisSet.C");
  buildList.push_back("MultiAnalysisFn.C");

  buildList.push_back("llh_discovery_potential.C");

  //buildList.push_back("Kent.C");
  buildList.push_back("BetaProf.C");
  buildList.push_back("EventLoaderExt.C");

  // For using vectors in CINT macros:
  buildList.push_back("VectorLoader.C");

  LOADSUCCESS = BasicBuild(projectSrcDir, buildList);
}
