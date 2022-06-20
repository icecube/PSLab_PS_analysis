
{
  gROOT->Macro( "$LAB_PRELOAD" );

  string this_project = "fluxusDev_PlotFlux";

  vector<string> project_dependencies;
  project_dependencies.push_back("rootExt");
  project_dependencies.push_back("rootDev");
  project_dependencies.push_back("fluxus");
  project_dependencies.push_back("llh");

  BasicLoad(this_project, project_dependencies);
}
