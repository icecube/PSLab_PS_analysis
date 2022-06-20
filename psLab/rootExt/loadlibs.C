
{
  gROOT->Macro( "$LAB_PRELOAD" );

  this_project = "rootExt";

  project_dependencies.clear();
  // none

  BasicLoad(this_project, project_dependencies);
}
