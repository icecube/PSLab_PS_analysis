
{
  gROOT->Macro( "$LAB_PRELOAD" );

  this_project = "coord_interface";

  project_dependencies.clear();
  // none

  // For now, only thing this project does is load coordinate-service
  // icetray library and provide public headers for interface

  BasicLoad(this_project, project_dependencies);
}
