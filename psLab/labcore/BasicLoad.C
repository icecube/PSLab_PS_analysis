#ifndef LABCORE_BASICLOAD_C_
#define LABCORE_BASICLOAD_C_

void BasicLoad(string this_project, vector<string> project_dependencies) {

  extern vector<string> PROJECTS_VECTOR;  // track projects already loaded
  extern bool LOADSUCCESS;  // set to false as soon as any project load fails

  if (!LOADSUCCESS) {
    cout << "Previous load error prevents current loading of " << 
      this_project << endl;
    return;
  }

  // 1. Check if project has already been loaded

  for (int i=0; i<PROJECTS_VECTOR.size(); ++i) {
    if (PROJECTS_VECTOR[i] == this_project) {
      return;  // don't continue, this project has already been loaded
    }
  }

  
  // 2. Load project dependencies (they perform own "already loaded" checks)

  for (int i=0; i<project_dependencies.size(); ++i) {
    if (LOADSUCCESS) {
      TString macro = "$LAB_MAIN_DIR/"+project_dependencies[i]+"/loadlibs.C";
      gROOT->Macro(macro); // execute macro
    }
  }


  // 3. Build this project

  if (LOADSUCCESS) {
    TString macro = ".x $LAB_MAIN_DIR/"+this_project+"/build.C";
    gROOT->ProcessLine(macro);
    if (LOADSUCCESS) {
      PROJECTS_VECTOR.push_back(this_project);
    }
  }

}

#endif // LABCORE_BASICLOAD_C_

