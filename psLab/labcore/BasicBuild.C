#ifndef LABCORE_BASICBUILD_C_
#define LABCORE_BASICBUILD_C_


// BasicBuild sets up include paths, build paths, creates directories
// if needed, compiles code (or loads, if no source code has changed 
// since last compile), and checks that the code was loaded.  If not,
// it returns false, and further compiling/loading should stop.

// The vector of compile code should be simply something like:
//   vector<string> buildList;
//   buildList.push_back("llh_general.C");
//   buildList.push_back("llh_random.C");
//   buildList.push_back("I3Event.h");  // if no corresponding C file exists


bool BasicBuild(TString projectSrcDir, vector<string> buildList) {

  cout << "Building/Loading:  " << projectSrcDir << endl;


  // Check that $LAB_MAIN_DIR is defined

  TString LAB_MAIN_DIR = getenv("LAB_MAIN_DIR");
  TString fullSrcDir = LAB_MAIN_DIR + "/" + projectSrcDir;

  // Include is now done globally by preload.C
  //  gSystem->AddIncludePath(" -I" + fullSrcDir + " ");



  // Check that $LAB_LIB_DIR is defined

  TString LAB_LIB_DIR = getenv("LAB_LIB_DIR");
  if (LAB_LIB_DIR == "") {
    cout << "ERROR: $LAB_LIB_DIR not defined in environment.\n";
    return false; // indicate script ended in failure
  }

  TString thisLibDir = LAB_LIB_DIR + "/" + projectSrcDir;

  // Create thisLibDir, if needed
  if (gSystem->AccessPathName(thisLibDir)) { // true means doesn't exist!!!
    cout << "** Making Directory " << thisLibDir << endl;
    // -p option: make parent directories as necessary
    gSystem->Exec(TString("mkdir -p ")+thisLibDir);
  }

  gSystem->SetBuildDir(thisLibDir);
  // Set the location where ACLiC will create libraries and use as
  // a scratch area.  Note that the libraries are actually stored in
  // sub-directories of 'build_dir' including the full pathname of the
  // to-be-compiled script.
  // e.g. If the code location is at /full/path/name/myCode.C
  // the library will be located at 'build_dir+/full/path/name/myCode_C.so'


  // Create project subdirectory of $LAB_LIB_DIR, if needed 

  if (gSystem->AccessPathName(thisLibDir)) { // true means doesn't exist!!!
    cout << "** Making Directory " << thisLibDir << endl;
    gSystem->mkdir(thisLibDir);
  }



  // COMPILE AND LOAD


  bool success = true;
  for (int i=0; i<int(buildList.size()); ++i) {
    TString fullPath = fullSrcDir + "/" + TString(buildList[i]);
    gROOT->ProcessLine(".L " + fullPath + "+");

    // check that libary was added correctly
    // (I wish the return value of gROOT->ProcessLine would flag this!!!)
    TString libraries = gSystem->GetLibraries();

    libraries.ReplaceAll("/./","/");
    libraries.ReplaceAll("//","/");
    libraries.ReplaceAll("_C."+TString(gSystem->GetSoExt()),".C");
    libraries.ReplaceAll("_h."+TString(gSystem->GetSoExt()),".h");
    if ( ! libraries.Contains(fullPath) ) {
      success = false;
      cout << "\nERROR.  " << fullPath << endl;
      cout << "does not appear to have been added to the libaries. Exit.\n";
      break;
    }
  }

  return success;     
}


#endif // LABCORE_BASICBUILD_C_
