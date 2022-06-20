//String LOADTREE  = "/home/lucarelf/data/ps_tracks/version-003-p02/";
TString LOADTREE  = "data/";

TTree* LoadTree_nugen_numu(vector <TString> nufiles) {

    TChain *tree  =new TChain("tree");

    for (unsigned int j=0;j<nufiles.size();j++){
        tree->Add(LOADTREE+nufiles[j]);
    }
    
    tree->GetEntries();


    tree->SetAlias("logMuexEn",  "logE"                    ); 
    tree->SetAlias("zenDeg"   ,  "zen*TMath::RadToDeg()"   );
    tree->SetAlias("aziDeg"   ,  "azi*TMath::RadToDeg()"   );
    tree->SetAlias("zenRad"   ,  "zen"                     );
    tree->SetAlias("aziRad"   ,  "azi"                     );
    tree->SetAlias("raRad"    ,  "ra"                      );
    tree->SetAlias("decRad"   ,  "dec"                     );
    tree->SetAlias("sigmaDeg" ,  "angErr*TMath::RadToDeg()");
    tree->SetAlias("RunID"    ,  "run"                     );
    tree->SetAlias("EventID"  ,  "event"                   );
    tree->SetAlias("timeMJD"  ,  "time"                   );
     
    //needed 
 
    return tree;
}


TTree* LoadTree_GoodRuns_Full(TString period, vector<double> &startVect, vector<double> &stopVect) {
  
  TChain *tree     = new TChain("tree");
  TChain *tree_GRL = new TChain("tree");

  if(period=="IC40") {
    TString files = LOADTREE+"IC40_exp.root";
    tree->Add(files);
    files = LOADTREE+"GRL/IC40_exp.root";
    tree_GRL->Add(files);
  }
  else if(period=="IC59") {
    TString files = LOADTREE+"IC59_exp.root";
    tree->Add(files);

    files = LOADTREE+"GRL/IC59_exp.root";
    tree_GRL->Add(files);
  }
  else if(period=="IC79") {
    TString files = LOADTREE+"IC79_exp.root";
    tree->Add(files);

    files = LOADTREE+"GRL/IC79_exp.root";
    tree_GRL->Add(files);
  }
  else if(period=="IC86_I") {
    TString files = LOADTREE+"IC86_2011_exp.root";
    tree->Add(files);

    files = LOADTREE+"GRL/IC86_2011_exp.root";
    tree_GRL->Add(files);
  }
  else if(period=="IC86_II_IV") {
    TString files2012 = LOADTREE+"IC86_2012_exp.root";
    TString files2013 = LOADTREE+"IC86_2013_exp.root";
    TString files2014 = LOADTREE+"IC86_2014_exp.root";
    tree->Add(files2012);
    tree->Add(files2013);
    tree->Add(files2014);

    files2012 = LOADTREE+"GRL/IC86_2012_exp.root";
    files2013 = LOADTREE+"GRL/IC86_2013_exp.root";
    files2014 = LOADTREE+"GRL/IC86_2014_exp.root";
    tree_GRL->Add(files2012);
    tree_GRL->Add(files2013);
    tree_GRL->Add(files2014);
  }
  else if(period=="IC86_V_VII") {
    TString files2015 = LOADTREE+"IC86_2015_exp.root";
    TString files2016 = LOADTREE+"IC86_2016_exp.root";
    TString files2017 = LOADTREE+"IC86_2017_exp.root";
    tree->Add(files2015);
    tree->Add(files2016);
    tree->Add(files2017);

    files2015 = LOADTREE+"GRL/IC86_2015_exp.root";
    files2016 = LOADTREE+"GRL/IC86_2016_exp.root";
    files2017 = LOADTREE+"GRL/IC86_2017_exp.root";
    tree_GRL->Add(files2015);
    tree_GRL->Add(files2016);
    tree_GRL->Add(files2017);
  }
  else if(period=="IC86_II") {
    TString files2012 = LOADTREE+"IC86_2012_exp.root";
    tree->Add(files2012);

    files2012 = LOADTREE+"GRL/IC86_2012_exp.root";
    tree_GRL->Add(files2012);
  }
  else if(period=="IC86_III") {
    TString files2013 = LOADTREE+"IC86_2013_exp.root";
    tree->Add(files2013);

    files2013 = LOADTREE+"GRL/IC86_2013_exp.root";
    tree_GRL->Add(files2013);
  }
  else if(period=="IC86_IV") {
    TString files2014 = LOADTREE+"IC86_2014_exp.root";
    tree->Add(files2014);

    files2014 = LOADTREE+"GRL/IC86_2014_exp.root";
    tree_GRL->Add(files2014);
  }
  else if(period=="IC86_V") {
    TString files2015 = LOADTREE+"IC86_2015_exp.root";
    tree->Add(files2015);

    files2015 = LOADTREE+"GRL/IC86_2015_exp.root";
    tree_GRL->Add(files2015);
  }
  else if(period=="IC86_VI") {
    TString files2016 = LOADTREE+"IC86_2016_exp.root";
    tree->Add(files2016);

    files2016 = LOADTREE+"GRL/IC86_2016_exp.root";
    tree_GRL->Add(files2016);
  }
  else if(period=="IC86_VII") {
    TString files2017 = LOADTREE+"IC86_2017_exp.root";
    tree->Add(files2017);

    files2017 = LOADTREE+"GRL/IC86_2017_exp.root";
    tree_GRL->Add(files2017);
  }
  else if(period=="IC86_II_VII") {
    TString files2012 = LOADTREE+"IC86_2012_exp.root";
    TString files2013 = LOADTREE+"IC86_2013_exp.root";
    TString files2014 = LOADTREE+"IC86_2014_exp.root";
    TString files2015 = LOADTREE+"IC86_2015_exp.root";
    TString files2016 = LOADTREE+"IC86_2016_exp.root";
    TString files2017 = LOADTREE+"IC86_2017_exp.root";
    tree->Add(files2012);
    tree->Add(files2013);
    tree->Add(files2014);
    tree->Add(files2015);
    tree->Add(files2016);
    tree->Add(files2017);

    files2012 = LOADTREE+"GRL/IC86_2012_exp.root";
    files2013 = LOADTREE+"GRL/IC86_2013_exp.root";
    files2014 = LOADTREE+"GRL/IC86_2014_exp.root";
    files2015 = LOADTREE+"GRL/IC86_2015_exp.root";
    files2016 = LOADTREE+"GRL/IC86_2016_exp.root";
    files2017 = LOADTREE+"GRL/IC86_2017_exp.root";
    tree_GRL->Add(files2012);
    tree_GRL->Add(files2013);
    tree_GRL->Add(files2014);
    tree_GRL->Add(files2015);
    tree_GRL->Add(files2016);
    tree_GRL->Add(files2017);
  }
  else if(period=="IC86-VIII-GFU") {
    LOADTREE  = "/data/ana/analyses/gfu/version-002-p04/";
    TString files = LOADTREE+"IC86_2018_data.root";
    tree->Add(files);

    files = LOADTREE+"GRL/IC86_2018_data.root";
    tree_GRL->Add(files);
  }
  else if(period=="IC86-VII-GFU") {
    LOADTREE  = "/data/ana/analyses/gfu/version-002-p04/";
    TString files = LOADTREE+"IC86_2017_data.root";
    tree->Add(files);

    files = LOADTREE+"GRL/IC86_2017_data.root";
    tree_GRL->Add(files);
  }

  tree->GetEntries();  // REQUIRED TO "INITIALIZE" CHAIN
  tree_GRL->GetEntries();  // REQUIRED TO "INITIALIZE" CHAIN

  // SET META-INFORMATION:
  // These aliases store meta-information for the tree
  double livetime = 0.;
  double tmin     = 999999.;
  double tmax     = 0.;
  double lt, start, stop;

  tree_GRL->SetBranchAddress("livetime",&lt);
  tree_GRL->SetBranchAddress("start",&start);
  tree_GRL->SetBranchAddress("stop",&stop);
  
  double prevStop = 0;
  startVect.clear();
  stopVect.clear();
  vector<double> t_start_sort, t_stop_sort;

  for(int i=0; i<tree_GRL->GetEntriesFast();i++) {
    tree_GRL->GetEntry(i);
    livetime += lt;
    t_start_sort.push_back(start);
    t_stop_sort.push_back(stop);
    
    //if(start>=tmax) prevStop = stop;
    if(start<tmin) tmin=start; 
    if(stop>tmax) tmax=stop; 
  }

  sort(t_start_sort.begin(), t_start_sort.end());
  sort(t_stop_sort.begin(), t_stop_sort.end()); 
  for(int i=0; i<t_start_sort.size(); ++i){
    if(i>0){
      start = t_start_sort[i];
      if(start-prevStop>=1./10.){
        startVect.push_back(prevStop);
        stopVect.push_back(start);
      }
    }
    prevStop = t_stop_sort[i];
  }
 
  Printf("From GRL tree, livetime is %.4f s",livetime*86400.);
  TString livetimeTotalStr = Form(" 1. * %.4f",livetime*86400.); // livetime in sec
  TString tminStr = Form(" 1. * %.4f",tmin); //in MJD
  TString tmaxStr = Form(" 1. * %.4f",tmax); //in MJD

  tree->SetAlias("livetimeTotal", livetimeTotalStr );
  tree->SetAlias("tmin"         , tminStr          );
  tree->SetAlias("tmax"         , tmaxStr          );

  //--------------------//
  tree->SetAlias("logMuexEn",  "logE"                    ); 
  tree->SetAlias("zenDeg"   ,  "zen*TMath::RadToDeg()"   );
  tree->SetAlias("aziDeg"   ,  "azi*TMath::RadToDeg()"   );
  tree->SetAlias("zenRad"   ,  "zen"                     );
  tree->SetAlias("aziRad"   ,  "azi"                     );
  tree->SetAlias("raRad"    ,  "ra"                      );
  tree->SetAlias("decRad"   ,  "dec"                     );
  tree->SetAlias("sigmaDeg" ,  "angErr*TMath::RadToDeg()");
  tree->SetAlias("RunID"    ,  "run"                     );
  tree->SetAlias("EventID"  ,  "event"                   );
  tree->SetAlias("timeMJD"  ,  "time"                    );
  
  return tree;
}
