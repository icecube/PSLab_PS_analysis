
{
  bool RELOAD;

  if (!RELOAD) {
    RELOAD = true;
    gROOT->Macro("$LAB_MAIN_DIR/llhTimeDep/loadlibs.C");
    if (!LOADSUCCESS) { return 0.; } // stop, and signal that script failed.

  figpaper1();

  //  TString macroPath = gROOT->GetMacroPath();
  //  macroPath += gSystem->ExpandPathName("$LAB_MAIN_DIR/macro_llh:");
  //  gROOT->SetMacroPath(macroPath);


    if (1) { initialize_ran1(-55); } // seed has to be a *NEGATIVE* integer
    else   { initialize_ran1(); }  // or, no seed = initialize to ran clock time
    int ranSeed = get_ran1_seed(); // if you want to know what the seed was

    gROOT->ProcessLine(".L Ark.C");
    bool OPT_USEREALDATA = false;
    
    
    
    I3Ark ark59;
    //gROOT->ProcessLine(".x load_ark_ic59_noIT.C(ark59, OPT_USEREALDATA)");
    //gROOT->ProcessLine(".x load_ark_ic59_IT2.C(ark59, OPT_USEREALDATA)");
    //gROOT->ProcessLine(".x load_ark_ic59.C(ark59, OPT_USEREALDATA)");
    gROOT->ProcessLine(".x load_ark_ic40.C(ark59, OPT_USEREALDATA)");
    
    /*
    NewLlhEnergy newllh59;
    newllh59.SetUseEnergy(true);
    newllh59.SetOptimizeAngleDeg(10.);
    newllh59.SetOptimizeTolerance(0.01);
    newllh59.SetMonitorLevel(0);
    newllh59.SetAnalysisSet(ark59.psData); */
    
    NewLlhGausTime newllh59;
    newllh59.SetUseEnergy(true);
    //newllh59.SetOptimizeAngleDeg(10.);
    newllh59.SetOptimizeTolerance(0);
    newllh59.SetMonitorLevel(0);
    newllh59.SetEMaxRatioWarnOnlyOnce(1);
    newllh59.close_ = 10.;
    newllh59.JimsTerm_ = true;
    newllh59.SpectralPenalty_ = false;
    newllh59.ndof_ = 3.;
    newllh59.SetLivetime(375.);
    newllh59.SetLocalCoordBkgProb(ark59.lcBkgProb);
    
    ZenithEnergyProb * zn = ark59.eProb;
    SimpleEnergyProb ee = zn->GetSimpleEnergyProb( zn->GetZenDegBand(90+16.148) );
    int nBand = ee.GetHistProbBkg().GetEntries();
    
  }
  
  //  EquatorialDeg testSearch(153.375, 11.375);
  EquatorialDeg testSearch(343.491, 16.148);

  //EquatorialDeg testSearch(224.361, -35.653);

  ark59.SetPointSource(testSearch, PowerLawFlux(1.,-2));
  newllh59.SetAnalysis(ark59.psData,testSearch);
  cout << "SetPointSource" << endl;
  newllh59.SetSearchCoord( testSearch );

  bool optEnergy = true;

  int nTrials = 10000;
  int nSrcEvents = 0;

  ark59.psData->GenerateDataSet_with_nSrcEvents(nSrcEvents);
  newllh59.PrepareAnalysis();
  eventVector = newllh59.GetEventVector(); 

//I3Analysis *DataSet = ark59.psData;//dynamic_cast<I3Analysis*> ark59.;

TH2D * ga = new TH2D("ga","ga",40,2.5,8.5,40,1.,4.);

double nsrc_best, gamma, eweight, w, w1, d;
int ntimes;
vector<double> dvect;
const Coord *cc;
EquatorialDeg *cceq;
EnergyProb * eProb1;

I3Event ev;
vector<I3Event> eventVector;// = baseEvents;
vector<I3Event> srcVect;
double value = 0.;
double valueb = 0.;
I3EventParameters evp;
evp.recoZenithDeg = 90+16.;
double bestg[40], x[40], y[40];
double bestgamma, bkgvalue;

for (int i=0;i<eventVector.size();i++) {
  if (eventVector[i].GetParams().recoZenithDeg > -34. && eventVector[i].GetParams().recoZenithDeg < -36.) {
    ev = eventVector[i];
    ev.SetEnergyProb( eventVector[i].GetEnergyProbFn() );
    break;
  }
}

for (int i=1;i<41;i++) {
  y[i] = ga->GetYaxis()->GetBinCenter(i);
  x[i] = ga->GetXaxis()->GetBinCenter(i);
}

for (int i=1;i<41;i++){ 
//for (int j=1;j<41;j++) {
  valueb = 0.;
  //for (int i=1;i<41;i++){ 
  for (int j=1;j<41;j++) {  
    evp.energyValue = ga->GetXaxis()->GetBinCenter(i);
    ev.SetParams(evp);
    gamma = ga->GetYaxis()->GetBinCenter(j);
    value = ark59.eProb->GetEnergyProbGamma(ev,gamma);
    bkgvalue = ark59.eProb->GetEnergyProbBkg(ev);
    if (!bkgvalue) bkgvalue = 1./nBand;
    if (bkgvalue) ga->SetBinContent(i,j, log10(value/bkgvalue) );
    //ga->SetBinContent(i,j,log10(value) );
    if (valueb < value) {valueb = value; bestgamma = gamma;}
    //cout << j << " " << evp.energyValue << " " << i << " " << gamma << " " << value << " " << valueb << " " << bkgvalue << " " << bestgamma << endl;
  }
  bestg[i-1] = bestgamma;
  //cout << endl << "   " << x[i-1] << " " << y[i-1] << " " << bestgamma << endl ;
}
    

TCanvas * c2 = new TCanvas();
c2->SetGrid(1,1);
ga->GetXaxis()->SetTitle("log_{10} (Reconstructed Energy (GeV))");
ga->GetYaxis()->SetTitle("Spectral index #Gamma");
ga->GetZaxis()->SetTitle("log_{10} Energy Weight");
ga->Draw("colz");

//TGraph * g = new TGraph(40,bestg,y);
TGraph * g = new TGraph(40,x,bestg);
g->SetLineWidth(3);
g->SetLineStyle(2);
g->Draw("LP");

TH1D * hh2 = ga->ProjectionX("hh2",14,14,"");
TH1D * hh3 = ga->ProjectionX("hh3",21,21,"");

  hh2->GetYaxis()->SetTitle("log_{10} Energy Weight");

  hh2->SetLineWidth(2);
  hh3->SetLineWidth(2);
  
  hh3->SetLineStyle(2);
  
  hh2->SetLineColor(9);
  hh3->SetLineColor(4);

TCanvas * c3 = new TCanvas();
hh2->Draw();
hh3->Draw("same");

leg = new TLegend(0.4,0.6,0.8,0.8);
leg->AddEntry(hh3,"E^{-3} Signal","l")-
leg->AddEntry(hh2,"E^{-2} Signal","l");
leg->SetBorderSize(0);
leg->Draw();

//lweights->Draw("same");

//  hTestStatistic.Draw();

//}
  
}

