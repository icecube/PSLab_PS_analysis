// This macro loads a specific dataset (loaded into ark59), and makes a number
// of scrambled datasets to plot a distribution of the individual event signal over
// background ratios used inside the likelihood method against the space angle
// of the event to the tested location. The idea is to provide a simple visual
// plot of the unbinned weights vs a binned analysis


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
    
  }
  
  //  EquatorialDeg testSearch(153.375, 11.375);
    EquatorialDeg testSearch(343.491, 16.148);
  //  EquatorialDeg testSearch(224.361, -35.653);

  ark59.SetPointSource(testSearch, PowerLawFlux(1.,-2));
  newllh59.SetAnalysis(ark59.psData,testSearch);
  cout << "SetPointSource" << endl;
  newllh59.SetSearchCoord( testSearch );

  bool optEnergy = true;

  int nTrials = 10000;
  int nSrcEvents = 0;


I3Analysis *DataSet = ark59.psData;//dynamic_cast<I3Analysis*> ark59.;
TH1D * weights = new TH1D("weights","weights",1800,0.,180.);
TH1D * lweights = new TH1D("lweights","lweights",400,0.,40.);
TH2D * tve = new TH2D("tve","tve",200,0.,20.,200,-20.,10.);

TH1D * weights1 = new TH1D("weights1","weights1",1800,0.,180.);
TH1D * lweights1 = new TH1D("lweights1","lweights1",400,0.,40.);
TH2D * tve1 = new TH2D("tve1","tve1",200,0.,20.,200,-20.,10.);

TH1D * evts = new TH1D("evts","evts",1800,0.,180.);

double nsrc_best, gamma, eweight, w, w1, d;
int ntimes;
vector<double> dvect;
const Coord *cc;
EquatorialDeg *cceq;
EnergyProb * eProb1;

I3Event ev;
vector<I3Event> eventVector;// = baseEvents;
vector<I3Event> srcVect;

  CountMonitor countMon(1., nTrials);

  for (int i=0;i<nTrials;i++) {

    ark59.psData->GenerateDataSet_with_nSrcEvents(nSrcEvents);
    
    newllh59.PrepareAnalysis();
    eventVector = newllh59.GetEventVector();  
      
    for (unsigned int j=0;j<eventVector.size();j++) {
    
      d = eventVector[j].GetCoord().DistanceTo(testSearch);
      eweight = eventVector[j].GetEnergyProbFn()->GetEnergyMaxRatio(eventVector[j]);

      w = ( eventVector[j].ProbFrom(testSearch) * eweight / eventVector[j].GetBkgSpaceProbFn()->GetBkgProbDensity(eventVector[j].GetCoord()) );
      
      w1 = ( eventVector[j].ProbFrom(testSearch) / eventVector[j].GetBkgSpaceProbFn()->GetBkgProbDensity(eventVector[j].GetCoord()) );
            
      weights->Fill( d, log10(w) );
      weights1->Fill( d, log10(w1) );
      tve->Fill(d, log10(w) );
      tve1->Fill(d, log10(w1) );
      
    }

    countMon.UpdateCount();

}

 cout << "end of loop?" << endl;

evts->SetLineColor(4);
evts->SetLineWidth(3);
weights->SetLineWidth(3);
weights1->SetLineWidth(3);
weights1->SetLineStyle(2);
lweights->SetLineWidth(3);

weights->GetXaxis()->SetTitle("#Delta#Psi (#circ)");

evts->Scale(1./nTrials);
weights->Scale(1./nTrials);
weights->Divide(evts);

weights1->Scale(1./nTrials);
weights1->Divide(evts);

for (int i=1;i<400;i++) {
  lweights->SetBinContent(i, log10(weights->GetBinContent(i)) );
  lweights1->SetBinContent(i, log10(weights1->GetBinContent(i)) );
}

TCanvas * c1 = new TCanvas();
//c1->SetLogy();
c1->SetGrid(1,1);

weights->Draw();
weights1->Draw("same");
evts->Draw("same");

  leg = new TLegend(0.1,0.7,0.48,0.9);
  leg->AddEntry(weights,"<Weight per event> per bin, use Energy","l");
  leg->AddEntry(weights1,"<Weight per event> per bin, no Energy","l");
  leg->AddEntry(evts,"<# Events> per bin","l");
  leg->Draw();

TCanvas * c2 = new TCanvas();
c2->SetGrid(1,1);
tve->GetYaxis()->SetTitle("Time Integrated Event Weight");
tve->GetXaxis()->SetTitle("#Delta#Psi (#circ)");
tve->Draw("colz");
lweights->Draw("same");

//  hTestStatistic.Draw();

//}
  
}

