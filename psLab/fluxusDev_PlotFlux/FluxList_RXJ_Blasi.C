
{

  FormulaFlux fluxRXJ_Blasi("1.33662e-11*2.25855e2*pow(x,0.21539-2)*exp(-sqrt(x/2659.37))");

  /*
  //RXJ Blasi

  Double_t blasi_loge[9] ={10.013,10.848,11.616,12.182,12.505,12.869,13.677,14.121,14.391};
  Double_t blasi_flux[9] ={0.679,0.777,0.875,0.851,0.777,0.606,-0.372,-1.522,-2.5};
  Int_t Nbin_blasi = 9;

  for (Int_t i = 0; i < 9; i++){
  //TeV
    blasi_loge[i]  = blasi_loge[i]-12.;
    blasi_flux[i] = blasi_flux[i]-12.;
  }
  TGraph *blasi_plot = new TGraph(Nbin_blasi,blasi_loge,blasi_flux);
  TMultiGraph *mg1 = new TMultiGraph();
  mg1->Add(blasi_plot, "P");

 
  TF1 *blasi_func = new TF1("blasi_func","log10([0]*pow(pow(10.,x),[1])*exp(-sqrt(pow(10.,x)/[2])))",-2.,2.5);
  //*exp(-pow(10.,x)/[2]))",10.,14.5);
  blasi_func->SetParameter(0,1e-12);
  //blasi_func->SetParLimits(0,10.,20.);
  blasi_func->SetParameter(1,0.3);
  //blasi_func->SetParLimits(1,0.,.5);
  blasi_func->SetParameter(2,1.3);
  //blasi_func->SetParLimits(2,0.5,3.); 

  blasi_plot->Fit(blasi_func,"MR+");

  Double_t par0 = blasi_func->GetParameter(0);
  Double_t par1 = blasi_func->GetParameter(1);
  Double_t par2 = blasi_func->GetParameter(2);
  Double_t newflux[9];

  blasi_func->SetParameter(0,par0);
  blasi_func->SetParameter(1,par1);
  blasi_func->SetParameter(2,par2);

  for (Int_t i = 0; i < 9; i++){
    cout << blasi_loge[i] << endl;
    newflux[i] = blasi_func->Eval(blasi_loge[i]);
    cout << newflux[i] << endl;
  }

  TGraph *blasi_plot1 = new TGraph(Nbin_blasi,blasi_loge,newflux);

  cout << " par0 " << par0 << " par1 " << par1 << " par2 " << par2 << endl;

  TCanvas *can_X = new TCanvas("can_X", "");
  can_X->cd();
  TAxis *aXk = (blasi_func->GetXaxis());
  //aXk->SetLimits(10.,14.5);
  aXk->SetTitle("Log(E)[eV]");
  TAxis *aYk = (blasi_func->GetYaxis());
  aYk->SetTitle("Log(E2 dN/dE [eV cm^{-2} s^{-1}");
  blasi_plot->SetMarkerStyle(20);
  blasi_plot1->SetMarkerStyle(20);
  blasi_plot1->SetLineColor(2);
  blasi_func->SetLineColor(3);
  blasi_plot->SetMinimum(-14.);
  blasi_plot->SetMaximum(-10.);

  blasi_plot->Draw("APL");
  blasi_plot1->Draw("PLsame");  
  */
}
