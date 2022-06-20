
// Welcome to MultiPeriodicAnalysisFn
// Here are the functions to maximize and find the proper seed values
// for periodic searches utilizing more than one dataset. It can also
// retrieve the I3Events for each dataset used.


#include "llhTimeDep/public/MultiPeriodicAnalysisFn.h"
#include "rootExt/public/FunctionsRoot.h"
#include "rootExt/public/generalfunctions.h"
#include "iostream"

void llhFncMP(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  f = 0.;
  vector<double> parVect;
  for(int i=0; i<npar; i++) parVect.push_back(par[i]);
  for(int i=0; i<int(ptr->analysisFnVect_.size()); ++i) {
    vector<double> individualParVect = ptr->parTrans_->Translate(i, parVect);
    f += ptr->analysisFnVect_[i]->EvalFCN(individualParVect);

    //Now we correct for the marginalization. 
    //Since it was applied N times on each individual llh
    //we have to remove it (N-1) times. 
    //Since result = -LogLikelihood it needs to be added (+). 
    //All marginalization weights are identical
  

    if (i>0 && ptr->analysisFnVect_[i]->JimsTerm_ && ptr->analysisFnVect_[i]->Get_MargWeight() < 0.) { 
      f += ptr->analysisFnVect_[i]->Get_MargWeight();
    }

  }
  return;
}



MultiPeriodicAnalysisFn::MultiPeriodicAnalysisFn() {
  
  nPar = 4;
  histoForProb_ = false;

  minuit_ = new TMinuit(nPar); 

  Double_t arglist[1];
  arglist[0] = 0.5; //0.5 for likelihood fit
  Int_t ierflg=0;
  minuit_->SetFCN(llhFncMP);
  minuit_->mnexcm("SET ERR",arglist,1,ierflg);
  minuit_->SetPrintLevel(-1); 
}


double MultiPeriodicAnalysisFn::EvalFCN(const vector<double>& parVect) const {
  double f;
  const int npar = parVect.size();
  double gin[npar];
  double par[npar];
  int iflag = 0;
  for(int i=0; i<npar; i++) par[i] = parVect[i];
  int npars=npar;
  llhFncMP(npars, gin, f, par, iflag);

  return f;
}

double MultiPeriodicAnalysisFn::EvaluateLlh(double *parValueArray) {
  vector<double> parVect;
  
  for (int i=0; i<nPar; ++i) {
    parVect.push_back(parValueArray[i]);
  }
  double minusLlh = EvalFCN(parVect);
  return -minusLlh;   // that is, max llh = - (minimizer result)
}


void MultiPeriodicAnalysisFn::StoreLogLambdaBest() 
{
  Double_t amin, edm, errdef;
  Int_t nvpar, nparx, istat;
  minuit_->mnstat(amin, edm, errdef, nvpar, nparx,istat); 
  logLambdaBest_ = -amin;  // that is, max llh = - (minimizer result)

  if ( analysisFnVect_[0]->monitorLevel_ > 1 ) { 
    printf("LogLambdaBest=%12.6lg\n",logLambdaBest_ );
  } 

  // The *worst* logLambdaBest should be zero (i.e. null hypothesis).
  // But minimizer will miss exact zero, leading logLambdaBest_ to be slightly
  // negative. Since this will cause probability calculation to choke, we 
  // fix it here
  
  if (logLambdaBest_ < 0.) { logLambdaBest_ = 0.;}
}


void MultiPeriodicAnalysisFn::SetNullTestStat(TH1D * inputhisto) {
  
  // This reads in a specific TH1D as the null test statistic
  // to use for p-values instead of using a chisquare distribution.
  // It also used to fit an exponential to the upper tail (default set to top 0.01%).
  
  // It takes the raw test statistic distribution (pdf) as input
  // and then makes the cdf to work with later.
  
  histoForProb_ = true;
  
  pvalHisto_ = new TH1D();
  nullTestStat_ = new TH1D();
  
  char sc[] = "Scale";
  //double firstfit=0;
  pvalHisto_ = DescendingCumulate(inputhisto,sc);
  nullTestStat_ = (TH1*)inputhisto->Clone("hnew");
  
  //This fitting bit is commented out since it started breaking everything.
  //Perhaps one day it will work again...
  
  /*
  int bins = inputhisto->GetNbinsX();
  for (int i=1;i<bins;i++) {
    if (pvalHisto->GetBinContent(i) < fitfrac_){
      firstfit = pvalHisto->GetBinCenter(i);
      break;
    }
  }
   
  double max = 2*inputhisto->GetBinLowEdge(bins) - inputhisto->GetBinLowEdge(bins-1);
  fitfn_ = new TF1("fitfn_","exp([0]+[1]*x)",firstfit,100);
   
  pvalHisto->Fit(fitfn_,"QO","",firstfit,max);
  */
  
}

double MultiPeriodicAnalysisFn::GetProbFromHisto(double teststat) const { //, bool useFit) {

  // this looks at a null teststat distribution
  // which we loaded, and if it's in the exponential
  // tail, should go with a fit expo->Eval(teststat) (top 1%)
  
  int bin = pvalHisto_->FindBin(teststat);
  double ptemp = pvalHisto_->GetBinContent(bin);
  
//  if (ptemp < fitfrac_ && useFit) { //small enough that we should do the fit stuff
//    ptemp = fitfn_->Eval(teststat);
//  }
  return ptemp;
  
  // in case we don't feel like using a fit and it is higher than the distribution
  // just return 1/trials (conservative?)
  if( !pvalHisto_->GetBinContent(bin) ) { return 1.0/nullTestStat_->GetEntries(); }
  
  return -1.; //if we get here something is wrong, 
              // and may as well hear about it trying to take the log
}


//    * Oh lordy, if we assume that we want to use something
//    * like the untriggered search over multiple years,
//    * that's a search which uses the event info to seed the
//    * first guess. I'm sure this will also be useful for plotting.
//    *
//    * This also assumes you don't want some sort of extra relative 
//    * weighting that comes from the detector config/cuts, could be
//    * cool to have that too (in the events?).
//    *


vector<I3Event> MultiPeriodicAnalysisFn::GetAllEvents() {

  vector<I3Event> allEvents;
  vector<I3Event> tempVect;

  for (int i=0; i<int(analysisFnVect_.size()); ++i) {
    NewLlhPeriodicTime* i3an = dynamic_cast<NewLlhPeriodicTime*>(analysisFnVect_[i]);
    
    //analysisFnVect_[i]->PrepareAnalysis();
    if (i3an->Get_nEvents()==0) { i3an->PrepareAnalysis(); }
    tempVect.clear();
    tempVect = i3an->GetEventVector();
    
    for (int j=0; j<int(tempVect.size()); j++) {
      allEvents.push_back(tempVect[j]);
    }
  }
  
  return allEvents;
  
}

// This function scans over all the events which could have an S/B ratio > 10,
// and tests consecutive doubles, triples... for compatability with and E^-2 flare.
// I tried testing for Gamma = 2 and like 3.5 before, but it doesn't seem to do much.
// This is configurable enough to do whatever you'd like.
// There is a final scan over Gamma at the end to round out the set parameters.

// Oh, I'm using this parameter close_ to set how far into the vector 1,2,3,4,5,10,15,20,25
// to test for flare compatability. You only need to go to 10, I usually go to 15 (close_=6).

//*/

void MultiPeriodicAnalysisFn::GetFlareGuessPeriodic(double & Guess_nsrc, 
                   double & Guess_gamma, double & Guess_mean, 
                     double & Guess_rms, double & sigmamin) {


  // Note: already have Coord * srcCoord_ set;

  vector<I3Event> evs = GetAllEvents();

  vector<double> tVectortemp;
  vector<double> tVectorclose;
  vector<double> tVector;

  double llhMax= -100.;
  double llhtemp, rms, avgtime;
  
  double sProb, bProb, eMaxRatio;
  I3Event event;
  
  double test[10];
  int tbin;

  NewLlhPeriodicTime* i3an = dynamic_cast<NewLlhPeriodicTime*>(analysisFnVect_[0]);
  
  for (int j=0; j<int(evs.size()); j++) {
  
    event = evs[j];
  
    sProb = event.ProbFrom(*srcCoord_);
    bProb = event.GetBkgSpaceProbFn()->GetBkgProbDensity(event);
    eMaxRatio = event.GetEnergyProbFn()->GetEnergyMaxRatio(event);
    
    tVector.push_back( fmod(event.GetMJD()-t0_,period_)/period_ ); // for sigmamin_
    
    double tTest=0.;
    if ( (eMaxRatio * sProb/bProb) > seedWtMin ) { // make seedWtMin a tunable parameter now.
                                                   // Was set at 10 for IC40 flare analysis,
                                                   // which took events up to ~12 deg away from
                                                   // the source, you should check what makes sense.
                                                   
      tbin = (int) 10.*fmod(event.GetMJD()-t0_,period_)/period_;
      test[tbin] += eMaxRatio * sProb/bProb;
      tTest = fmod(event.GetMJD()-t0_,period_)/period_;
      tVectorclose.push_back( tTest );
      tVectorclose.push_back(tTest+1 ); // add in the event twice at phase+1 to 
                                        // try and help deal with the period better...
    }
  }
  
  int tbinmax=0;
  double max=0;
  for (int j=0;j<10;j++) {
    if (test[j]>max) { max=test[j]; tbinmax=j; }
    if (i3an->monitorLevel_ > 1) { cout << test[j] << " " << flush; }
  } if (i3an->monitorLevel_ > 1) { cout << " : " << tbinmax << endl; }
  
  sort(tVectorclose.begin(),tVectorclose.end());

  sigmamin = 10.;
  double timed;
  int sp = 1;

  for (unsigned int i=sp; i<tVectorclose.size(); i++) { //This assumes events are time-ordered.
    timed = tVectorclose[i] - tVectorclose[i-sp];
    if (fabs(timed) < sigmamin) {sigmamin = fabs(timed);}
  }

  sigmamin /= sqrt(2.);
  if (sigmamin < 5e-3) { sigmamin=5e-3; } // for the periodic analysis, we aren't interested
                                          // in extremely short values.
  
  // I did some checks to see if you might find more flares if we tested a soft
  // spectrum for the seed as well, it didn't seem to be that useful.

  const int ngamma=2;
  double g[ngamma] = {2.0, 3.5};
//  const int ngamma=1;
//  double g[ngamma] = {2.0};

  int kmax = (int) i3an->close_;
  
  if (i3an->monitorLevel_ > 0) { cout << tVectorclose.size() << " events with S/B > " << seedWtMin << ". sigmamin is " << sigmamin << endl; }

  if (10 < kmax) { kmax = 10; }

  int n[] =  {int(tVectorclose.size())/66,
              int(tVectorclose.size())/50, int(tVectorclose.size())/33,
              int(tVectorclose.size())/25, int(tVectorclose.size())/20,
              int(tVectorclose.size())/10, int(tVectorclose.size())/7,
              int(tVectorclose.size())/5,  int(tVectorclose.size())/4,
              int(tVectorclose.size())/3,  int(tVectorclose.size())/2};

  //int n[] =  {1,2,3,4,5,10,15,20,35,tVectorclose.size()/3,tVectorclose.size()/2};
  int ns[] = {1,2,2,2,2, 2, 2, 2, 2, 2, 2};
  int kmin=0;
    
  for (int m=0;m<10;m++) {
    
    avgtime = m/10. + 1./20.;
    rms = 0.2;
  
    llhtemp = EvaluateLlh( 2., 2., avgtime, rms );

    if (llhtemp > llhMax) {
      llhMax = llhtemp;
      Guess_mean = avgtime;
      Guess_rms = rms;
      Guess_nsrc = 2.;
    }
  }

  for (int k=kmin; k<kmax; k++){ // taking a sub-sample of events to test for seed
  
   if ( n[k] > int(tVectorclose.size()) ) continue;
  
    for (unsigned int i=0; i<(tVectorclose.size()-n[k]); i++) {
      for (int j=0;j<=n[k];j++){
        tVectortemp.push_back(tVectorclose[i+j]);
      }
            
   //get mean of the times to seed this Gaussian
      double sum=0, sumsq=0;
      llhtemp = 0.;
      for (unsigned int ii=0;ii<tVectortemp.size();ii++) {
        sum += tVectortemp[ii];
      }
      avgtime = sum/tVectortemp.size();

      for (unsigned int ii=0;ii<tVectortemp.size();ii++) {
        sumsq += (tVectortemp[ii] - avgtime)*(tVectortemp[ii] - avgtime);
      }

    // calculate the rms to use as the sigma
      rms = sqrt( sumsq/tVectortemp.size() );
    
    // This prevents very narrow seeds if we don't use the weighting term
    //  if ( !analysisFnVect_[0]->JimsTerm_ && rms < 0.01 ) {
      if ( rms < 0.005 ) {
        tVectortemp.clear();
        continue;
      }

      for (int h=0; h<ngamma; h++) {
        llhtemp = EvaluateLlh( ns[k]+1., g[h], avgtime, rms );
        if (llhtemp > llhMax) {
          llhMax = llhtemp;
          Guess_mean = avgtime;
          Guess_rms = rms;
          Guess_nsrc = ns[k]+1.;
        }
      }
      tVectortemp.clear();
    }
  }
  
  tVectorclose.clear();
  
  double sllhmax = -100;
  double sllhtemp;
  for (double d=1.; d<4.; d+=0.2) { // loop over gamma with 15 steps for best spectrum seed
    sllhtemp = EvaluateLlh( Guess_nsrc, d, Guess_mean, Guess_rms );
    if (sllhtemp > sllhmax ) {
      Guess_gamma = d;
      sllhmax = sllhtemp;
    }
  }

  double nllhMax = -100;
  for (int m=1;m<40;m++) { // finally a loop over nsrc for the seed value
    
    llhtemp = EvaluateLlh( m/2., Guess_gamma, Guess_mean, Guess_rms );

    if (llhtemp > nllhMax) {
      Guess_nsrc = m/2.;
      nllhMax = llhtemp;
    }
  }

  if ( Guess_mean > 1. ) { Guess_mean -= 1.; }
  
  if (i3an->monitorLevel_ > 0) { printf("Guess Params: %f %f %f %f %f \n",Guess_nsrc,Guess_gamma,Guess_mean,Guess_rms,sigmamin); }
  
}
