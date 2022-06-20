
#include "llhTimeDep/public/MultiBlockAnalysisFn.h"
#include "iostream"

void llhFncMSBlock(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  f=0;
  double result = 0.;
  
  vector<double> parVect;
  for(int i=0; i<npar; i++) parVect.push_back(par[i]);

  for (int i=0; i<int(tdMLlhBlock->GetAnalysisFn().size()); ++i) {
    vector<double> individualParVect = tdMLlhBlock->GetParTranslator()->Translate(i, parVect);
    result = tdMLlhBlock->GetAnalysisFn().at(i)->EvalFCN(individualParVect);
    f += result;
  }
  return;
}

MultiBlockAnalysisFn::MultiBlockAnalysisFn() {
  tdMLlhBlock = this;

  nPar = 4;
  laglimit_ = 0.5;
  Ndof=-1; 
}

void MultiBlockAnalysisFn::MaximizeLlh() {

    minuit_ = new TMinuit(nPar);
    Double_t arglist[1];
    arglist[0] = 0.5; //0.5 for likelihood fit
    Int_t ierflg=0;
    minuit_->SetFCN(llhFncMSBlock);
    minuit_->mnexcm("SET ERR",arglist,1,ierflg);
    minuit_->SetPrintLevel(-1);

    PrepareAnalysis();

    parDefVect_.clear();
    parDefVect_.push_back( MinuitParDef("nSrc",2.5,1., 0.,100.) );
    parDefVect_.push_back( MinuitParDef("gamma",2.5,0.1, 1., 4.) );

    double lagGuess = SearchForLag(laglimit_);
    parDefVect_.push_back( MinuitParDef("lag", lagGuess, laglimit_/10., -1.0*laglimit_, laglimit_) );

    double threshMax = ( GetHighestBlock(blocksFile_) +
                         GetSecondHighestBlock(blocksFile_) ) /2.;

    double guessThresh;
    SearchBlockSpace(blocksFile_, lagGuess, guessThresh, threshMax);
    parDefVect_.push_back( MinuitParDef("thresh",guessThresh, threshMax/20., 0., threshMax) );

    for (int i=0; i<int(parDefVect_.size()); ++i) {

      const MinuitParDef& pd = parDefVect_[i];
      minuit_->mnparm(i, pd.name.c_str(), pd.initValue, pd.initStepSize,
                            pd.lowLimit, pd.upLimit, ierflg);
      //cout << pd.name.c_str() << " Init: " << pd.initValue << " max: " << pd.upLimit << endl;
    }

    double arglist2[2];
    arglist2[0] = 500;
    arglist2[1] = 0.1;
    minuit_->mnexcm("MIGRAD", arglist2, 2, minuitOut_);

    StoreLogLambdaBest();
  }

double MultiBlockAnalysisFn::EvalFCN(const vector<double>& parVect) const {
  double f;
  const int npar = parVect.size();
  double gin[npar];
  double par[npar];
  int iflag = 0;
  for(int i=0; i<npar; i++) par[i] = parVect[i];
  int npars=npar;
  llhFncMSBlock(npars, gin, f, par, iflag);

  return f;
}

double MultiBlockAnalysisFn::EvaluateLlh(double *parValueArray) {
  vector<double> parVect;

  //for (int i=0; i<int(parDefVect_.size()); ++i) {
  for (int i=0; i<nPar; ++i) {
    parVect.push_back(parValueArray[i]);
  }
  double minusLlh = EvalFCN(parVect);
  return -minusLlh;   // that is, max llh = - (minimizer result)
}


void MultiBlockAnalysisFn::StoreLogLambdaBest()
{
  Double_t amin, edm, errdef;
  Int_t nvpar, nparx, istat;
  minuit_->mnstat(amin, edm, errdef, nvpar, nparx,istat);;
  logLambdaBest_ = -amin;  // that is, max llh = - (minimizer result)

  // The *worst* logLambdaBest should be zero (i.e. null hypothesis).
  // But minimizer will miss exact zero, leading logLambdaBest_ to be slightly
  // negative. Since this will cause probability calculation to choke, we
  // fix it here
  if (logLambdaBest_ < 0.) { logLambdaBest_ = 0.;}
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

vector<I3Event> MultiBlockAnalysisFn::GetAllEvents() {
	    vector<I3Event> allEvents;
	    vector<I3Event> tempVect;

	    for (int i=0; i<int(analysisFnVect_.size()); ++i) {

	        NewLlhBlockTime* i3an = dynamic_cast<NewLlhBlockTime*>(analysisFnVect_[i]);

	        if (i3an->Get_nEvents()==0) {i3an->PrepareAnalysis(); }

	        tempVect.clear();
	        tempVect = i3an->GetEventVector();

	        for (int j=0; j<int(tempVect.size()); j++) {
	        allEvents.push_back(tempVect[j]);
	        }
	    }
	    return allEvents;
 }


vector<I3Event> MultiBlockAnalysisFn::GetAllSrcEvents() {
	    vector<I3Event> allSrcEvents;
	    vector<I3Event> tempVect;

	    for (int i=0; i<int(analysisFnVect_.size()); ++i) {

	        NewLlhBlockTime* i3an = dynamic_cast<NewLlhBlockTime*>(analysisFnVect_[i]);

	        if (i3an->Get_nEvents()==0) {i3an->PrepareAnalysis(); }

	        tempVect.clear();
	        tempVect = i3an->GetEventVectorSrc();

	        for (int j=0; j<int(tempVect.size()); j++) {
	        allSrcEvents.push_back(tempVect[j]);
	        }
	    }
	    return allSrcEvents;
 }


/*
// This function scans over all the events which could have an S/B ratio > 10,
// and tests consecutive doubles, triples... for compatability with and E^-2 flare.
// I tried testing for Gamma = 2 and like 3.5 before, but it doesn't seem to do much.
// This is configurable enough to do whatever you'd like.
// There is a final scan over Gamma at the end to round out the set parameters.

// Oh, I'm using this parameter close_ to set how far into the vector 1,2,3,4,5,10,15,20,25
// to test for flare compatability. You only need to go to 10, I usually go to 15 (close_=6).

void MultiBlockAnalysisFn::GetFlareGuess(double & Guess_nsrc, double & Guess_gamma, double & Guess_mean, double & Guess_rms) {


  // Note: already have Coord * srcCoord_ set;

  vector<I3Event> evs = GetAllEvents();

  vector<double> tVectortemp;
  vector<double> tVectorclose;

  double llhMax= -100.;
  double llhtemp, rms, avgtime;

  double sProb, eMaxRatio;
  const EnergyProb* eProb(NULL);
  I3Event event;

  for (int j=0; j<int(evs.size()); j++) {

    event = evs[j];

    sProb = event->ProbFrom(*srcCoord_);
    bProb = event->GetBkgSpaceProbFn()->GetBkgProbDensity(*event); //HMMM
    eProb = event->GetEnergyProbFn();
    eMaxRatio = eProb->GetEnergyMaxRatio(*event);


    if ( (eMaxRatio * sProb/bProb) > 10. ) {
        tVectorclose.push_back(event.GetMJD());
    }
  }

  sort(tVectorclose.begin(),tVectorclose.end());

  // I did some checks to see if you might find more flares if we tested a soft
  // spectrum for the seed as well, it didn't seem to be that useful.

//  const int ngamma=2;
//  double g[ngamma] = {2.0, 3.5};
  const int ngamma=1;
  double g[ngamma] = {2.0};

  int kmax = (int) close_;
  if (8 < kmax) { kmax = 8; }

  int n[] =  {1,2,3,4,5,10,15,20,25};
  int ns[] = {1,2,3,4,5, 5, 5, 5, 5};
  double p[4];

  int kmin=0;

  while ( tVectorclose.size()*1.0 < n[kmax] ) { kmax--; }

  for (int k=kmin; k<kmax; k++){ //k==1 for pairs, k==2 for triples, etc.

    for (unsigned int i=0; i<(tVectorclose.size()-n[k]); i++) {
      for (int j=0;j<=n[k];j++){
        tVectortemp.push_back(tVectorclose[i+j]);
      }

   //get mean of the times to seed this gaussian
      double sum=0, sumsq=0;
      llhtemp = 0.;
      for (unsigned int ii=0;ii<tVectortemp.size();ii++) {
        sum += tVectortemp[ii];
      }
      avgtime = sum/tVectortemp.size();

      for (unsigned int ii=0;ii<tVectortemp.size();ii++) {
        sumsq += (tVectortemp[ii] - avgtime)*(tVectortemp[ii] - avgtime);
      }

    //calculate the rms to use as the sigma
      rms = sqrt( sumsq/tVectortemp.size() );

//      if (useE) {
        for (int h=0; h<ngamma; h++) {
          p[0] = ns[k]+1.;
          p[1] = g[h];
          p[2] = avgtime;
          p[3] = log10(rms);
          llhtemp = EvaluateLlh( p );
          if (llhtemp > llhMax) {
            llhMax = llhtemp;
            Guess_mean = avgtime;
            Guess_rms = rms;
            Guess_nsrc = ns[k]+1.;
          }
        }
        tVectortemp.clear();
//        } else { // I have no idea why you'd want to do the analysis without energy, but...
//         llhtemp = EvaluateLlh( ns[k]+1., 0., avgtime, log10(rms));
//         if (llhtemp > llhMax) {
//           llhMax = llhtemp;
//           Guess_mean = avgtime;
//           Guess_rms = rms;
//           Guess_nsrc = ns[k]+1.;
//         }
//         tVectortemp.clear();
//       }
    }
  }

  tVectorclose.clear();

  double sllhmax = -100;
  double sllhtemp;
  for (double d=1.; d<4.; d+=0.2) { // loop over gamma with 15 steps for best seed
    sllhtemp = EvaluateLlh( Guess_nsrc, d, Guess_mean, log10(Guess_rms));
    if (sllhtemp > sllhmax ) {
      Guess_gamma = d;
      sllhmax = sllhtemp;
    }
  }
  } */

// Here are two quick functions for doing a lightcurve-based timedep
// analysis on multiple years/datasets. The initial conditions need to be
// calculated with all the data, the first one calculates a best possible
// lag btw photons and neutrinos and the second gets the max height of the
// lc and initial seed.


double MultiBlockAnalysisFn::SearchForLag(double laglimit) {

  double llhMax=-100.;
  double lagb=0., llhtemp;

  double p[] = {2, 2., 0., 0.};

  for (double d=-1.0*laglimit; d<laglimit; d=d+laglimit/10.){
    p[2] = d;
    llhtemp = EvaluateLlh( p );
    if (llhtemp > llhMax) {
      llhMax = llhtemp;
      lagb = d;
    }
  }

  return lagb;
}

void MultiBlockAnalysisFn::SearchBlockSpace(string blocksFile, double lag, double & initV, double & maxT) {

  //string BlocksTimeFile_ = analysisFnVect_[0]->GetBlocksFile();

  //maxT = (GetHighestBlock(blocksFile.c_str()) + GetSecondHighestBlock(blocksFile_.c_str()))/2.;

  ifstream fin;
  fin.open(blocksFile.c_str());
  double blockbegin, blockdur, blocklev;
  double max=-1;
  double secondmax=-2;
  while (fin >> blockbegin) {
    fin >> blocklev >> blockdur;
    if (blocklev > secondmax) {
      if (blocklev > max) {
        secondmax = max;
        max = blocklev;
      } else {
        secondmax = blocklev;
      }
    }
  }
  fin.close();

  maxT = (max + secondmax) / 2.;

  double step = maxT/20.;
  double llhMax=-100.;
  double llhtemp;

  double p[] = {2, 2., lag, 0.};

  for (double d=0.;d<maxT;d+=step) {
    p[3] = d;
    llhtemp = EvaluateLlh( p );
    if (llhtemp > llhMax) {
      llhMax = llhtemp;
      initV = d;
    }
  }

} //*/
