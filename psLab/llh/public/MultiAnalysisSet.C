#include "llh/public/MultiAnalysisSet.h"

#include "rootExt/public/randomfunctions.h"

#include "llh/public/TimePdf.h"
#include "llh/public/I3Event.h"
#include "llh/public/I3SignalGenerator.h"
#include "llh/public/I3Analysis.h"


void MultiAnalysisSet::ConstructEventPtrList() {
  evPtrListMulti_.Clear();
  for (int i=0; i<int(aSetVect_.size()); ++i) {
    const EventPtrList* tempList = aSetVect_[i]->GetEventPtrList();
    for (int j=0; j<tempList->GetSize(); ++j) {
      evPtrListMulti_.AddEvent(tempList->GetEvent(j));
    }
  }
}



void MultiAnalysisSet::ConstructEventPtrListSrc() {
  evPtrListSrcMulti_.Clear();
  for (int i=0; i<int(aSetVect_.size()); ++i) {
    const EventPtrList* tempList = aSetVect_[i]->GetEventPtrListSrc();
    //Printf("SOURCE EVENTS IN SET %d = %d", i, tempList->GetSize());
    for (int j=0; j<tempList->GetSize(); ++j) {
      evPtrListSrcMulti_.AddEvent(tempList->GetEvent(j));
    }
  }
}


void MultiAnalysisSet::ReCollocateEvents() {
  const int nsets = aSetVect_.size();
  EventPtrList tmpEvList[nsets];
  double tmin, tmax;
  //First create temporary lists with correctly recollocated events...
  for (int i=0; i<nsets; ++i) {
    const EventPtrList* evList = aSetVect_[i]->GetEventPtrList();
    for (int j=0; j<evList->GetSize(); ++j) {
      for(int set=0; set<nsets; set++){
        if(aSetVect_[set]->GetSource()->GetTimePdf()){
          tmin = aSetVect_[set]->GetSource()->GetTimePdf()->GetTmin();
          tmax = aSetVect_[set]->GetSource()->GetTimePdf()->GetTmax();
        }
        else if (aSetVect_[set]->GetSource()->GetTimePdfVect().at(0)){
          tmin = aSetVect_[set]->GetSource()->GetTimePdfVect().at(0)->GetTmin();
          tmax = aSetVect_[set]->GetSource()->GetTimePdfVect().at(0)->GetTmax();
        }
        else log_fatal("ERROR: tmin and tmax not set. No time PDF found.\n");
        if(evList->GetEvent(j)->GetTime().GetMJD()>tmin && evList->GetEvent(j)->GetTime().GetMJD()<tmax){
          tmpEvList[set].AddEvent(evList->GetEvent(j));
          break;
        } 
      }
    }
  }

  //...then clear the original lists of events and re-fill them with the temporary lists
  for(int i=0; i<nsets; ++i) {
    const_cast<EventPtrList*>(dynamic_cast<I3Analysis*> (aSetVect_[i])->GetEventPtrList())->Clear();
    for(int j=0; j<tmpEvList[i].GetSize(); ++j){
      const_cast<EventPtrList*>(dynamic_cast<I3Analysis*>(aSetVect_[i])->GetEventPtrList())->AddEvent(tmpEvList[i].GetEvent(j));
    }
  }
}



void MultiAnalysisSet::AddAnalysisSet(AnalysisSet* aSet) { 
  aSetVect_.push_back(aSet); 
  ConstructEventPtrList();
  ConstructEventPtrListSrc();
}


double MultiAnalysisSet::BkgNumberDensity(const Coord& coord) const {
  double sum = 0.;
  for (int i=0; i<int(aSetVect_.size()); ++i) {
    sum += aSetVect_[i]->BkgNumberDensity(coord);
  }
  return sum;
}


double MultiAnalysisSet::GetMeanSrcNev() const {
  double sum = 0.;
  for (int i=0; i<int(aSetVect_.size()); ++i) {
    sum += aSetVect_[i]->GetMeanSrcNev();
  }
  return sum;
}

double MultiAnalysisSet::GetMeanSrcNevTime() const {
  double sum = 0.;
  for (int i=0; i<int(aSetVect_.size()); ++i) {
    if ( aSetVect_[0]->GetSource()->GetTimePdf() ) {
      sum += aSetVect_[i]->GetSource()->GetTimePdf()->GetNorm() * aSetVect_[i]->GetMeanSrcNev();
    }
    else if(aSetVect_[0]->GetSource()->GetTimePdfVect().at(0)){
      for(unsigned int iFl=0; iFl<aSetVect_[0]->GetSource()->GetTimePdfVect().size(); iFl++){
        sum += aSetVect_[i]->GetSource()->GetTimePdfVect().at(iFl)->GetNorm() * aSetVect_[i]->GetMeanSrcNev();
      }
    }
  }        // Woof!
  return sum;
}

double MultiAnalysisSet::GetMeanSrcNevStacking() const {
  double sum = 0.;
  for (int iset=0; iset<int(aSetVect_.size()); ++iset) {
    I3MultiSignalGenerator* multiSigGen = const_cast<I3MultiSignalGenerator*>(dynamic_cast<const I3MultiSignalGenerator*>(aSetVect_[iset]->GetSource()));
    sum += multiSigGen->GetMeanSrcNev();
  }
  return sum;
}

double MultiAnalysisSet::GetFluxTimeIntegralSet(int iset){
  double sum=0;
  for(unsigned int iFl=0; iFl<aSetVect_[0]->GetSource()->GetTimePdfVect().size(); iFl++){
    sum += aSetVect_[iset]->GetSource()->GetTimePdfVect().at(iFl)->GetNorm();
  }
  return sum;
}

double MultiAnalysisSet::GetFluxEnergyIntegralSet(int iset){
  return aSetVect_[iset]->GetMeanSrcNev();
}

double MultiAnalysisSet::
GetMeanSrcNevForFluxModel(const FluxBase& fluxModel) const {
  double sum = 0.;
  for (int i=0; i<int(aSetVect_.size()); ++i) {
    sum += aSetVect_[i]->GetMeanSrcNevForFluxModel(fluxModel);
  }
  return sum;
}


/*double MultiAnalysisSet::GetEnergyQuantile(double prob) const 
{
  int nSets = aSetVect_.size();
  
  int nbins = 100;

  double lowedge = 2.;
  double highedge = 10.;
  
  TH1D * srcHisto = new TH1D("srcHisto","srcHisto", nbins, lowedge, highedge);
  
  for (int i = 0; i < nSets; ++i)
    {
      
      const SourceModule *srcModule = aSetVect_[i]->GetSource();
      
      const I3PointGenerator *i3point = dynamic_cast<const I3PointGenerator*>(srcModule);
	
      
      for (int j=0; j<int(i3point->candidateEvents_.size()); ++j) {
	const I3Event& e = i3point->candidateEvents_[j];
	I3MCParameters mc = e.GetMCParams();
	
	
	if(log10(mc.mcEnergy) < lowedge || log10(mc.mcEnergy) > highedge)
	  {
	    log_warn("Energy out histogram range: %f\n", log10(mc.mcEnergy));
	  }
	
	srcHisto->Fill(log10(mc.mcEnergy),mc.srcWeight);   
      }
      
    }
  
  double quantile;
  
  srcHisto->GetQuantiles(1, &quantile, &prob);
  
  delete srcHisto;
  
  return pow(10., quantile);
    
}


double MultiAnalysisSet::GetEnergyMinSet(int iset){
  int nSets = aSetVect_.size();
  if(iset>(nSets-1)){
    log_error("ERROR: set index too high in MultiAnalysisSet::GetEnergyMinSet");
    return 0;
  }
  const SourceModule *srcModule = aSetVect_[iset]->GetSource();
  const I3PointGenerator *i3point = dynamic_cast<const I3PointGenerator*>(srcModule);
  return i3point->GetEnergyMin();
}

double MultiAnalysisSet::GetEnergyMaxSet(int iset){
  int nSets = aSetVect_.size();
  if(iset>(nSets-1)){
    log_error("ERROR: set index too high in MultiAnalysisSet::GetEnergyMaxSet");
    return 0;
  }
  const SourceModule *srcModule = aSetVect_[iset]->GetSource();
  const I3PointGenerator *i3point = dynamic_cast<const I3PointGenerator*>(srcModule);
  return i3point->GetEnergyMax();
}*/

double MultiAnalysisSet::GetFluxScaleForNev(double nev) const
{
  // Use contribution from a single, non-zero set to determine flux
  // (This assumes that same flux was used for all sets... 
  //  otherwise this question is nonsensical to begin with)

  double totalMeanSrcNev = GetMeanSrcNev();
  int nSets = aSetVect_.size();
  for (int i=0; i<nSets; ++i) {
    double partialMeanSrcNev = aSetVect_[i]->GetMeanSrcNev();

    // check that this contribution is non-zero and non-negligible
    if (partialMeanSrcNev > totalMeanSrcNev/(nSets*2.)) {
      // at least one set has to pass this test, unless all contribute zero

      double thisSetNev = nev * (partialMeanSrcNev / totalMeanSrcNev);
      return aSetVect_[i]->GetFluxScaleForNev(thisSetNev);
    }
  }
  return 0.;  // evidently, sources are not generating any signal events
}


void MultiAnalysisSet::GenerateDataSet() {
  for (int i=0; i<int(aSetVect_.size()); ++i) {
    aSetVect_[i]->GenerateDataSet();
  }
  ConstructEventPtrList();
  ConstructEventPtrListSrc();
}


void MultiAnalysisSet::GenerateDataSet_with_nSrcEvents(int nSrcEvents) {

  const int nSets = aSetVect_.size();
  vector<int> nSrcPerSet(nSets,0); // initialize all values to zero

  double meanSrcNev_Sum = GetMeanSrcNev();
  double meanSrcNev_TimeSum = 0.;
  
  if ( aSetVect_[0]->GetSource()->GetTimePdf() ) {
    meanSrcNev_TimeSum = GetMeanSrcNevTime();  
  }
  // randomly decide which DataSet to generate each signal event in,
  // according to the relative weight of meanSrcEvents from each set  
  if ( aSetVect_[0]->GetSource()->GetTimePdf() ) {
    // This is the parallel way to get the nsrc if there is a TimePdf set (null by default)
    for (int nev=0; nev<nSrcEvents; ++nev) {
      double sumT = 0.;
      double ranNumT = random_uniform(0., meanSrcNev_TimeSum);
      for (int set=0; set<nSets; ++set) {
        sumT += aSetVect_[set]->GetSource()->GetTimePdf()->GetNorm() * aSetVect_[set]->GetMeanSrcNev();
        if (ranNumT <= sumT) {
  	  // this is the dataset which will get this src event
	      nSrcPerSet[set] += 1;
          break;
        }
      }
    }
    if((dynamic_cast<I3Analysis*>(aSetVect_[0]))->GetIsRndOnlyTime()) PrepareTimeOnlyScramble(nSrcEvents, nSrcPerSet, 0);
    for (int set=0; set<nSets; ++set) {aSetVect_[set]->GenerateDataSet_with_nSrcEvents(nSrcPerSet[set]);Printf("ns to be generated in sample %d: %d", set, nSrcPerSet[set]);}
  }
  else if(aSetVect_[0]->GetSource()->GetTimePdfVect().at(0)){
    for(unsigned int iFl=0; iFl<aSetVect_[0]->GetSource()->GetTimePdfVect().size(); iFl++){
      for(int set=0; set<nSets; ++set){
        const_cast<I3SignalGenerator*>(dynamic_cast<const I3SignalGenerator*>(aSetVect_[set]->GetSource()))->SetTimePdf(aSetVect_[set]->GetSource()->GetTimePdfVect().at(iFl));
        nSrcPerSet[set]=0;
        (dynamic_cast<I3Analysis*>(aSetVect_[set]))->SetInjFlare(iFl);
      }
      if ( aSetVect_[0]->GetSource()->GetTimePdf() ) {
        meanSrcNev_TimeSum = GetMeanSrcNevTime();
      }
      // randomly decide which DataSet to generate each signal event in,
      // according to the relative weight of meanSrcEvents from each set  
      if ( aSetVect_[0]->GetSource()->GetTimePdf() ) {
        // This is the parallel way to get the nsrc if there is a TimePdf set (null by default)
        for (int nev=0; nev<nSrcEvents; ++nev) {
          double sumT = 0.;
          double ranNumT = random_uniform(0., meanSrcNev_TimeSum);
          for (int set=0; set<nSets; ++set) {
            sumT += aSetVect_[set]->GetSource()->GetTimePdf()->GetNorm() * aSetVect_[set]->GetMeanSrcNev();
            if (ranNumT <= sumT) {
              // this is the dataset which will get this src event
              nSrcPerSet[set] += 1;
              break;
            }
          }
        }
      }
      if((dynamic_cast<I3Analysis*>(aSetVect_[0]))->GetIsRndOnlyTime()) PrepareTimeOnlyScramble(nSrcEvents, nSrcPerSet, iFl);
      for (int set=0; set<nSets; ++set){
        Printf("\t events to be injected in sample %d = %d", set, nSrcPerSet[set]);
        aSetVect_[set]->GenerateDataSet_with_nSrcEvents(nSrcPerSet[set]);
        const_cast<I3SignalGenerator*>(dynamic_cast<const I3SignalGenerator*>(aSetVect_[set]->GetSource()))->ResetTimePdf();
      }
    } 
  }
  else {
    for (int nev=0; nev<nSrcEvents; ++nev) {
      double sum = 0.;
      double ranNum = random_uniform(0., meanSrcNev_Sum);
      Printf("I'M IN LOOP WITH NO TIME PDF");
      for (int set=0; set<nSets; ++set) {
        sum += aSetVect_[set]->GetMeanSrcNev();
        if (ranNum <= sum) {
	        // this is the dataset which will get this src event
          nSrcPerSet[set] += 1;
	        break;
        }
      }
    }
    if((dynamic_cast<I3Analysis*>(aSetVect_[0]))->GetIsRndOnlyTime()) PrepareTimeOnlyScramble(nSrcEvents, nSrcPerSet, 0);
    for (int set=0; set<nSets; ++set) aSetVect_[set]->GenerateDataSet_with_nSrcEvents(nSrcPerSet[set]);
  } 
  ConstructEventPtrList();
  ConstructEventPtrListSrc();
  if((dynamic_cast<I3Analysis*>(aSetVect_[0]))->GetIsRndOnlyTime()) ReCollocateEvents();
}

void MultiAnalysisSet::GenerateDataSetStacking_with_nSrcEvents(int nSrcEvents) {

  const int nSets = aSetVect_.size();
  vector<int> nSrcPerSet(nSets,0); // initialize all values to zero

  double meanSrcNev_Sum = GetMeanSrcNevStacking();

  for (int nev=0; nev<nSrcEvents; ++nev) {
    double sum = 0.;
    double ranNum = random_uniform(0., meanSrcNev_Sum);
    for (int iset=0; iset<nSets; ++iset) {
      I3MultiSignalGenerator* multiSigGen = const_cast<I3MultiSignalGenerator*>(dynamic_cast<const I3MultiSignalGenerator*>(aSetVect_[iset]->GetSource()));
      sum += multiSigGen->GetMeanSrcNev(); 
      if (ranNum <= sum) {
            // this is the dataset which will get this src event
        nSrcPerSet[iset] += 1;
        break;
      }
    }
  }
 
  for (int set=0; set<nSets; ++set) {Printf("Generating %d events in sample %d", nSrcPerSet[set], set); aSetVect_[set]->GenerateDataSet_with_nSrcEvents(nSrcPerSet[set]);}
  
  ConstructEventPtrList();
  ConstructEventPtrListSrc();
}

void MultiAnalysisSet::PrepareTimeOnlyScramble(int nSrcEvents, vector<int> nSrcPerSet, int iFlare){  
  //if randomizeonlytime: create indeces of events to be injected as signal from tree
  int rndNumber;
  int selEvents = (dynamic_cast<I3Analysis*>(aSetVect_[0]))->GetSelEvNo();
  int nFlares = (dynamic_cast<I3Analysis*>(aSetVect_[0]))->GetNflares();
  int nSets = aSetVect_.size();
  vector<int> sgnIdxVec;
  for (int set=0; set<nSets; ++set) {
    Printf("Sample %d: ns to be generated = %d", set, nSrcPerSet[set]);
    if(set==0 && iFlare==0){
      //create a vector of indeces: they will be used to randomly pick events up from source tree
      for(int j=0; j<selEvents; ++j) sgnIdxVec.push_back(j);
      //remove some indeces: only nSrcEvents will be injected
      for(int j=0; j<selEvents-nFlares*nSrcEvents; ++j){
        rndNumber=random_uniform(0.,sgnIdxVec.size());
        sgnIdxVec.erase(sgnIdxVec.begin()+rndNumber);
      }
      (dynamic_cast<I3Analysis*> (aSetVect_[set]))->SetModifiableSrcIdx(sgnIdxVec);
      (dynamic_cast<I3Analysis*> (aSetVect_[set]))->SetSrcIdx(sgnIdxVec);
    }
    else{
      if(set!=0){
        (dynamic_cast<I3Analysis*> (aSetVect_[set]))->SetModifiableSrcIdx((dynamic_cast<I3Analysis*> (aSetVect_[set-1]))->GetModifiableSrcIdxVec());
        (dynamic_cast<I3Analysis*> (aSetVect_[set]))->SetSrcIdx((dynamic_cast<I3Analysis*> (aSetVect_[set-1]))->GetSrcIdxVec());
      }
      else{
        (dynamic_cast<I3Analysis*> (aSetVect_[set]))->SetModifiableSrcIdx((dynamic_cast<I3Analysis*> (aSetVect_[nSets-1]))->GetModifiableSrcIdxVec());
        (dynamic_cast<I3Analysis*> (aSetVect_[set]))->SetSrcIdx((dynamic_cast<I3Analysis*> (aSetVect_[nSets-1]))->GetSrcIdxVec());
      }
    }
  }
}


void MultiAnalysisSet::UseRealData() {
  for (int i=0; i<int(aSetVect_.size()); ++i) {
    aSetVect_[i]->UseRealData();
  }
  ConstructEventPtrList();
}
