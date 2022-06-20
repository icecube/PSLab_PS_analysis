#include "llh/public/LlhTest.h"

#include "rootExt/public/CountMonitor.h"

#include "llh/public/classes.h"
#include "llh/public/CoordClasses.h"
#include "llh/public/LlhFunctionsBase.h"


bool LlhTest::Ready() {
  if (llh_ && aSet_ && srcCoord_) { return true; }
  cout << "Warning: LlhTest not ready, not completely initialized.\n";
  // log_warn
  return false;
}


void LlhTest::TestRealData() {
  if (!Ready()) { return;}

  aSet_->UseRealData();

  // this loads all of the event likelihood values into llh_
  llh_->SetAnalysis(*aSet_, *srcCoord_);

  llh_->MaximizeLlh();

}


void LlhTest::TestGenerateData_nSrcEvents(int nSrcEvents) {
  if (!Ready()) { return;}

  aSet_->GenerateDataSet_with_nSrcEvents(nSrcEvents);

  // this loads all of the event likelihood values into llh_
  llh_->SetAnalysis(*aSet_, *srcCoord_);

  llh_->MaximizeLlh();

}
  

void LlhTest::Stats_nSrcEvents(int nSrcEvents, int nTrials, 
			       TH1 *hTestStatistic,
			       TH1 *hLogEstProb,
			       TH1 *hPar )
{
  CountMonitor countMon(monitorPercent_, nTrials);

  for (int iTrial=0; iTrial < nTrials; ++iTrial) {

    TestGenerateData_nSrcEvents(nSrcEvents);

    if (hTestStatistic) { hTestStatistic->Fill( llh_->GetTestStatistic());}

    if (hLogEstProb) { hLogEstProb->Fill( log10(llh_->GetEstProb() ));}

    if (hPar) { hPar->Fill( llh_->GetPar(0) ); }

    countMon.UpdateCount();
  }

  if (monitorPercent_) { cout << endl; }
}


  
