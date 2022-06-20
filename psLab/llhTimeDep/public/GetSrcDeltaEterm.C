#include "llhTimeDep/public/GetSrcDeltaEterm.h"

#include "TH1D.h"


double GetSrcDeltaEterm(vector<I3Event>* evVectPtr, double pcmin, double pcmax){
 // The source events and weights need to be set.
 // This returns the ln(emax/emin) where emin and max
 // are by default the 90% containment
 
  TH1D srcHisto; // much easier to clean-up and avoid collisions than new TH1D*
  srcHisto.SetBins(100,2,9);
 
  for (vector<I3Event>::iterator e = evVectPtr->begin();
       e != evVectPtr->end();
       e++) 
  {
    I3MCParameters mc = e->GetMCParams();
    //double weight = fluxModel.GetFlux(mc.mcEnergy) * mc.PS_FlatSpectrumRate; 
    srcHisto.Fill(log10(mc.mcEnergy),mc.srcWeight);
  }

//  return srcHisto;
  
  double a[2] = {pcmin,pcmax};
  double b[2];
  
  srcHisto.GetQuantiles(2,b,a);
  
  cout << " " << b[0] << " " << b[1] << endl;
  
  return log(pow(10.,b[1])/pow(10.,b[0]));
}
