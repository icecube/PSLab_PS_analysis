

double DifferentialToIntegralFluxConstant(double EnLow, double EnReference,
					  double spectralIndex)
{
  return -1./(spectralIndex+1) * pow(EnLow/EnReference, spectralIndex+1);
}


double IntegralToDifferentialFluxConstant(double EnLow, double EnReference, 
					  double spectralIndex)
{
  return -(spectralIndex+1) / pow(EnLow/EnReference, spectralIndex+1);
}
//                     Inf
//     Phi_Int = Integral    Phi_0 * (E/EnReference)^gamma dE
//                     EnLow
//
//     Phi_Int = Phi_0 * (-1)/(gamma+1) * (EnLow/EnReference)^(gamma+1)
//
// ==>   Phi_0 = - (gamma+1) / (EnLow/EnRefernece)^(gamma+1) * Phi_Int
//




double ScaleFlux(double old_EnUnits, double new_EnUnits,
		 double old_EnReference, double new_EnReference,
		 double spectralIndex)
{
  double scaleFactor = (new_EnUnits/old_EnUnits) * 
    pow(new_EnReference/old_EnReference , -fabs(spectralIndex));
  return scaleFactor;
}


// Example:  you have a flux:

//   d Phi                                E
//   ----- = N_0  GeV^-1 cm^-2 s^-1    ( --- )^-3
//    dE                                 GeV

// where the first GeV is the old_EnUnits, 
// and the second GeV is the old_EnReference 
// (i.e. sets the energy scale at which d Phi / dE = N_0 * (1)

// To convert from GeV to TeV in both places, 
// you propagate through all the factors of 0.001:

//   (0.001)^-1 * ( 1 / 0.001 )^ -3  =  1000 * (1000)^-3 = 10^-6

// and that is your scale factor

//   d Phi                                         E
//   ----- = (10^-6 * N_0)   GeV^-1 cm^-2 s^-1  ( --- )^-3
//    dE                                          GeV



// this works if histogram x axis is spectral index, and histogram
// y-axis is old flux

TH1D* RescaleIndexHistogram(double old_EnUnits, double new_EnUnits,
			    double old_EnReference, double new_EnReference,
			    TH1D *h)
{
  TString name = h->GetName();
  TH1D *hnew = dynamic_cast<TH1D*>(h->Clone(name+"_rescaled"));
  for (int i=1; i<= h->GetNbinsX(); ++i) {
    double gamma = h->GetBinCenter(i);
    double scaleFactor = ScaleFlux(old_EnUnits, new_EnUnits,
				   old_EnReference, new_EnReference,
				   gamma);
    hnew->SetBinContent(i, scaleFactor * h->GetBinContent(i));
  }
  return hnew;
}
       



double ConvertFlux(double old_EnUnits, double new_EnUnits,
		   double old_EnReference, double new_EnReference,
		   double spectralIndex, double old_Flux_N)
{
  return old_Flux_N * ScaleFlux(old_EnUnits, new_EnUnits,
				old_EnReference, new_EnReference,
				spectralIndex);
}


