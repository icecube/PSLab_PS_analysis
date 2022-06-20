// Variables declared outside the functions will have global scope
// and can be accessed after function exit.
//
// Note that if they are declared again (either by re-running this 
// script or another one which declares the same variables)
// their values may be reset, and should not be relied upon.



// Default macro, gives user instructions

void macro_SetPointSource_PowerLaw() {
  cout << "Usage:\n";
  cout << ".x macro_SetPointSource_PowerLaw.C(double raDeg, double decDeg,\n";
  cout << "                             double index, double fluxConstant=1.)";
  cout << "\n or:\n";
  cout << ".x macro_SetPointSource_PowerLaw.C(EquatorialDeg srcCoord,\n";
  cout << "                             double index, double fluxConstant=1.)";

  cout << "\n";
  cout << "Create source with a power-law and set psData with this source.\n";
}


void macro_SetPointSource_PowerLaw(double raDeg, double decDeg, double index,
		       double fluxConstant = 1.)
{
  EquatorialDeg srcCoord(raDeg,decDeg);
  macro_SetPointSource_PowerLaw(srcCoord, index, fluxConstant);
}


void macro_SetPointSource_PowerLaw(EquatorialDeg srcCoord, double index,
		       double fluxConstant = 1.)
{
  cout << "macro_SetPointSource_PowerLaw:\n";

  // Must have already defined these elsewhere
  extern EventLoader evLoader;
  extern double livetime;
  extern I3Analysis psData;
  extern mySrcLocation;

  mySrcLocation = srcCoord;

  cout << "!! mySrcLocation set to:  " << mySrcLocation.GetRa()<<" r.a., ";
  cout << mySrcLocation.GetDec() << " dec.\n";
  cout << "    Index: " << index;
  cout << "    FluxConstant: " << fluxConstant << "  (GeV^-1 cm^-2 s^-1)\n";

  vector<I3Event> sourceEvents;
  evLoader.LoadSourceEvents(sourceEvents, mySrcLocation);

  PowerLawFlux flux(fluxConstant, index); 

  i3point= I3PointGenerator(sourceEvents, flux, mySrcLocation, livetime);

  psData.SetSource(i3point);
  cout << "!! psData was set with the source; accessible via: psData.GetSource()\n";
}
