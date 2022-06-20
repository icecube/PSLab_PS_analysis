// Variables declared outside the functions will have global scope
// and can be accessed after function exit.
//
// Note that if they are declared again (either by re-running this 
// script or another one which declares the same variables)
// their values may be reset, and should not be relied upon.



// Default macro, gives user instructions

void macro_SetPointSource_FormulaFlux() {
  cout<<"Create source with formula flux and set psData with this source.\n";
  cout<<" Usage:\n";
  cout<<".x macro_SetPointSource_FormulaFlux.C(double raDeg, double decDeg,\n";
  cout<<"                                      const char* formula)\n";
  cout<<" or:\n";
  cout<<".x macro_SetPointSource_FormulaFlux.C(double raDeg, double decDeg,\n";
  cout<<"                                      const FormulaFlux& flux)\n";
}


void macro_SetPointSource_FormulaFlux(double raDeg, double decDeg,
				      const char* formula)
{
  FormulaFlux flux(formula);
  macro_SetPointSource_FormulaFlux(raDeg, decDeg, flux);
}


void macro_SetPointSource_FormulaFlux(double raDeg, double decDeg,
				      const FormulaFlux& flux)
{
  cout << "macro_SetPointSource_FormulaFlux:\n";

  // Must have already defined these elsewhere
  extern EventLoader evLoader;
  extern double livetime;
  extern I3Analysis psData;

  extern mySrcLocation = EquatorialDeg(raDeg, decDeg);

  cout << "  mySrcLocation set to:  "<<mySrcLocation.GetRa()<<" r.a., ";
  cout << mySrcLocation.GetDec() << " dec.\n";
  cout << "  Formula: " << flux.GetFormula() << endl;
  cout << "  where x is GeV, and flux is GeV^-1 cm^-2 s^-1\n";

  vector<I3Event> sourceEvents;
  evLoader.LoadSourceEvents(sourceEvents, mySrcLocation);

  i3point= I3PointGenerator(sourceEvents, flux, mySrcLocation, livetime);

  psData.SetSource(i3point);
  cout << "Source assigned to psData; accessible via: psData.GetSource()\n";
}
