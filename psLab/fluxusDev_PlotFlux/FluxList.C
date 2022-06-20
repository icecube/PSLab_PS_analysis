 
{
  FormulaFlux fluxMGRO_J1908_old("1.2e-9*pow(x,-2)*exp(-sqrt(x/3e5))");
  EquatorialDeg srcMGRO_J1908_old( 287.27 , 6.28 );
  
  FormulaFlux fluxCrab("3.01e-7*pow(x,-2.39)*exp(-x/7e3)");
  EquatorialDeg srcCrab( 83.633 , 22.014 );

  FormulaFlux flux3C279("1e-20*pow(x/1e6,-.25)*(x<1e6)+1e-8*pow(x,-2)*(x>=1e6 && x<3e8)");
  EquatorialDeg src3C279( 194.1 , -5.8 );
}
