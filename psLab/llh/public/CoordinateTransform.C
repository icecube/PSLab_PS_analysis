#include "llh/public/CoordinateTransform.h"
#include "TMath.h"
#include <star/pal.h>

void LocalToEq(double zenDeg, double azDeg, double timeMJD,
               double& raDeg, double & decDeg)
{

  double azRad  = azDeg  * TMath::DegToRad();
  double zenRad = zenDeg * TMath::DegToRad();

  double PalAz = M_PI/2. - azRad - IC_LONG; //change azimuth zero as required by PAL
  double raRadApp  = NAN;
  double decRadApp = NAN;
  char coordtype[] = {'A'};
  palOap (coordtype,              //type of coordinates – ‘A’ for azimuth/zenith
          PalAz,                  //observed azimuth (radians; azimuth is N=0, E=90∘)
          zenRad,                 //observed zenith (radians)
          timeMJD,                //date/time (Modified Julian Date, JD−2400000.5)
          0,                      //UT1−UTC (UTC seconds)
          IC_LONG,                //IceCube longitude (radians, east +ve)
          IC_LAT,                 //IceCube geodetic latitude (radians)
          0,0,0,0,0,0,0,0,        //very detailed information (humidity, temperature, ...)
          &raRadApp,&decRadApp    //Returned geocentric apparent [α,δ] (radians)
          );

  double raRad  = NAN;
  double decRad = NAN;
  palAmp(raRadApp,decRadApp, // apparent [ α, δ ] (radians)
         timeMJD,            // TDB for apparent place (JD−2400000.5)
         2000.0,             // equinox: Julian epoch of mean place
         &raRad, &decRad     // [ α, δ ] (radians)
         );

  raDeg  = raRad  * TMath::RadToDeg();
  decDeg = decRad * TMath::RadToDeg();
}

void EqToLocal(double raDeg, double decDeg, double timeMJD,
               double& zenDeg, double& azDeg)
{

  double raRad  = raDeg  *TMath::DegToRad();
  double decRad = decDeg *TMath::DegToRad();

  double tt = timeMJD+palDtt(timeMJD)/86400.;
  double appRA  = NAN;
  double appDEC = NAN;
  palMap( raRad, decRad,      // RA, dec (radians) 
          0.0, 0.0, 0.0, 0.0, // additional information (proper motion, parallax, ...)
          2000.0,             //epoch and equinox of data (Julian)
          tt,                 //TDB for apparent place (JD−2400000.5)
          &appRA, &appDEC     //apparent RA, apparent DEC
          );

  double PalAz=NAN;
  double PalZen=NAN;
  double PalHA=NAN;
  double PalDec=NAN;
  double PalRA=NAN;

  // SLA_AOP : Apparent to Observed
  // ACTION  : Apparent to observed place, for sources distant from the solar system.
  palAop(appRA, appDEC,     //geocentric apparent [α,δ] (radians)
           timeMJD,         // date/time (Modified Julian Date, JD−2400000.5)
           0,               //UT1−UTC (UTC seconds)
           IC_LONG,         //IceCube longitude (radians, east +ve)
           IC_LAT,          //IceCube geodetic latitude (radians)
           0,0,0,0,0,0,0,0, //very detailed information (humidity, temperature, ...)
           &PalAz,           //observed azimuth (radians: N=0, E=90∘)
           &PalZen,          //observed zenith (radians)
           &PalHA,           //observed Hour Angle (radians)
           &PalDec,          //observed δ (radians)
           &PalRA            //observed α (radians)
           );

  //convert to IceCube Coordinates
  double azRad  = M_PI/2. - PalAz - IC_LONG;

  //normaize azimuth to [0, 360)
  azRad = fmod(fmod(azRad,2*M_PI)+2*M_PI,2*M_PI);

  assert ( azRad  <= 2*M_PI);
  assert ( azRad  >= 0 );
  assert ( PalZen <= M_PI);
  assert ( PalZen >= 0);

  azDeg  = azRad  * TMath::RadToDeg();
  zenDeg = PalZen * TMath::RadToDeg();
}
