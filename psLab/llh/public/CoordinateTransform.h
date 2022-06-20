#ifndef LLH_COORDINATE_TRANSFORM_H_
#define LLH_COORDINATE_TRANSFORM_H_

#define IC_LONG -63.453056*TMath::DegToRad()
#define IC_LAT -89.99*TMath::DegToRad()
#define LENGTH_OF_SIDEREAL_DAY_SECONDS 86164.128

void LocalToEq(double zenDeg, double azDeg, double timeMJD,
	       double& raDeg, double & decDeg);

void EqToLocal(double raDeg, double decDeg, double timeMJD,
	       double& zenDeg, double& azDeg);

double LengthOfSiderealDayInSeconds() { return LENGTH_OF_SIDEREAL_DAY_SECONDS; }

#endif // LLH_COORDINATE_TRANSFORM_H_
