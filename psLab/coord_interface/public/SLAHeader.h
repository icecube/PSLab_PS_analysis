#ifndef COORD_INTERFACE_SLAHEADER_H_
#define COORD_INTERFACE_SLAHEADER_H_


// All this header file does is provide a reference for include
// statements in the lab code.  These match exactly the coordinate-service
// declarations in icetray.


namespace SLACoordinateTransform {

#define EPOCH_ 2000

double Local2RA(
		double ang_zenith,
		double ang_azimuth,
		double time_jd,
		double epoch=EPOCH_);

double Local2Dec(
		double ang_zenith,
		double ang_azimuth,
		double time_jd,
		double epoch=EPOCH_);

double Equa2LocalZenith(
		double ra,
		double dec,
		double time_jd,
		double epoch=EPOCH_);

double Equa2LocalAzimuth(
		double ra,
		double dec,
		double time_jd,
		double epoch=EPOCH_);

double Local2GalacticLongitude(
        double ang_zenith,
        double ang_azimuth,
        double time_jd,
        double epoch=EPOCH_);

double Local2GalacticLatitude(
        double ang_zenith,
        double ang_azimuth,
        double time_jd,
        double epoch=EPOCH_);

double Galactic2LocalZenith(
        double gal_longitude,
        double gal_latitude,
        double time,
        double epoch=EPOCH_);

double Galactic2LocalAzimuth(
        double gal_longitude,
        double gal_latitude,
        double time,
        double epoch=EPOCH_);

double Equatorial2GalacticLongitude(
        double right_ascension,
        double declination,
        double time_jd,
        double epoch=EPOCH_);

double Equatorial2GalacticLatitude(
        double right_ascension,
        double declination,
        double time_jd,
        double epoch=EPOCH_);

double Galactic2RA(
        double ga_longitude,
        double ga_latitude,
        double time_jd,
        double epoch=EPOCH_);

double Galactic2Dec(
        double ga_longitude,
        double ga_latitude,
        double time_jd,
        double epoch=EPOCH_);


}

#endif // COORD_INTERFACE_SLAHEADER_H_
