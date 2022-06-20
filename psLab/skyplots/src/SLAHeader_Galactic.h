#ifndef SKYPLOT_SLAHEADER_GALACTIC_H_
#define SKYPLOT_SLAHEADER_GALACTIC_H_


namespace SLACoordinateTransform {

#define EPOCH_ 2000

double Equatorial2GalacticLongitude(
				    double right_ascension, 
				    double declination, 
				    double time_jd, 
				    double epoch = EPOCH_);

double Equatorial2GalacticLatitude(
				    double right_ascension, 
				    double declination, 
				    double time_jd, 
				    double epoch = EPOCH_);

double Galactic2RA(
		   double ga_longitude, 
		   double ga_latitude, 
		   double time_jd, 
		   double epoch = EPOCH_);

double Galactic2Dec(
		   double ga_longitude, 
		   double ga_latitude, 
		   double time_jd, 
		   double epoch = EPOCH_);
}

#endif // SKYPLOT_SLAHEADER_GALACTIC_H_
