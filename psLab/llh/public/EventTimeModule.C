#include "llh/public/EventTimeModule.h"

#include "rootExt/public/log_report.h"
#include "rootExt/public/randomfunctions.h"

#include "llh/public/CoordinateTransform.h"


// CONSTRUCTORS //


EventTimeModule::EventTimeModule() {
  // arbitrary day, 2008/02/04
  //  double mjd1 = 54500.;

  // arbitrary, start of J2000 epoch
  double mjd1 = 51544.5;

  // add one *sidereal* day, so that distribution of random times will be
  // exactly uniform in right ascension
  // (whereas a full earth day would have a slight overlap in r.a.)
  double mjd2 = mjd1 + LengthOfSiderealDayInSeconds()/86400.;
  
  SetTimesWithinRange(Time(mjd1), Time(mjd2));
}


// FUNCTIONS


void EventTimeModule::SetTimesWithinRange(Time timeMin, Time timeMax) {
  timeMin_ = timeMin;
  timeMax_ = timeMax;
}


Time EventTimeModule::GetRandomTime() {
  if ( ! CompareTime(timeMin_, timeMax_) ) {
    log_fatal("Error: timeMin_ not less than timeMax_ in EventTimeModule::"
	      "GetRandomTime.\n");
    exit(1);
  }
  double mjd = random_uniform( timeMin_.GetMJD(), timeMax_.GetMJD() );
  return Time(mjd);
}

//arrival time of event is taken from a Gaussian centered on the source flaring time
Time EventTimeModule::GetGaussianTime(double mean, double sigma) {

  if ( ! CompareTime(timeMin_, timeMax_) ) {
    log_fatal("Error: timeMin_ not less than timeMax_ in EventTimeModule::"
	      "GetRandomTime.\n");
    exit(1);
  }
  
  double mjd = random_gaussian(mean, sigma);
  return Time(mjd);
}


// We assume in this module that dec doesn't change with time, which
// ignores precesion... should deal with this later in a more careful way,
// but for *scrambled* events this should not be so important

void EventTimeModule::RandomizeEventGaussian(I3Event& ev, double mean, double sigma) {
  double azDeg = ev.GetParams().recoAzimuthDeg;

  // calculate the r.a. value for a new random time (implem. in derived class)
  double mjd = GetGaussianTime(mean, sigma).GetMJD();
  
  // use some magic numbers as a fast short-cut here instead of full
  // Coordinate Transformation
  double raDeg = 145.6453 + fmod(mjd, 0.997269566)/0.997269566*360. - azDeg;
  raDeg = fmod(raDeg, 360.);
  if (raDeg < 0.) { raDeg += 360.; }

  ev.SetMJD(mjd);
  ev.SetRaDeg(raDeg);
}

void EventTimeModule::RandomizeEvent(I3Event& ev, bool rndOnlyTimes) {

  double mjd = GetRandomTime().GetMJD();
  ev.SetMJD(mjd);
  if(!rndOnlyTimes) {
    double azDeg = ev.GetParams().recoAzimuthDeg;
    // calculate the r.a. value for a new random time (implem. in derived class)
    // use some magic numbers as a fast short-cut here instead of full
    // Coordinate Transformation
    double raDeg = 145.6453 + fmod(mjd, 0.997269566)/0.997269566*360. - azDeg;
    raDeg = fmod(raDeg, 360.);
    if (raDeg < 0.) { raDeg += 360.; }
    ev.SetRaDeg(raDeg);
  }  
}
