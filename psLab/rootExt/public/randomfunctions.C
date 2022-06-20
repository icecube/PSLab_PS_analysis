#include "rootExt/public/randomfunctions.h"

#include "rootExt/public/generalfunctions.h"
// this is needed for poidev



// static external variable 'my_idum' is defined here so user of these
// functions does not have to know about it and scope is restricted to
// subsequent functions in this call.

static long my_idum;
static bool ran1_init = false;

static long my_seed = 0;  // for look-back purposes, if you want to reproduce something


/*** CALL FUNCTIONS ***/


/* i.e. the functions you will actually call */


double random_uniform(double min, double max)
{
  if (!ran1_init)
    initialize_ran1();
  return min+(max-min)*ran1(&my_idum);
}

int random_poisson(double mean)
{
  if (!ran1_init)
    initialize_ran1();
  return (int) (0.5+poidev(mean, &my_idum));
  // poidev actually returns a double, so we carefully convert to integer
}

double random_gaussian(double mean, double sigma)
{
  if (!ran1_init)
    initialize_ran1();
  return mean+sigma*gasdev(&my_idum);
}



/*** MISC FUNCTIONS ***/


void initialize_ran1()
{
  srand((unsigned)time(NULL));
  initialize_ran1( -rand() );
}


void initialize_ran1(long ranSeed) {
  // Better to force quit on seeding error than let a batch job
  // generate non-randomized data

  //TO BE DECOMMENTED!!!
  /*if ( ran1_init) {
    printf("Error: ran1 was already initialized.\n"
	   "Re-initialization not implemented.\n");
    exit(1);
  }*/

  if (ranSeed>=0) {
    printf("Error: random seed must be negative (non-zero) integer.\n");
    exit(1);
  }

  my_seed = ranSeed;
  my_idum = my_seed;
  ran1(&my_idum);
  printf("\nInitializing ran1 with random seed = %ld\n",my_seed);
  ran1_init = true;
}



long get_ran1_seed() {
  if (my_seed==0) { printf("ran1 not yet initialized.\n"); }
  return my_seed;
}



/*** WORK FUNCTIONS ***/

/* i.e. the functions which actually do the work */


/*   RANDOM NUMBER GENERATOR ran1()  */

/* Recommended random number generator from Numerical Recipes in C
   Intialize with initialize_ran1(&idum)
   then use ran1(&idum) repeatedly to get random number between 0 and 1
   For more than 10^8 random numbers in one calculation, 
   it is recommended to use ran2() */

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
float ran1(long *idum)
{
  int j;
  long k;
  static long iy=0;
  static long iv[NTAB];
  float temp;

  if (*idum <=0 || !iy) {
    if (-(*idum) < 1) *idum = 1;
    else *idum = -(*idum);
    for (j=NTAB+7; j>=0; j--) {
      k=(*idum)/IQ;
      *idum=IA*(*idum-k*IQ)-IR*k;
      if (*idum <0) *idum += IM;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ;
  *idum=IA*(*idum-k*IQ)-IR*k;
  if (*idum < 0) * idum += IM;
  j=iy/NDIV;
  iy=iv[j];
  iv[j] = *idum;
  if ((temp=AM*iy) > RNMX) return RNMX;
  else return temp;
} 
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX



/*** POISSON RANDOM DEVIATE ***/

/* From Numerical Recipes in C
   Returns as a double precision number an integer value that is a random 
   deviate drawn from a Poisson distribution of mean xm, using ran1(idum) 
   as a source of uniform random deviates. */

#define PI 3.14159265358979312
double poidev(double xm, long *idum)
{
  double gammln(double xx);
  float ran1(long *idum);
  static double sq,alxm,g,oldm=(-1.0); /* oldm is a flag for whether xm has */ 
  double em,t,y;                       /* changed since last  call. */

  if (xm < 12.0) {         /* Use direct method. */
    if (xm != oldm) {
      oldm=xm;
      g=exp(-xm);          /* If xm is new, compute the exponential. */
    }
    em = -1;
    t=1.0;
    do {         /* Instead of adding exponential deviates it is equivalent
		    to multiply uniform deviates. We never
		    actually have to take the log, merely compare
		    to the pre-computed exponential. */
      ++em;
      t *= ran1(idum);
    } while (t > g);
  } else {      /* Use rejection method. */
    if (xm != oldm) { /* If xm has changed since the last call, then pre- */
      oldm=xm;        /* compute some functions that occur below. */
      sq=sqrt(2.0*xm);
      alxm=log(xm);
      g=xm*alxm-gammln(xm+1.0);
      /* The function gammln is the natural log of the gamma function, 
	 as given in §6.1. */
    }
    do {
      do {   /* y is a deviate from a Lorentzian comparison function. */ 
	y=tan(PI*ran1(idum));
	em=sq*y+xm;            /* em is y, shifted and scaled. */
      } while (em < 0.0);      /* Reject if in regime of zero probability. */
      em=floor(em);         /* The trick for integer-valued distributions. */
      t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
      /* The ratio of the desired distribution to the comparison function; we 
	 accept or reject by comparing it to another uniform deviate. 
	 The factor 0.9 is chosen so that t never exceeds 1. */
    } while (ran1(idum) > t);
  }
  return em;
}
#undef PI


/*** GAUSSIAN NORMAL DEVIATE ***/

/* From Numerical Recipes In C */
/* Returns a normally distributed deviate with zero mean and unit variance,
   using ran1(idum) as the source of uniform deviates. */

double gasdev(long *idum)
{
  float ran1(long *idum);
  static int iset=0;
  static double gset;
  double fac,rsq,v1,v2;

  if (*idum < 0) iset=0;  /* Reinitialize */
  if (iset == 0) {
    do {
      v1=2.0*ran1(idum)-1.0;
      v2=2.0*ran1(idum)-1.0;
      rsq=v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac=sqrt(-2.0*log(rsq)/rsq);
    /* Now make the Box-Muller transformation to get two normal
       deviates.  Return one and save the other for next time. */
    gset=v1*fac;
    iset=1;
    return v2*fac;
  } else {
    iset=0;
    return gset;
  }
}
