#include "rootExt/public/generalfunctions.h"

#include "rootExt/public/log_report.h"

void sort_array_from_zero(int nentries, double array[])
{
  double *ra = new double [nentries+1];
  for (int i=0; i<nentries; i++)
    ra[i+1] = array[i];
  hpsort(nentries, ra);
  for (int i=0; i<nentries; i++)
    array[i] = ra[i+1];
  delete [] ra;
}


double quantile_from_sorted_array(int nentries, double ra[], double quantile)
{
  return ra[ (int) (quantile*nentries)];
}

//  prob = ( t choose s) p^s (1-p)^(t-s)
//       = exp[ ln t! - ln s! - ln (t-s)! + s ln p + (t-s) ln (1-p) ]

double binomial_prob(int ntotal, int nselected, double p) {
  if (ntotal<nselected || ntotal<0 || nselected<0 || p<0. || p>1.) {
    log_warn("Warning: binomial_prob inputs not valid.\n");
    return 0.;
  }
  if (p==0.) { return (nselected==0)      ? 1. : 0.; }
  if (p==1.) { return (nselected==ntotal) ? 1. : 0.; }
  double lnChooseTerm = factln(ntotal) - factln(nselected)
    - factln(ntotal-nselected);
  double lnProbTerm = nselected*log(p) + (ntotal-nselected)*log(1.-p);
  return exp(lnChooseTerm+lnProbTerm);
}


/* Poisson probability of <relation> n_successes, for given mean */
/* valid char relations: "lt","le","eq","ge","gt"  */

double poisson_prob(double mean, const char* relation, int n_successes)
{
  double logprob;
  double cumulative_prob, decumulative_prob;

  double prob=-1.;  // Makes compiler happy to initialize here

  if (mean>0.) {
    logprob = -mean + n_successes*log(mean) - factln(n_successes);
    prob = exp(logprob);
  } else if (mean==0. && n_successes==0) {
    prob = 1.;
  } else {
    log_fatal("Error: poisson_prob called either with mean (%lg) <0 or with mean==0 and counts (%d) !=0\n",mean,n_successes);
  }

  /* see explanation for poisson function in numerical recipes */
  if (n_successes==0) {
    cumulative_prob = 0.;
    decumulative_prob = 1.;
  } else {
    cumulative_prob = gammq((double) n_successes,mean);
    decumulative_prob = gammp((double) n_successes,mean);
  }

  if (strcmp(relation,"lt")==0)
    return cumulative_prob;
  else if (strcmp(relation,"le")==0)
    return cumulative_prob+prob;
  else if (strcmp(relation,"eq")==0)
    return prob;
  else if (strcmp(relation,"ge")==0)
    return decumulative_prob;
  else if (strcmp(relation,"gt")==0)
    return decumulative_prob-prob;
  else {
    log_fatal("Relation '%s' in poisson() not understood.  Exit.\n",
           relation);
    exit(1); // Makes compiler happy since function should return double.
   }
}



/*  N U M E R I C A L   R E C I P E S   I N   C  */


/* Given dof, returns: p = probability of smaller chisq
                       q = probability of greater chisq    */
void chisq_prob(double chisq, double dof, double *p, double *q)
{
  *p = gammp(dof/2.,chisq/2.);
  *q = gammq(dof/2.,chisq/2.);
}


/* Returns the value of ln[Gamma(xx)] for xx>0. */

double gammln(double xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j <= 5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}


/* Returns ln(n!) */

double factln(int n)
{
  double gammln(double xx);
  static float a[101];/* A static array is automatically initialized to zero. */
  if (n < 0) {
    log_fatal("Negative factorial (%d) in routine factln",n);
  }
  if (n <= 1) return 0.0;
  if (n <= 100) return a[n] ? a[n] : (a[n]=gammln(n+1.0)); 
                            /* In range of table. */
  else return gammln(n+1.0); /* Out of range of table. */
}


/* Returns the incomplete gamma function P(a,x) */

double gammp(double a, double x)
{
	void gcf(double *gammcf, double a, double x, double *gln);
	void gser(double *gamser, double a, double x, double *gln);
	double gamser,gammcf,gln;

	if (x < 0.0 || a <= 0.0) {
	  log_fatal("Invalid arguments in routine gammp:\n(x=%lg)<0.0 or (a=%lg)<=0.0\n",x,a);
	}

	if (x < (a+1.0)) {
		gser(&gamser,a,x,&gln);
		return gamser;
	} else {
		gcf(&gammcf,a,x,&gln);
		return 1.0-gammcf;
	}
}


/* Returns the incomplete gamma function Q(a,x) */

double gammq(double a, double x)
{
	void gcf(double *gammcf, double a, double x, double *gln);
	void gser(double *gamser, double a, double x, double *gln);
	double gamser,gammcf,gln;

	if (x < 0.0 || a <= 0.0) {
	  log_fatal("Invalid arguments in routine gammq:\n(x=%lg)<0.0 or (a=%lg)<=0.0\n",x,a);
	}
	if (x < (a+1.0)) {
		gser(&gamser,a,x,&gln);
		return 1.0-gamser;
	} else {
		gcf(&gammcf,a,x,&gln);
		return gammcf;
	}
}


/* INCOMPLETE GAMMA FUNCTIONS FOR gammp and gammq ROUTINES */
// Raising ITMAX from default 100 to 1000 // cbf 2007/04/19
// (typically, need at most ~ few * sqrt(a) iterations, and only for
// when x is close to a.  ITMAX of 1000 is good for a~50,000.


#define ITMAX 1000
#define EPS 3.0e-7
void gser(double *gamser, double a, double x, double *gln)
{
	double gammln(double xx);
	int n;
	double sum,del,ap;

	*gln=gammln(a);
	if (x <= 0.0) {
		if (x < 0.0) {
		  log_fatal("x (%lg) < 0.0 in routine gser",x);
		}
		*gamser=0.0;
		return;
	} else {
		ap=a;
		del=sum=1.0/a;
		for (n=1;n<=ITMAX;n++) {
			++ap;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS) {
				*gamser=sum*exp(-x+a*log(x)-(*gln));
				return;
			}
		}
		log_fatal("a (%lg) too large, ITMAX (%d) too small in routine gser\n",a,ITMAX);
	}
}
#undef ITMAX
#undef EPS

// Raising ITMAX from default 100 to 1000 // cbf 2007/04/19
#define ITMAX 1000
#define EPS 3.0e-7
#define FPMIN 1.0e-30
void gcf(double *gammcf, double a, double x, double *gln)
{
	double gammln(double xx);
	int i;
	double an,b,c,d,del,h;

	*gln=gammln(a);
	b=x+1.0-a;
	c=1.0/FPMIN;
	d=1.0/b;
	h=d;
	for (i=1;i<=ITMAX;i++) {
		an = -i*(i-a);
		b += 2.0;
		d=an*d+b;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=b+an/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS) break;
	}
	if (i > ITMAX) {
	  log_fatal("a (%lg) too large, ITMAX (%d) too small in gcf\n",a,ITMAX);
        }
	*gammcf=exp(-x+a*log(x)-(*gln))*h;
}
#undef ITMAX
#undef EPS
#undef FPMIN



/* A HEAPSORT ROUTINE */

/* returns ra[] sorted in ascending order */

/* ! ! ! NOTE: array must start at 1, hpsort skips 0 !!! */

void hpsort(unsigned long n, double ra[])
{
  unsigned long i, ir, j, l;
  double rra;

  if (n < 2) return;
  l=(n >> 1) +1;
  ir=n;

  for (;;) {
    if (l > 1) {
      rra=ra[--l];
    } else {
      rra=ra[ir];
      ra[ir]=ra[1];
      if (--ir == 1) {
        ra[1]=rra;
        break;
      }
    }
    i=l;
    j=l+l;
    while (j <= ir) {
      if (j < ir && ra[j] < ra[j+1]) j++;
      if (rra < ra[j]) {
        ra[i]=ra[j];
        i=j;
        j <<= 1;
      } else break;
    }
    ra[i]=rra;
  }
}



/* Numerical Recipes standard error handler */

void nrerror(char error_text[])
{
  fprintf(stderr,"Numerical Recipes run-time error...\n");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"...now exiting to system...\n");
  //  exit(1);
  log_fatal("\n");
}
