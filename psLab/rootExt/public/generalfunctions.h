#ifndef ROOTEXT_GENERAL_FUNCTIONS_H_
#define ROOTEXT_GENERAL_FUNCTIONS_H_

double binomial_prob(int ntotal, int nselected, double p);

void sort_array_from_zero(int nentries, double array[]);

double quantile_from_sorted_array(int nentries, double ra[], double quantile);

double poisson_prob(double mean, const char* relation, int n_successes);


/*  N U M E R I C A L   R E C I P E S   I N   C  */

void chisq_prob(double chisq, double dof, double *p, double *q);

/* RETURNS THE VALUE OF ln[Gamma(xx)] for xx>0. */
double gammln(double xx);

/* Returns ln(n!) */
double factln(int n);

/* Returns the incomplete gamma function P(a,x) */
double gammp(double a, double x);

/* Returns the incomplete gamma function Q(a,x) */
double gammq(double a, double x);

/* INCOMPLETE GAMMA FUNCTIONS FOR gammp and gammq ROUTINES */
void gser(double *gamser, double a, double x, double *gln);
void gcf(double *gammcf, double a, double x, double *gln);

/* HEAPSORT ROUTINE */
void hpsort(unsigned long n, double ra[]);

/* Numerical Recipes standard error handler */
void nrerror(char error_text[]);

#endif // ROOTEXT_GENERALFUNCTIONS_H_
