#ifndef ROOTEXT_RANDOMFUNCTIONS_H_
#define ROOTEXT_RANDOMFUNCTIONS_H_

#include <vector>





/*** CALL FUNCTIONS ***/

/* i.e. the functions you actually will call */

double random_uniform(double min, double max);

int random_poisson(double mean);

double random_gaussian(double mean, double sigma);


/*** MISC FUNCTIONS ***/

/*  INITIALIZATION SUBROUTINE FOR ran1()  */
/* calls ran1 with external static my_idum negative to intialize. 
   then do not change my idum */
void initialize_ran1();

void initialize_ran1(long ranSeed);

long get_ran1_seed();


/*** WORK FUNCTIONS ***/

/* i.e. the functions that actually do the work */

/*   RANDOM NUMBER GENERATOR ran1()  */
/* recommended random number generator from Numerical Recipes in C */
/* initialize with initialize_ran1() */
float ran1(long *idum);

/* POISSON RANDOM DEVIATE */
/* From Numerical Recipes In C */
double poidev(double xm, long *idum);

/* GAUSSIAN NORMAL DEVIATE */
/* From Numerical Recipes In C */
double gasdev(long *idum);


/*** C++ FUNCTIONS ***/

// For reasons not understood, this full function definition has to be 
// *here*, in the header, in order for other compilable functions to
// find it.  (Do not try moving the definition to the .C file.)

template <class T>
void MergeVectorsRandomly (const vector<T>& v1, const vector<T>& v2,
			   vector<T>& m)
{
  m.clear();
  int n1 = v1.size();
  int n2 = v2.size();
  while ( n1>0 || n2>0 ) {
    // odds of next item being from vector 1 or 2 is weighted by
    // the remaining number of elements in these vectors
    double prob = random_uniform (0., n1+n2);
    if (prob<n1) {
      --n1;
      m.push_back(v1[n1]);
    } else {
      --n2;
      m.push_back(v2[n2]);
    }
  }
}


#endif // ROOTEXT_RANDOMFUNCTIONS_H_


