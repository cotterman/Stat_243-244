#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include "timer.h"

void mydesk(double *x, int length_x, double *mean, double *var) {
  double length_x_dbl = length_x;
  double sum_x = 0;
  double sum_x2 = 0;
  int i = 0;
  for(i=0; i<length_x; i++) {
     sum_x2 += pow(x[i],2);
     sum_x += x[i];
  }
  *mean = sum_x/length_x;
  *var = (1/(length_x_dbl-1))*(sum_x2 - length_x_dbl * pow(*mean,2));
}

void provisional(double *x, int length_x, double *mean, double *var) {
  double length_x_dbl = length_x;
  double t_mean_old = 0;
  double t_mean_new = 0;
  double t_var = 0;
  int i = 0;
  for(i=0; i<length_x; i++){
     t_mean_new = t_mean_old + (1/(double)(i+1))*(x[i]-t_mean_old);
     t_var = t_var + (x[i]-t_mean_old)*(x[i]-t_mean_new);
     t_mean_old = t_mean_new;
  }
  *mean = t_mean_new;
  *var = t_var / (length_x_dbl-1);
}

void sas(double *x, int length_x, double *mean, double *var) {
  double length_x_dbl = length_x;
  double sum_x = 0;
  double sum_x2 = 0;
  int i = 0;
  double xi = 0;
  for(i=0; i<length_x; i++) {
	 xi = x[i] - x[0];
     sum_x2 += pow(xi,2);
     sum_x += xi;
  }
  *mean = sum_x/length_x + x[0];
  double t_mean = sum_x/length_x;
  *var = (1/(length_x_dbl-1))*(sum_x2 - length_x_dbl * pow(t_mean,2));
}

int is_d_equal(double d1, double d2, double diff_allowed){
   if (fabs(d1-d2) < diff_allowed) return 1;
   else return 0; /*return zero if d1 differs from d2 by more than diff_allow*/
}

void test_algorithms(double orig_mean, double orig_var, double *ystart, int length_y, int c_start, int c_end, int c_step, double diff_allowed, int errs_to_stop){
  int c = 0; /*all variables not passed to the function must be declared within the function*/
  int j = 0;
  double desk_mean;
  double desk_var;
  double prov_mean;
  double prov_var;
  double sas_mean;
  double sas_var;
  double correct_mean;
  double my_cv;
  double ynew[length_y];
  int error_count = 0;
  for(c=c_start; c<=c_end; c+=c_step) {
	  desk_mean = 0;
      desk_var = 0;
 	  prov_mean = 0;
      prov_var = 0;
  	  sas_mean = 0;
      sas_var = 0;
      error_count = 0;
      correct_mean = orig_mean+c;
      my_cv = sqrt(orig_var)/(orig_mean+c);
	  for(j=0; j<length_y; j++) {
		  ynew[j] = ystart[j]+c;
      }
      printf("Adding %d generates CV of %12.11lf (mean:%lf, var: %lf)\n", c, my_cv, correct_mean, orig_var);
      mydesk(ynew, length_y, &desk_mean, &desk_var);
        printf("Desk Calc: mean is %lf and variance is %lf\n", desk_mean, desk_var);
      provisional(ynew, length_y, &prov_mean, &prov_var);
        printf("Provisional Means Calc: mean is %lf and variance is %lf\n", prov_mean, prov_var);
      sas(ynew, length_y, &sas_mean, &sas_var);
        printf("Centering Calc: mean is %lf and variance is %lf\n", sas_mean, sas_var);

      //test_results(correct_mean, orig_var, desk_mean, desk_var, prov_mean, prov_var, sas_mean, sas_var);
      if (!is_d_equal(correct_mean, desk_mean, diff_allowed)) {
		  printf("Desk mean is wrong\n");
		  error_count++;
	  }
      if (!is_d_equal(orig_var, desk_var, diff_allowed)) {
		  printf("Desk variance is wrong\n");
		  error_count++;
	  }
      if (!is_d_equal(correct_mean, prov_mean, diff_allowed)) {
		  printf("Prov mean is wrong\n");
		  error_count++;
	  }
      if (!is_d_equal(orig_var, prov_var, diff_allowed)) {
		  printf("Prov variance is wrong\n");
		  error_count++;
	  }
      if (!is_d_equal(correct_mean, sas_mean, diff_allowed)) {
		  printf("SAS mean is wrong\n");
		  error_count++;
	  }
      if (!is_d_equal(orig_var, sas_var, diff_allowed)) {
		  printf("SAS variance is wrong\n");
		  error_count++;
	  }
      printf("Total errors: %d \n", error_count);
      if (error_count>=errs_to_stop) {
		  printf("\n");
		  break;
	  }
  }
}

uint64_t uniform_randoms (uint64_t seed, const int a, uint64_t modulus, int n, int N_to_gen, int printme, double *random_array, int suppress_out){
   uint64_t x_new = seed;
   double sas_mean = 0;
   double sas_var = 0;
   if (!suppress_out) printf("\n***Start of array*** \n");
   for(n=0; n<N_to_gen; n++) {
      x_new = (a*x_new) % modulus;
      random_array[n] = x_new/(double)modulus;
      if (printme) printf("%d: random_new is %15.14lf \n", n, random_array[n]);
    }
   if (!suppress_out) printf("***End of array*** \n");

   //calculate mean and variance for several different sequences;
     sas(random_array, N_to_gen, &sas_mean, &sas_var);
     if (!suppress_out) printf("Mean of my Uni(0,1) R.V. is %lf and variance is %lf\n", sas_mean, sas_var);

   return x_new;
}

void max_min(double *array_polar, int N_obs, double *mymax, double *mymin) {
	int myloop = 0;
	*mymax = array_polar[0];
	*mymin = array_polar[0];
    for(myloop=1;myloop<N_obs;myloop++){
		if (array_polar[myloop]>*mymax) *mymax=array_polar[myloop];
		else if (array_polar[myloop]<*mymin) *mymin=array_polar[myloop];
	}
}

void array_polar(uint64_t modulus, int a, uint64_t seed, int N_obs, double *array_polar_a, double *array_polar_b, int printvals, int N_samples, int presults){
	int myr = 0;
	for(myr=0;myr<N_samples;myr++){
	   double array_unif1[1] = { 0 };
	   double array_unif2[1] = { 0 };
	   double v1 = 0;
	   double v2 = 0;
	   double s  = 10;
	   double x1 = 0;
	   double x2 = 0;
       double intermediate = 0;
	   uint64_t x_new = 0;
       int n = 0;
       int myc = 0;
	   for(myc=0;myc<N_obs;myc++){
          while (s>=1){
             seed = uniform_randoms (seed, a, modulus, n, 1, 0, array_unif1, 1);
             seed = uniform_randoms (seed, a, modulus, n, 1, 0, array_unif2, 1);
             v1 = 2*array_unif1[0] - 1;
             v2 = 2*array_unif2[0] - 1;
             s = pow(v1,2) + pow(v2,2);
          }
          intermediate = sqrt( -2*log(s)/s );
          array_polar_a[myc] = v1*intermediate;
          array_polar_b[myc] = v2*intermediate;
          if (printvals) printf("x%d a is %lf, x%d b is %lf \n", myc, array_polar_a[myc], myc, array_polar_b[myc]);
          s = 10;
       }
       //find the means, variance, mins, and maxes for each sample;
          double sas_mean = 999;
          double sas_var = 999;
          double mymax = 999;
          double mymin = 999;
          sas(array_polar_a, N_obs, &sas_mean, &sas_var);
          max_min(array_polar_a, N_obs, &mymax, &mymin);
          if (presults) printf("M1, Sample %da: mean=%lf, var=%lf, max=%lf, min=%lf \n", myr, sas_mean, sas_var, mymax, mymin);
          sas_mean = 999;
          sas_var = 999;
          mymax = 999;
          mymin = 999;
          sas(array_polar_b, N_obs, &sas_mean, &sas_var);
          max_min(array_polar_b, N_obs, &mymax, &mymin);
          if (presults) printf("M1, Sample %db: mean=%lf, var=%lf, max=%lf, min=%lf \n", myr, sas_mean, sas_var, mymax, mymin);
   }
}

   void array_sumunif(uint64_t modulus, int a, uint64_t seed, int toadd, int N_obs, int N_samples, double *array_sum, int presults) {
       int iadd = 0;
       int n = 0;
       int iobs = 0;
       int isamp = 0;
       int N_for_seed = 100;
       double seed_array[N_for_seed];
       double array_current[N_obs];
       double myscale = sqrt( (3./toadd) );
       for(isamp=0;isamp<N_samples;isamp++){
          for(iadd=0;iadd<toadd;iadd++){
             uniform_randoms (seed, a, modulus, n, N_obs, 0, array_current, 1);
             for(iobs=0;iobs<N_obs;iobs++){
			    array_current[iobs] = array_current[iobs]*2 - 1; /*convert array_current to unif(-1,1) -- each of these has var of (1/3)*/
                if (iadd==0) array_sum[iobs] = (myscale) * array_current[iobs]; /*array_sum is the avg of our toadd array_currents, matching by element*/
                else array_sum[iobs] += (myscale) * array_current[iobs]; /*must scale in order for resulting R.N.s to have var of 1*/
	     	 }
             seed = uniform_randoms (seed, a, modulus, n, N_for_seed, 0, seed_array, 1);
	      }
          double sas_mean = 999;
          double sas_var = 999;
          double mymax = 999;
          double mymin = 999;
          sas(array_sum, N_obs, &sas_mean, &sas_var);
          max_min(array_sum, N_obs, &mymax, &mymin);
          if (presults) printf("M2, Sample %d: mean=%lf, var=%lf, max=%lf, min=%lf \n", isamp, sas_mean, sas_var, mymax, mymin);
       }
   }


/***********************************************************************/
/***********************************************************************/



int main(int argc, char *argv[]){

  //exercise control over which parts of program are run;
  int run_Q1 = 1;
  int run_Q2 = 1;
  int run_Q3 = 1;
  int run_Q4 = 1;


  /*************************************************************/
  /* Calculate the mean and variance using 3 different methods */
  /*************************************************************/
  printf("\n*************** Question 1 *****************\n");


  if (run_Q1) {

  //simple test data to use with algorithms
  double x[] = { 1, 2, 3, -4, -8, 25}; /*define an array called x*/
  int length_x = sizeof(x) / sizeof(x[0]); /*number of elements in x*/

  //method 1: Desk Calculator Algorithm
  double desk_mean = 0;
  double desk_var = 0;
  mydesk(x, length_x, &desk_mean, &desk_var); /* x refers to the adress of element zero*/
  printf("Desk Calc: mean is %lf and variance is %lf\n", desk_mean, desk_var);

  //method 2: Method of Provisional Means
  double prov_mean = 0;
  double prov_var = 0;
  provisional(x, length_x, &prov_mean, &prov_var);
  printf("Provisional Means Calc: mean is %lf and variance is %lf\n", prov_mean, prov_var);

  //method 3: Centering around the First Observation (used by SAS)
  double sas_mean = 0;
  double sas_var = 0;
  sas(x, length_x, &sas_mean, &sas_var);
  printf("Centering Calc: mean is %lf and variance is %lf\n", sas_mean, sas_var);

  }

  /************************/
  /* Ill-Conditioned Data */
  /************************/
  printf("\n\n*************** Question 2 *****************\n");


  if (run_Q2) {

  //create starting data, to be manipulated
  double ystart[] = {5, 3, 9, 6, 3, 7, 2, 4, 8, 2, 3, 5, 6};
  const int length_y = sizeof(ystart) / sizeof(ystart[0]);
  double orig_mean = 0;
  double orig_var = 0;
  double orig_cv = 0;
  double my_cv = 0;
  double correct_mean = 0;
  sas(ystart, length_y, &orig_mean, &orig_var);
  orig_cv = sqrt(orig_var)/orig_mean;
  printf("Baseline: mean is %lf, variance is %lf, and CV is %lf\n \n", orig_mean, orig_var, orig_cv);

  //loop over increasing constant to be added to sequence of numbers.  test each algorithm for mean and variance.
  int c_start;
  int c_end;
  int c_step;
  double diff_allowed;
  int errs_to_stop = 1;
     //detect when rounding errors create unequal answers;
        diff_allowed = 0;
        c_start=1;
        c_end=10000000;
        c_step=1;
        errs_to_stop = 1;
        test_algorithms (orig_mean, orig_var, ystart, length_y, c_start, c_end, c_step, diff_allowed, errs_to_stop);
     //detect when differences become visible;
        diff_allowed = .00001;
        c_start=261900;
    	c_end=10000000;
        c_step=50;
        errs_to_stop = 1;
        test_algorithms (orig_mean, orig_var, ystart, length_y, c_start, c_end, c_step, diff_allowed, errs_to_stop);
     //detect when algorithm(s) break down even more severeley;
        diff_allowed = 1;
        c_start=63000000;
	    c_end=1000000000;
        c_step=100000;
        errs_to_stop = 1;
        test_algorithms (orig_mean, orig_var, ystart, length_y, c_start, c_end, c_step, diff_allowed, errs_to_stop);
     //detect when algorithm(s) break down completely;
        diff_allowed = 5;
        c_start= 150000000;
	    c_end=  1000000000;
        c_step=   10000000;
        errs_to_stop = 1;
        test_algorithms (orig_mean, orig_var, ystart, length_y, c_start, c_end, c_step, diff_allowed, errs_to_stop);
     //detect when more than one algorithm breaks down;
        /*diff_allowed = .00001;
        c_start=63000000;
	    c_end=  100000000;
        c_step=100000;
        errs_to_stop = 2;
        test_algorithms (orig_mean, orig_var, ystart, length_y, c_start, c_end, c_step, diff_allowed, errs_to_stop);*/
        //the provisional means and SAS methods generated correct results even when a constant of 100,000,000 was added to data
  }


  /**************************/
  /* Uniform Random Numbers */
  /**************************/
  printf("\n*************** Question 3 *****************\n");
  if (run_Q3) {

  int n = 0;
  uint64_t modulus = 1; /*creates an unsigned integer that consumes 64 bits rather than 32 like normal*/
  modulus <<= 31; /*this bit-shift operator will essentially add 31 to the exponent of m, which is currently expressed as 2^0*/
  modulus*=2; /*we had to do this in 2 steps b/c addly enough, I couldn't get the bit shift to move 32 places*/
  const int a = 8003;
  uint64_t seed = 3498722111U;
  uint64_t x_new = 0;
  int N_to_gen = 0;

  N_to_gen = 10;
  double random_array1[N_to_gen];
  seed = uniform_randoms (seed, a, modulus, n, N_to_gen, 1, random_array1, 0);

  N_to_gen = 30;
  double random_array2[N_to_gen];
  seed = uniform_randoms (seed, a, modulus, n, N_to_gen, 1, random_array2, 0);

  N_to_gen = 1000;
  double random_array3[N_to_gen];
  seed = uniform_randoms (seed, a, modulus, n, N_to_gen, 0, random_array3, 0);

  N_to_gen = 50000;
  double random_array4[N_to_gen];
  seed = uniform_randoms (seed, a, modulus, n, N_to_gen, 0, random_array4, 0);

  }


  /***************************/
  /* Normal Random Variables */
  /***************************/
  printf("\n\n*************** Question 4 *****************\n");
  if (run_Q4) {

    uint64_t modulus = 1; /*creates an unsigned integer that consumes 64 bits rather than 32 like normal*/
	modulus <<= 31; /*this bit-shift operator will essentially add 31 to the exponent of m, which is currently expressed as 2^0*/
	modulus*=2; /*we had to do this in 2 steps b/c addly enough, I couldn't get the bit shift to move 32 places*/
	const int a = 8003;
	uint64_t seed = 3498725311U;

	//implement 2 methods for generating N(0,1) random numbers and examine results of each
       int N_samples = 10; /*number of samples to create -- make this an even number*/
       int N_obs = 20; /*number of obs in each sample*/
       printf("\nCreate and Analyze %d Samples of %d observations for each method (M1 and M2). \n", N_samples, N_obs);

       //generate N(0,1) random numbers using the Polar Method;
          double array_polar_1a[N_obs];
    	  double array_polar_1b[N_obs];
          array_polar(modulus, a, seed, N_obs, array_polar_1a, array_polar_1b, 0, N_samples/2, 1);

       //generate N(0,1) random numbers by adding together size Unif();
          int toadd = 6; /*number of vectors of Unif(-1,1) to combine to create N(0,1)*/
          double array_sum[N_obs];
          array_sumunif(modulus, a, seed, toadd, N_obs, N_samples, array_sum, 1);

    //generate N(0,1) randon numbers again, but now time how long it takes to complete this task many times under our 2 methods
       N_samples = 10000; /*number of samples to create -- make this an even number*/
       N_obs = 10000; /*number of obs in each sample*/
       printf("\nTime the creation of %d samples of %d observations each for methods M1 and M2. \n", N_samples, N_obs);

       //generate N(0,1) random numbers using the Polar Method;
          TIME0
            double array_polar_biga[N_obs];
    	    double array_polar_bigb[N_obs];
            array_polar(modulus, a, seed, N_obs, array_polar_biga, array_polar_bigb, 0, N_samples/2, 0);
 		  TIME1("Time for Polar Method")

       //generate N(0,1) random numbers by adding together size Unif();
          TIME0
            toadd = 6; /*number of vectors of Unif(-1,1) to combine to create N(0,1)*/
            double array_sum_big[N_obs];
            array_sumunif(modulus, a, seed, toadd, N_obs, N_samples, array_sum_big, 0);
          TIME1("Time for Summing-of-Uniforms Method")

  } /*end of Q4*/

  return 0;
}

