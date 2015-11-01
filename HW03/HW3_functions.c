#include <stdint.h>
#include <stdio.h>
#include <math.h>

//test to see if 2 doubles are within diff_allowed units of one another
int is_d_equal(double d1, double d2, double diff_allowed){
   if (fabs(d1-d2) < diff_allowed) return 1;
   else return 0; /*return zero if d1 differs from d2 by more than diff_allow*/
}

//obtain the mean and variance of a vector using the SAS algorithm.
void var_mean(double *x, int length_x, double *mean, double *var) {
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

//obtain the min and max values contained in a vector
void max_min(double *array, int N_obs, double *mymax, double *mymin) {
	int myloop = 0;
	*mymax = array[0];
	*mymin = array[0];
    for(myloop=1;myloop<N_obs;myloop++){
		if (array[myloop]>*mymax) *mymax=array[myloop];
		else if (array[myloop]<*mymin) *mymin=array[myloop];
	}
}

//generate array of length N_obs containing U(0,1) random numbers;
uint64_t uniform_randoms(uint64_t seed, const int a, uint64_t modulus, int N_obs, double *random_array, int printme, int givestats){
   uint64_t x_new = seed;
   double sas_mean = 0;
   double sas_var = 0;
   int n = 0;
   if (printme) printf("\n***Start of array*** \n");
   for(n=0; n<N_obs; n++) {
      x_new = (a*x_new) % modulus;
      random_array[n] = x_new/(double)modulus;
      if (printme) printf("element %d: %15.14lf \n", n, random_array[n]);
    }
   if (printme) printf("***End of array*** \n");
   if (givestats) {
          double sas_mean = 999;
          double sas_var = 999;
          double mymax = 999;
          double mymin = 999;
          var_mean(random_array, N_obs, &sas_mean, &sas_var);
          max_min(random_array, N_obs, &mymax, &mymin);
          printf("Summary: mean=%lf, var=%lf, max=%lf, min=%lf \n \n", sas_mean, sas_var, mymax, mymin);
   }
   return x_new;
}

// dmalloc
double *dmalloc(unsigned long n) {
	double *x;
	x=(double*)malloc((size_t)n*sizeof(double));
	if(x==NULL){
		printf("Could not allocate %ld doubles \n",n);
		exit(1);
	}
	return(x);
}

//create array of N_obs each that contain N(0,1) random numbers (array_sumunif)
uint64_t normal_randoms (uint64_t modulus, int a, uint64_t seed, int N_obs, double *random_array, int printme, int givestats) {
    int toadd = 6; /*number of vectors of Unif(-1,1) to combine to create N(0,1)*/
    int iadd = 0;
    int n = 0;
    int iobs = 0;
    int N_for_seed = 100;
    double seed_array[N_for_seed];
    double array_current[N_obs];
    double myscale = sqrt( (3./toadd) );
    for(iadd=0;iadd<toadd;iadd++){
       uniform_randoms (seed, a, modulus, N_obs, array_current, 0, 0);
       for(iobs=0;iobs<N_obs;iobs++){
	      array_current[iobs] = array_current[iobs]*2 - 1; /*convert array_current to unif(-1,1) -- each of these has var of (1/3)*/
          if (iadd==0) random_array[iobs] = (myscale) * array_current[iobs]; /*random_array is the avg of our toadd array_currents, matching by element*/
          else random_array[iobs] += (myscale) * array_current[iobs]; /*must scale in order for resulting R.N.s to have var of 1*/
       }
       seed = uniform_randoms (seed, a, modulus, N_for_seed, seed_array, 0, 0);
    }
    if (printme) {
		printf("\n***Start of array*** \n");
		for(n=0; n<N_obs; n++) {
		     printf("element %d: %15.14lf \n", n, random_array[n]);
		}
        printf("***End of array*** \n");
	}
	if (givestats) {
       double sas_mean = 999;
       double sas_var = 999;
       double mymax = 999;
       double mymin = 999;
       var_mean(random_array, N_obs, &sas_mean, &sas_var);
       max_min(random_array, N_obs, &mymax, &mymin);
       printf("Summary: mean=%lf, var=%lf, max=%lf, min=%lf \n \n", sas_mean, sas_var, mymax, mymin);
    }
    return seed;
}


