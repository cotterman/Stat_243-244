#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define LOC(a,b)  ((a) * p + (b))
#include "Stat244_hw3_header1.h" //this basically means to copy and paste code from the .h file into the top of this code document
#include "timer.h"
FILE *infile;

/* This program calls many functions located in other .c and .f files.  To compile and run this code, see the file called make for instructions. */


/******************************************************************************/
/******************************************************************************/


int main(int argc, char *argv[]){

	//exercise control over which parts of program are run;
	  int run_Q1 = 1;
	  int run_Q2 = 1;
	  int run_downhill = 1; //set to one if you wish to use the downhill simplex minimization method for ARMA estimation
	  int run_powell = 1; //set to one if you wish to use Powell's minimization method for ARMA estimation
	  int timeme, myi;
	  int mytimes=1000; //number times to run each ARMA estimation procedure (for time estimation purposes)
	  int mytimescov = 1000; //number times to run the covariance procedures (for time estimation purposes)
	  int printme=0; //display output (choose 0 if only execution times are desired)

	int i, k, t, j, r, c, p, q, maxk, goodparams, mydata, nparams, nobs;
	double myd;
	char *filename = NULL;
	nobs = 1000; ////number of observations (same for all datasets used here)
	double sdata[nobs]; //temp place to store each column

  //specify which dataset you are using (numbered 1 thru 5), or loop thru them all
  for(mydata=1;mydata<=5;mydata++){

	//specify parameters that are specific to the data
		if (mydata==1) filename = "sdata1";
		if (mydata==2) filename = "sdata2";
		if (mydata==3) filename = "sdata3";
		if (mydata==4) filename = "sdata4";
		if (mydata==5) filename = "sdata5";
		if (mydata==1 || mydata==4 || mydata==5) {
			p = 1;
			q = 1;
		}
		else if (mydata==2) {
			p = 3;
			q = 1;
		}
		else if (mydata==3) {
			p = 2;
			q = 2;
		}
		nparams = p + q;
		double* alpha = (double*)malloc(p*sizeof(double));
		double* beta = (double*)malloc(q*sizeof(double));
		double* vertices = (double*)malloc((nparams+1)*nparams*sizeof(double)); //create starting vertices for downhill method
		double* copyvertices = (double*)malloc((nparams+1)*nparams*sizeof(double));
		double* startpt = (double*)malloc(nparams*sizeof(double)); //create vector of starting values for Powell's method
		double* copystartpt = (double*)malloc(nparams*sizeof(double));
		double* directions = (double*)malloc(nparams*nparams*sizeof(double)); //enter numbers to be read by row into "double-dimension matrix" for Powell's method
		if (mydata==1 || mydata==5) {
			alpha[0] = .3; //actual parameters used to create data (mydata=5 was seeded with a larger error variance)
			beta[0] = .3;
		}
		else if (mydata==2) {
			alpha[0] = .3; //actual parameters used to create data
			alpha[1] = .2;
			alpha[2] = .1;
			beta[0] = .3;
		}
		else if (mydata==3) {
			alpha[0] = .3; //actual parameters used to create data
			alpha[1] = .2;
			beta[0] = .3;
			beta[1] = .2;
		}
		else if (mydata==4) {
			alpha[0] = .5; //actual parameters used to create data
			beta[0] = .5;
		}
		if (mydata==1 || mydata==4 || mydata==5){
			vertices[0] = .4; //should be length (nparams+1)*nparams
			vertices[1] = .4;
			vertices[2] = .2;
			vertices[3] = .2;
			vertices[4] = .1;
			vertices[5] = .1;
			startpt[0] = .4; //should be length nparams
			startpt[1] = .4;
			directions[0] = .1; //should be length nparams * nparams
			directions[1] = 0;
			directions[2] = 0;
			directions[3] = .1;
		}
		else if (mydata==2 | mydata==3){
			double tmp[] = {.3,.3,.3,.3,.4,.3,.2,.1,.2,.2,.2,.2,.3,.1,.1,.1,.2,.1,.2,.3}; //should be length (nparams+1)*nparams
			memcpy(vertices, tmp, sizeof(tmp));
			startpt[0] = .4; //should be length nparams
			startpt[1] = .3;
			startpt[2] = .2;
			startpt[3] = .1;
			memset(directions, 0, sizeof(double)*nparams*nparams); //sets everything to zero
			for(i=0;i<nparams;i++)
			    directions[i*nparams+i] = .1;
		}
	//initialize even more variables;
		double t_alpha[p];
		double t_beta[q];
		double params[p + q]; //should contain list p AR coefs followed by q MA coefs
		for(i=0;i<p+q;i++){
			if (i<p) params[i] = alpha[i];
			else params[i] = beta[i-p];
		}

	printf("\n********************************************************\n");
	printf("\n********************** USING %s *******************\n", filename);
	printf("\n********************************************************\n \n");
		goodparams = 1 - chkparms(params,p,q); //chkparms returns 1 if parameters are bad; 0 if they are ok
			printf("\n Verification parameters alpha=");
			for(i=0;i<p;i++){
				printf("%lf, ",alpha[i]);
			}
			printf("beta=");
			for(i=0;i<q;i++){
				printf("%lf, ",beta[i]);
			}
			printf("are good = %d. \n \n", goodparams);

	//Read in the time-series data generated in R
		infile=fopen(filename, "r");
			if (!infile) {
				printf("Cannot Open File");
				return 1;
			}
			for (i=0; i<nobs; i++) {
					fscanf(infile, "%lf", &myd); //fscanf reads row-by-row
					sdata[i] = myd;
		    }
			fclose(infile);
		if (printme==1) printf("Original Data (from R), first 25 observations \n");
		for (i=0; i<25; i++) {
			if (printme==1) printf("%lf  ",sdata[i]);
		}

if (run_Q1) {
if (printme==1) printf("\n\n*************** Q1: ARMA Parameter Estimation ****************\n\n");

		double lik = 0;
		double likstar[(nparams+1)];
		zeromat (likstar, (nparams+1));
		double ftol = .00000000000001;


	//Estimate parameters of the ARMA process using the Downhill Simplex Method to maximize likelihood (Nelder and Mead)
	if (run_downhill==1){
	TIME0
	//create copy of startpt so we can run this multiple times (must re-inialize to starting values after each run);
		for(myi=0;myi<(nparams+1)*nparams;myi++){
			copyvertices[myi] = vertices[myi];
		}
	for (timeme=0; timeme<mytimes; timeme++){

		//My modified amoeba function
			//void arma_amoeba(double sdata[], int nobs, double **p, double y[], int ar_p, int ma_q, double ftol, int *nfunk)
			//p is a (ndim+1) x (ndim) input matrix, containing initial guesses.  Each row is one vertex of the starting simplex.
			//y is an ndim+1 vector containing the evaluations at each vertex of p
			//ftol is the fractional convergence tolerence (can  be set to something close to the machine epsilon)
			//nfunk returns the number of function evaluations taken
			//OUTPUT: new p matrix containing optimal parameter values and y vector containing corresponding likstar
		//create matrix of starting vertices
			  if (printme==1) printf("Number of parameters: %d \n \n", nparams);
			  double **vptarray; //vptarray will be an array of (nparam+1) pointers;
			  double *verticespt = vertices; //vertices will be a matrix.  each row of which will be pointed to by the pointers contained in vptarray.
			  /* By allocating all the memory for vptarray in a single malloc call and then setting the row pointers, a dereferenced version
		    	of the matrix (either *vptarray or vptarray[0]) can be passed to routines which expect the matrix to be stored as a vector.
		    	The first allocation only provides memory for the row pointers, *not* for the actual data values which will be stored in vptarray. */
		        vptarray = (double**)malloc((unsigned)((nparams+1) * sizeof(double*)));
 		 	  //here we tell vptarray[i] to start at each row i of vertices
		 	    for(i=0;i<(nparams+1);i++,verticespt += nparams)
		        	vptarray[i] = verticespt; //assumes vertices is stored  by rows
		 	  if (printme==1) printf("Original vertices matrix:\n");
		    	for(i=0;i<(nparams+1);i++) { //row loop
		     		if (printme==1) printf("Row %d: ", i);
		    		for(j=0;j<nparams;j++) { //col loop
		        	    if (printme==1) printf("%lf  ", vptarray[i][j]); //note that *vptarray is the same as vptarray[0] -- points to the 1st element of vertices
			    	}
			    	if (printme==1) printf("\n");
		    	}
		    	if (printme==1) printf("\n");
		  //find the likstars that correspond to each of the "rows" of "matrix vptarray" and place them in vector y
		  		for(r=0;r<(nparams+1);r++){ //iterate through "rows" of matrix "vptarray"
					for(c=0;c<nparams;c++){
						if (c<p) t_alpha[c] = vptarray[r][c];
						else t_beta[c-p] = vptarray[r][c];
					}
					/*printf("Starting value for alpha:");
						for(j=0;j<p;j++) {
							printf("%lf ", t_alpha[j]);
						}
						printf("\n");
					printf("Starting value for beta:");
						for(j=0;j<q;j++) {
							printf("%lf ", t_beta[j]);
						}
						printf("\n");*/
		  			doarma(sdata, t_alpha, t_beta, &nobs, &p, &q, &lik);
		  			likstar[r] = lik;
		  			if (printme==1) printf("Starting likstar: %lf \n", lik);
		  			if (printme==1) printf("\n");
				}

		  int nfunk = 0;

		  if (printme==1) printf("AMOEBA MINIMIZATION PROCEDURE BEGINS \n");
		  arma_amoeba(sdata, nobs, vptarray, likstar, p, q, ftol, &nfunk);
		  if (printme==1) printf("\n");

		  //Print our ending vertices matrix and corresponding L* values
		  		if (printme==1) printf("ARMA parameter values:\n");
		  		for(i=0;i<(nparams+1);i++) { //row loop
		  		    if (printme==1) printf("Row %d: ", i);
		  		    for(j=0;j<nparams;j++) { //col loop
		  		    	if (printme==1) printf("%lf  ", vptarray[i][j]); //note that *vptarray is the same as vptarray[0] -- points to the 1st element of vertices
		  			}
			    	if (printme==1) printf("\n");
				}
			    if (printme==1) printf("Minimized likstar values:\n");
			    for(i=0;i<(nparams+1);i++){
				    if (printme==1) printf("%lf ",likstar[i]);
				}
				if (printme==1) printf("\n");
				if (printme==1) printf("Number of function evaluations: %d \n \n",nfunk);
		//reinitialize starting values
			for(myi=0;myi<(nparams+1)*nparams;myi++){
				vertices[myi] = copyvertices[myi];
			}
	} //end time loop
	printf("Time to run the Downhill Simplex Method %i times", timeme);
	TIME1(" ")
	printf("\n");
	} //end of downhill simplex method

	//Estimate parameters of the ARMA process using Powell's Method to maximize likelihood
		//powell(p,xi,n,ftol,iter,fret,func)
		//arma_powell(double sdata[], double p[], int nobs, double **xi, int ar_p, int ma_q, double ftol, int *iter, double *fret)
			//p is the initial starting point ([1,...n])
			//n is the number of parameters
			//xi is the n x n pointer-array matrix whose columns are the set of directions n_i (usually unit vectors)
			//ftol is the fractional tolerance in the function value such that failure to decrease by more than this amount on one iteration signals doneness.
			//OUTPUT: p is set to the best point found
				//fret is the returned function value at p
				//iter is the number of iterations taken

	if (run_powell==1){
	TIME0;
	//create copy of startpt so we can run this multiple times (must re-inialize to starting values after each run);
		for(myi=0;myi<nparams;myi++){
			copystartpt[myi] = startpt[myi];
		}
	for (timeme=0; timeme<mytimes; timeme++){
		int niter = 0;
		double plik = 0;
		//create array of pointers for set of directions
			  double **dptarray; //vptarray will be an array of (nparam) pointers;
			  double *directpt = directions;
		        dptarray = (double**)malloc((unsigned)((nparams) * sizeof(double*)));
 		 	  //here we tell dptarray[i] to start at each row i of directions
		 	    for(i=0;i<(nparams);i++,directpt += nparams)
		        	dptarray[i] = directpt; //assumes directions is stored  by rows

		if (printme==1) printf("POWELL'S MINIMIZATION PROCEDURE BEGINS \n");
		arma_powell(sdata, startpt, nobs, dptarray, p, q, ftol, &niter, &plik);
		if (printme==1) printf("\n");
		//Print estimated parameters and corresponding L* value
		  		if (printme==1) printf("ARMA parameter values (Powell):\n");
		  		for(i=0;i<(nparams);i++) { //row loop
		  		    if (printme==1) printf("%lf  ", startpt[i]);
		  		}
			    if (printme==1) printf("\n");
			    if (printme==1) printf("Minimized likstar value: %lf \n", plik);
				if (printme==1) printf("Number of function evaluations: %d \n \n",niter);
		//reinitialize starting values
			for(myi=0;myi<nparams;myi++){
				startpt[myi] = copystartpt[myi];
			}
	} //end time loop
	printf("Time to run Powell's Method %i times", timeme);
	TIME1(" ")
	printf("\n");
	} //end of powell method

  } //end looping over datasets

} /*end run_Q1*/

if (run_Q2) {
if (printme==1) printf("\n\n*************** Q2: Autocovariance Estimation ****************\n \n");

	int maxk=20; //number of autocovariances we wish to calculate.
	double covar[maxk]; //vector in which we will place the autocovariances
	double t_sum; //temp value holder
	double mean = 0; //mean value of data
	double var = 0;
	int newobs = nobs+maxk;
	double t_sdat[nobs+maxk];
	double m_sdat[nobs+maxk];
	double covarfft[newobs];

	var_mean(sdata, nobs, &mean, &var);
	if (printme==1) printf("Mean is %lf and variance is %lf \n \n",mean, var);

 TIME0;
	for (timeme=0; timeme<mytimescov; timeme++){

	//Calculate autocovariances using the basic formula
	zeromat(t_sdat,nobs+maxk);
	for(k=0; k<maxk; k++){
		t_sum = 0;
		for(t=0; t<nobs-k; t++){
			t_sum += (sdata[t]-mean)*(sdata[t+k]-mean);
		}
		covar[k] = t_sum/(nobs);
	}
	if (printme==1) printf("The first %d covariances (using basic formula) are: \n",maxk);
	for (i=0; i<maxk; i++) {
		if (printme==1) printf("%lf  ",covar[i]);
	}
 } //end time loop
 printf("\n");
 printf("Time to get covariances with basic method %i times", timeme);
 TIME1(" ")
 printf("\n \n");

 TIME0;
 for (timeme=0; timeme<mytimescov; timeme++){

	//Calculate autocovariances using the Fast Fourier Transform
	   //Substract the mean from each observation.  Leave K zeros the end of the vector
	   		for (i=0; i<nobs; i++) {
				t_sdat[i] = sdata[i] - mean;
			}
			//printf("After subtracting mean:\n");
			//printmat(t_sdat,1,maxk);
	   //Calculate the Fast Fourier Transform -- use drffti.f (for forward) and
	   		double waf[2*newobs+15];
			drffti_(&newobs, waf); //initialization
			drfftf_(&newobs, t_sdat, waf); //transform -- output will be in t_sdat, with a different element for the real and complex components.
			//printf("After FFT:\n");
			//printmat(t_sdat,1,maxk);
	   //Calculate the squared modulus.  (This means summing the real part squared with the complex part squared.)
			for(i=0;i<newobs;i++){
				if (i==0) m_sdat[i] = 0;
				else if (2*i < newobs) {
					m_sdat[i] = (t_sdat[2*i] * t_sdat[2*i]) + (t_sdat[2*i-1] * t_sdat[2*i-1]);
					m_sdat[newobs - i] = (t_sdat[2*i] * t_sdat[2*i]) + (t_sdat[2*i-1] * t_sdat[2*i-1]);
				}
				else if (2*i == newobs) m_sdat[i] = (t_sdat[2*i-1] * t_sdat[2*i-1]);
				//m_sdat[i] = (t_sdat[i] * t_sdat[i]);
			}
			//printf("After squaring each term:\n");
			//printmat(m_sdat,1,maxk); //This matrix matches with the one I generated in R
	   //Perform an inverse Fast Fourier Transform and divide each element by (N+k)*N
	   		double wab[2*newobs+15];
	   		drffti_(&newobs, wab); //initialization
			drfftb_(&newobs, m_sdat, wab); //transform -- output will be in m_sdat
			for(i=0;i<newobs;i++){
				covarfft[i] = m_sdat[i]/(double)(newobs*nobs);
			}
		//the kth element of the resulting vector will be the kth covariance value
			if (printme==1) printf("Covariances using FFT:\n");
			if (printme==1) printmat(covarfft,1,maxk);

 } //end time loop
 printf("\n");
 printf("Time to get covariances using the FFT %i times", timeme);
 TIME1(" ")
 printf("\n \n");

} /*end run_Q2*/


}
