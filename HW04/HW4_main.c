#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include "HW4_header.h" //this basically means to copy and paste code from HW4_header.h into the top of this code document
FILE *infile;

//This program may use functions located in HW4_functions.c.  Here are the corresponding function protocols:
   int is_d_equal (double, double, double);
   void var_mean (double *, int, double *, double *);
   void max_min (double *, int, double *, double *);
   uint64_t uniform_randoms (uint64_t, const int, uint64_t, int, double *, int, int);
   uint64_t normal_randoms (uint64_t, int, uint64_t, int, double *, int, int);
   void GS_ortho(double *, double *, double *, int, int, int);
   void cholesky(int, double *, double *, double *, int );
   void GS_regress(double *, double *, double *, int, int, int);
   void InvUTMat (double *, double *, int);
//Will also need to access the beatonsweep.c file when compiling
	void beatonsweep(double *, int , int *, int);

//Note: it seems to make more sense to instead include functions that I want to access in a "header file"
	//this is what I did for the remaining functions that I may want to use in this code -- they are in HW4_header.h
	//the exception is that I am including here the new functions that I wrote for the non-linear regression


//Calculate the Jacobian Matrix, J, using the functional form of the prob1 data
	//inputs: parameter test values, data values, number of observations
void Jacobian1(double *X, double *Beta, int n, double *jacobian){
	int i;
	for (i=0; i<n; i++){
		//1st column of jacobian matrix will contain values for the partial derivative of F with respect to beta[0]
			jacobian[i] = - X[i] / (Beta[1] + X[i]);
		//2nd column of jacobian matrix will contain values for the partial derivative of F with respect to beta[1]
			jacobian[i+n] =  (Beta[0] * X[i]) / ((Beta[1] + X[i]) * (Beta[1] + X[i]));
	}
	printf("X matrix \n");
	printmat(X,n,1);
	printf("Beta matrix \n");
	printmat(Beta,2,1);
	printf("Jacobian matrix \n");
	printmat(jacobian,n,2);
}

//Calculate the Jacobian Matrix, J, using the functional form of the prob2 data
	//inputs: parameter test values, data values, number of observations
void Jacobian2(double *X, double *Beta, int n, double *jacobian){
	int i;
	//note that 1st row of X is u, 2nd row is v, and 3rd row is w
	for (i=0; i<n; i++){
		//1st column of jacobian matrix will contain values for the partial derivative of F with respect to beta[0]
			jacobian[i] = 1;
		//2nd column of jacobian matrix will contain values for the partial derivative of F with respect to beta[1]
			jacobian[i+n] = (-X[i] * X[i+n])/ (pow( (Beta[1]*X[i+n] + Beta[2]*X[i+2*n]),2 ));
		//3nd column of jacobian matrix will contain values for the partial derivative of F with respect to beta[2]
			jacobian[i+2*n] = (-X[i] * X[i+2*n])/(pow( (Beta[1]*X[i+n] + Beta[2]*X[i+2*n]),2));
	}

}

//Calculate the residuals using the functional form of prob1 data
	//inputs: parameter test values, data values, number of observations
double MyFunc1(double *X, double *Y, double *Beta, int n, double *resids){
	int i;
	double Y_hat[n];
	double RSS;
	zeromat(Y_hat,n);
	for (i=0; i<n; i++){
		Y_hat[i] = (Beta[0]*X[i])/(Beta[1]+X[i]); //2nd column of Beta matrix contains estimates for k
		resids[i] = Y_hat[i]-Y[i];
	}
	printf("Residuals \n");
	printmat(resids, n, 1);
	RSS = dots(resids, resids, 1, 1, n); //note that dot product of the residual vector with itself is the RSS
	printf("RSS is %lf \n", RSS);
	return RSS;
}

//Calculate the residuals using the functional form of prob2 data
	//inputs: parameter test values, data values, number of observations
double MyFunc2(double *X, double *Y, double *Beta, int n, double *resids){
	int i;
	double Y_hat[n];
	double RSS;
	zeromat(Y_hat,n);
	for (i=0; i<n; i++){
		Y_hat[i] = Beta[0] +  X[i] / (Beta[1]*X[i+n] + Beta[2]*X[i+2*n]);
		resids[i] = Y_hat[i]-Y[i];
	}
	printf("Residuals \n");
	printmat(resids, n, 1);
	RSS = dots(resids, resids, 1, 1, n);
	printf("RSS is %lf \n", RSS);
	return RSS;
}

//Apply the Gauss-Newton algorithm to solve non-linear regression problems
void NonLinReg(double *X, double *Y, int n, int p, double *Beta, void (*gradfunc)(double*, double*, int, double*), double (*func)(double*, double*, double*, int, double*)) {
	int i, j;
	double RSS = 0;
	double RSS_new = 0;
	int done = 0;
	double stepfrac = 1;
	double eps1 = 0.00000001; //will use in stopping rule
	double eps2 = 10*eps1;

	while (done==0) {

	//use parameter starting values to find Jacobian matrix, J
		double jacobian[n*p];
		zeromat (jacobian, n*p);
		(*gradfunc)(X, Beta, n, jacobian);

	//calculate the residuals (and RSS)
		double resids[n];
		zeromat (resids, n);
		RSS = (*func)(X, Y, Beta, n, resids);

	//calculate step = (J'J)^{-1}J'r using GS linear regression program from assignment 3
		//Use J (nxp) where X normally goes, and residual vector where Y vector normally goes
		double stepsize[p];
		zeromat (stepsize, p);
		GS_regress(stepsize, resids, jacobian, n, p, 1);
		printf("steps to take \n");
		printmat(stepsize,p,1);

	//calculate new parameter test values using Gauss-Newton equation
		double Beta_new[p];
		for(j=0;j<p;j++){
			Beta_new[j] = Beta[j]-stepsize[j];
		}
		printf("updated parameter values \n");
		printmat(Beta_new,p,1);

	//calculate RSS under new parameter test values
		double resids_new[n];
		zeromat (resids_new, n);
		RSS_new = (*func)(X, Y, Beta_new, n, resids_new);

	//compare RSS under new parameter values vs. previous parameter values
		//if divergence, do step-halving until no divergence, or until we have cut step in half 5 times
		int nhalfs = 0;
		while ((RSS_new>RSS) && nhalfs<=5){
			for(j=0; j<p; j++)  Beta_new[j] = Beta[j] - stepfrac * stepsize[j];
			zeromat (resids_new, n);
			RSS_new = (*func)(X, Y, Beta_new, n, resids_new);
			stepfrac /= 2;
			done = 1;
			nhalfs++;
		}
		//if no divergence, see if stopping rule says to stop
		for(j=0; j<p; j++) {
			if(done==0  &&  abs(Beta[j] - Beta_new[j]) > eps1 * (abs(Beta[j]) + eps2) ) done = 0;
			else done = 1;
			Beta[j] = Beta_new[j];
		}

	} /*while loop*/

}

/******************************************************************************/
/******************************************************************************/


int main(int argc, char *argv[]){

int i, j, k;

	//exercise control over which parts of program are run;
	  int run_Q1 = 1;
	  int run_Q2 = 1;

if (run_Q1) {
printf("\n**************** Sweep *****************\n\n");

	int p=3; //dimension of square matrix
	int nsweep=3; //number of columns that we want to sweep
	double matA[] = { 11, 15, 2, 8 , 30, 4, 5, 7, 13}; //matrix to sweep
	int indices[] = { 1, 2, 3 }; //indices of columns we want to seep

	printf("Before the sweep:\n");
	printmat(matA, p, p);

	beatonsweep(matA, p, indices, nsweep);
	printf("\nAfter the sweep:\n");
	printmat(matA, p, p);

	beatonsweep(matA, p, indices, nsweep);
	printf("\nAfter 2 sweeps (should look like before sweep):\n");
	printmat(matA, p, p);

} /*end of run_Q1*/



if (run_Q2) {
printf("\n********* Non-Linear Regression ********\n\n");


//note that static variables declared here are stored as globals in the sense of being stored in memory on the global data segment
	// but static variables declared here will NOT be global in scope.
	// must pass parameters as arguments to functions

//I will be specify which gradient function and which funcform function to call in my parameter list given to the NonLinReg function
	void (*gradient)(double *, double *, int, double *);
	double (*funcform)(double *, double *, double *, int, double *);


//non-linear regression for prob1.dat
	//inputs obtained from background and data files
	int n1 = 8; //number of observations
	int m1 = 1; //number of independent variables
	int p1 = 2; //number of parameters
	double first_dbl, second_dbl;
	infile=fopen("prob1.dat", "r"); //1st col is X, 2nd col is Y
		if (!infile) {
			printf("Cannot Open File");
			return 1;
		}
        double matX1[n1*m1];
        double matY1[n1*1];
		for (i=0; i<n1; i++) {
			fscanf(infile, "%lf  %lf", &first_dbl, &second_dbl); //fscanf reads row-by-row
			matX1[i] = first_dbl;
			matY1[i] = second_dbl;
	    }
		fclose(infile);
		printf("Matrix X1 (from file), non-linear regression original\n");
		printmat(matX1,n1,m1);
		printf("Matrix Y1 (from file), non-linear regression original\n");
		printmat(matY1,n1,1);
	//starting values for parameter estimates, based on result from grid search in R
		double beta1[] =  { 10, 5 };
	//run the non-linear regression(!)
		//NonLinReg(matX1, matY1, n1, p1, beta1, Jacobian1, MyFunc1);


//inputs for prob2.dat (1st col is Y, 2-4 cols are X)
	int n2 = 15; //number of observations
	int m2 = 3;  //number of independent variables
	int p2 = 3;  //number of parameters
	double dbl_1, dbl_2, dbl_3, dbl_4;
	infile=fopen("prob2.dat", "r");
		if (!infile) {
			printf("Cannot Open File");
			return 1;
		}
        double matX2[n2*m2];
        double matY2[n2*1];
		for (i=0; i<n2; i++) {
			fscanf(infile, "%lf %lf %lf %lf", &dbl_1, &dbl_2, &dbl_3, &dbl_4);
			matY2[i] = dbl_1;
			matX2[i] = dbl_2;
			matX2[i+1*n2] = dbl_3;
			matX2[i+2*n2] = dbl_4;
	    }
		fclose(infile);
		printf("Matrix X2 (from file), non-linear regression original\n");
		printmat(matX2,n2,m2);
		printf("Matrix Y2 (from file), non-linear regression original\n");
		printmat(matY2,n2,1);
	//starting values for parameter estimates, based on result from grid search in R
		double beta2[] =  { 0, 0, 10 };
	//run the non-linear regression(!)
		//NonLinReg(matX2, matY2, n2, p2, beta2, Jacobian2, MyFunc2);


} /*end run_Q2*/


}
