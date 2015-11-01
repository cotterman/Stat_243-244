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

void cholesky(int dimA, double *matA, double *matL, double *matU, int printme) {
  double sumin = 0;
  double diag[dimA];
  int i, j, k;
  for(i=0; i<dimA; i++) {
	  for(j=i; j<dimA; j++) {
		  for(sumin=matA[j*dimA+i], k=i-1; k>=0; k--) sumin -= matL[k*dimA+i] * matL[k*dimA+j];
		  if (i==j) {
			  if (sumin<0.0) printf("matA is not positive definite = fail.\n");
			  diag[i]=sqrt(sumin);
			  matL[i*dimA+j]=diag[i];
			  matU[j*dimA+i]=diag[i];
		  }
		  else {
			  matL[i*dimA+j]=sumin/diag[i];
		      matU[j*dimA+i]=sumin/diag[i];
	      }
	  }
  }
  if (printme) {
	  printf("The original matrix A:\n");
	  printmat(matA, dimA, dimA);
	  printf("The lower diagonal matrix L:\n");
	  printmat(matL, dimA, dimA);
	  printf("The upper diagonal matrix U:\n");
	  printmat(matU, dimA, dimA);
  }
}

//Gram-Schmidt Orthogonalization
	//this function takes a matrix, X,  consisting of columns x1, x2,....xp where each column has n elements
	//this function returns:
		//matrix, Q, containing the orthonormal equivalents q1, q1, ....qp of matrix X's columns
		//an upper triangular matrix R (optionally)
void GS_ortho(double *matX, double *matQ, double *matR, int n, int p, int getresids){
	int i, j;
	if (getresids) p++;
	double mylen[p]; //store length the lengths of our vectors w_j before rescaling in here
	double mydot; //store the dot project of w_(j-1), v_j for current column vector j in here
	double cumproj[n]; //we will store the projections in here (cumulative 1D column array)
		zeromat(cumproj, n);

	//Calculate the Q matrix, column by column
	for(j=0; j<p; j++){
		mydot = dots(matQ+((j-1)*n),matX+j*n, 1,1,n); //dot product of w_(j-1), v_j
  		for(i=0; i<n; i++){
			if (j==0) {
				matQ[j*n+i] = matX[j*n+i]; //1st column should be dealt with differently
			}
			else {
				cumproj[i] += mydot * matQ[(j-1)*n+i];
				matQ[j*n+i] = matX[j*n+i] - cumproj[i];
			}
		}
		//scale vectors
			mylen[j] = dots(matQ+j*n, matQ+j*n, 1, 1, n); /*length of vector wj is the dot product of jth col of matQ with itself*/
  			for(i=0; i<n; i++) {
				if (!(getresids && j==p)) matQ[j*n+i] = matQ[j*n+i] / mylen[j]; //do this everytime except for last col when getresid is specified (containing y values)
			}
	}
	//Calculate the R matrix
	if (getresids) p--; //we want this matrix created the same way even when we use the getresids option (must counteract earlier p++)
	for(j=0; j<p; j++){
  		for(i=0; i<p; i++){
			if (i==j) matR[j*p+i] = mylen[j]; //diagonal values
			else if (i<j) {
				matR[j*p+i] = dots(matQ+(i*n),matX+(j*n), 1,1,n); //dot product of col w_i and col v_j
			}
			else matR[j*p+i] = 0;
		}
	}

}

//find inverse of an upper-triangular matrix
void InvUTMat (double *matA, double *matI, int p){
	int i, j, k;
	double tmp;
	double product;
	double divisor;
	for (j=p-1; j>=0; j--){
		matI[j*p+j] = 1.0 / matA[j*p+j];
		for(k=j-1; k>=0; k--){
			for(i=k+1; i<=j; i++){
				product = (matA[i*p+k]*matI[j*p+i]);
				divisor = matA[k*p+k];
				tmp =  product / divisor;
				matI[j*p+k] -= tmp;
			}
		}
	}
}

//Regression using Gram-Schmidt Orthogonalization
	//slight modification since appearance in HW3 = now we are passing in a vector for Beta coefficients (so we can use them outside of function)
void GS_regress(double *Betas, double *matY, double *matX, int nrowX, int ncolX, int outtest){
	int getresid = 1;
	int c,i,j;
	int n = nrowX;
	int p = ncolX;

	//create an augmented X matrix so that we can use GS_ortho to get residuals
		int cellsX = nrowX*ncolX;
		int cellsR = ncolX*ncolX;
		int cellsXY = nrowX*(ncolX+1);
		double matXY[cellsXY];
		for(c=0;c<cellsXY;c++){
			if (c<cellsX) matXY[c]=matX[c];
			else matXY[c]=matY[c-cellsX];
		}
		if (outtest) printf("Augmented Matrix = XY combo\n");
	  	if (outtest) printmat(matXY,nrowX,ncolX+1);
	//run GS_ortho to obtain matrix Q, matrix R, and our residuals
		double matQresid[cellsXY]; //this will contain Q plus a column with the residuals
		double matR[cellsR];
		double matQ[cellsX];
		double Resids[nrowX];
		zeromat (matQresid, cellsXY);
	    zeromat (matR, cellsR);
		GS_ortho(matXY, matQresid, matR, nrowX, ncolX,1);

	//separate matrix Q from vector of residuals
			for(c=0;c<cellsXY;c++){
				if (c<cellsX) matQ[c]=matQresid[c];
				else Resids[c-cellsX]=matQresid[c];
			}
			if (outtest) printf("Matrix Qresid should contain matrix Q plus residuals\n");
    		if (outtest) printmat(matQresid,nrowX,ncolX+1);
    		if (outtest) printf("****Residuals\n");
    		if (outtest) printmat(Resids,nrowX,1);

	//examine and test output of GS_ortho
    	if (outtest) printf("Matrix Q (should be orthonormal)\n");
    	if (outtest) printmat(matQ,nrowX,ncolX);
    	if (outtest) printf("Matrix R (should be upper triangular)\n");
    	if (outtest) printmat(matR,ncolX,ncolX);
    	double matX_remake[cellsX];
    	multmat(matQ, matR, matX_remake, nrowX, ncolX, ncolX);
    	if (outtest) printf("Matrix X, recreated\n");
    	if (outtest) printmat(matX_remake,nrowX,ncolX);

	//Find the inverse of R
			double matRI[cellsR];
			zeromat(matRI, cellsR);
			InvUTMat (matR, matRI, ncolX);
			if (outtest) printf("This should be the inverse of R\n");
			if (outtest) printmat(matRI,ncolX,ncolX);
				//test inverse function
				double matID[cellsR];
				zeromat (matID, cellsR);
				multmat(matR,matRI,matID,ncolX,ncolX,ncolX);
				if (outtest) printf("This should be the Identity matrix\n");
				if (outtest) printmat(matID,ncolX,ncolX);

	//Beta estimates are R^(-1)trans(Q)Y
			double RIQtrans[cellsX]; //this will be a pxn matrix
			for (i=0; i<p; i++) {
				for(j=0; j<n; j++) {
		   			/*note: the i,j element of RIQtrans is the dot product of the ith row of RI with the jth col of transQ
		   		        and the jth col of transQ is the same as the jth row of Q*/
					RIQtrans[j*p+i] = dots(matRI+i, matQ+j, ncolX, nrowX, ncolX);
				}
   			}
			multmat(RIQtrans, matY, Betas, ncolX, nrowX, 1);
			if (outtest) printf("****Beta Coefficients \n");
			if (outtest) printmat(Betas,ncolX,1);

	//MSE is the sum of the squares of the residuals, divided by (n-p)
			double MSE = 0;
			for(i=0; i<nrowX; i++){
				MSE += (Resids[i] * Resids[i])/(nrowX-ncolX);
			}
			if (outtest) printf("****MSE: %lf\n\n",MSE);

	//Std errs for Betas are found using S = R^(-1) and the MSE
			double covBetas[cellsR]; //this will be a pxp matrix containing our Beta variance info
			double SStrans[cellsR]; //this will be a pxp matrix
			double sdBetas[ncolX];
			for (i=0; i<p; i++) {
				for(j=0; j<p; j++) {
		   			/*note: the i,j element of SStrans is the dot product of the ith row of RI with the jth col of transRI
		   		        and the jth col of transRI is the same as the jth row of RI*/
					SStrans[j*p+i] = dots(matRI+i, matRI+j, ncolX, ncolX, ncolX);
					covBetas[j*p+i] = sqrt(SStrans[j*p+i] * MSE);
					if (i==j) sdBetas[i]=covBetas[j*p+i]; //take diagnoal elements of covariance matrix
				}
   			}
   			if (outtest) printf("****Std Errs for Beta Coefficients \n");
			if (outtest) printmat(sdBetas,ncolX,1);

   	//Print results in sensible order
   		if (!outtest) {
   			printf("****Beta Coefficients \n");
   			printmat(Betas,ncolX,1);
   	   		printf("****Std Errs for Beta Coefficients \n");
			printmat(sdBetas,ncolX,1);
		    printf("****Residuals\n");
    		printmat(Resids,nrowX,1);
    		printf("****MSE: %lf\n\n",MSE);
		}
}

