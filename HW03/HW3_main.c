#include <stdint.h>
#include <stdio.h>
#include <math.h>
FILE *infile;

//This program will use functions located in HW3_functions.c.  I must include corresponding function protocols here
   int is_d_equal (double, double, double);
   void var_mean (double *, int, double *, double *);
   void max_min (double *, int, double *, double *);
   uint64_t uniform_randoms (uint64_t, const int, uint64_t, int, double *, int, int);
   uint64_t normal_randoms (uint64_t, int, uint64_t, int, double *, int, int);


//dot product of vectors;
   //ix and iy are strides. use stride of nrow to use a row of the input matrix, and stride of 1 to use a column;
   //n is the number of elements in each of the vectors (i.e., the row or col of the input matrices);
double dots(double *x, double *y, int ix, int iy, int n){
	double sum = 0;
	int i;
	for(i=0; i<n; i++, x += ix, y += iy) sum += *x * *y;
	return(sum);
}

//matrix multiplication (generates a matrix c, which will be n x p);
   //matrix a is n x m
   //matrix b is m x p
void multmat(double *a, double *b, double *c, int n, int m, int p) {
   int i, j;
   for (i=0; i<n; i++) {
	   for(j=0; j<p; j++) {
		   /*note that a+i is the ith row of a, and  b+j*m is the jth col of b*/
		   c[j*n+i] = dots(a+i, b+j*m, n, 1, m); /*element C_ij is the dot product of the ith row of A with the jth col of B*/
	   }
   }
}

//multiplication of the inverse of a matrix with itself (generates a matrix c, which will be m x m);
   //matrix trans(x) is m x n
   //matrix x is n x m
void xtransx (double *x, double *c, int n, int m) {
   int i, j, t;
   for (i=0; i<m; i++) {
	   for(j=0; j<=i; j++) {
		   /*note: the i,j element of c is the dot product of the ith row of trans(x) with the jth col of x
		           amd the ith row of trans(x) is the same as the ith col of x*/
		   c[j*m+i] = dots(x+(i*n), x+(j*n), 1, 1, n);
	   }
   }
}

//prints an n x m matrix;
void printmat(double *a, int n, int m){
	printf("This %d by %d matrix contains the following elements:\n", n, m);
	int i, j;
	for (i=0; i<n; i++){
		printf("Row %d: ", i);
		for (j=0; j<m; j++){
			printf("%lf  ",a[j*n+i]);
		}
		printf("\n");
	}
	printf("\n");
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

void zeromat (double *matA, int ncells) {
	int i = 0;
	for(i=0;i<ncells;i++){
		matA[i] = 0;
	}
}

//F-statistic (nrow is the number of grps, and ncol is the number of obs per grp)
double F_stat(double *mat, int nrow, int ncol, int printme){
	double fstat = 0;
	int N = ncol*nrow;
	int mycol;
	int i, j, grp;
	double grpmeans[nrow];
	double allmeans = 0;
	for(grp=0; grp<nrow; grp++){
		double t_mean=0;
		for(mycol=0; mycol<ncol; mycol++) t_mean+=(1.0/ncol)*mat[mycol*nrow+grp];
		grpmeans[grp]=t_mean;
		allmeans += (1.0/nrow)*t_mean;
	}
	if (printme) printf("Means for each grp: \n");
	if (printme) printmat(grpmeans, 1, grp);
	if (printme) printf("Overall mean: %lf \n", allmeans);
	//use same lettering as in instructions
		int k = nrow;
		int n = ncol;
		double numerator = 0;
		double denominator = 0;
	for(i=0;i<k;i++){
		numerator += n*pow((grpmeans[i]-allmeans),2)/(k-1);
		for(j=0;j<n;j++){
			denominator += pow((mat[j*nrow+i]-grpmeans[i]),2)/(N-k);
		}
	}
	fstat = numerator/denominator;
	return fstat;
}

void simulation(uint64_t modulus, int a, uint64_t seed, int nrow, int ncol, double mycorr, int nsim, double *fstats_cor, double *fstats_unc, int printme) {
  int s;
  for(s=0; s<nsim; s++){

	//generate a matrix containing N(0,1) random numbers;
      int cell_total = nrow * ncol; //for efficiency and compatibility with other programs, we will store these matrices as 1D arrays
	  double rmat1[cell_total];
      seed = normal_randoms (modulus, a, seed, cell_total, rmat1, 0, 0);
      if (printme) printf("%d uncorrelated observations for %d groups\n", ncol, nrow);
	  if (printme) printmat(rmat1, nrow, ncol);

	//generate a matrix containing 5 rows consisting of 6 correlated observations, with specified level of correlation
	  //the correlation matrix for each of our groups of data will be as follows
	     double corrmat[]={
	    	 1, mycorr, pow(mycorr,2), pow(mycorr,3), pow(mycorr,4), pow(mycorr,5),
			 mycorr, 1, mycorr, pow(mycorr,2), pow(mycorr,2), pow(mycorr,4),
	    	 pow(mycorr,2), mycorr, 1, mycorr, pow(mycorr,2), pow(mycorr,3),
	    	 pow(mycorr,3), pow(mycorr,2), mycorr, 1, mycorr, pow(mycorr,2),
	    	 pow(mycorr,4), pow(mycorr,3), pow(mycorr,2), mycorr, 1, mycorr,
	    	 pow(mycorr,5), pow(mycorr,4), pow(mycorr,3), pow(mycorr,2), mycorr, 1};
	    	 if (printme) printf("Correlation matrix for each group\n");
	    	 if (printme) printmat(corrmat, ncol, ncol);
	  //find the cholesky decomposition
		  int dimC = ncol;
		  int corrcells = dimC*dimC;
		  double matLcorr[corrcells];
		  double matUcorr[corrcells];
		  zeromat (matLcorr, corrcells);
		  if (printme) printmat(matLcorr, dimC, dimC);
  		  zeromat (matUcorr, corrcells);
		  if (printme) printmat(matUcorr, dimC, dimC);
		  cholesky(dimC, corrmat, matLcorr, matUcorr, 0);
		  //verify result
			     double matcheck2[dimC*dimC];
			     multmat(matLcorr, matUcorr, matcheck2, dimC, dimC, dimC);
			     if (printme) printf("This should be our original matrix corrmat:\n");
    			 if (printme) printmat(matcheck2, dimC, dimC);
	  //multiply L (lower triangular from cholesky decomp) by the vectors of random normals for each group;
	  	  double rmatcorr[cell_total]; //this matrix will contain my correlated observations for each group
	      double t_rmat1[ncol];
	      double t_rmatcorr[ncol];
	      int grp, newi, obs;
	      for(grp=0; grp<nrow; grp++){
	    	 newi=0;
	    	 //get vector of the uncorrelated obs for group grp
	    	 	for(obs=grp; obs<cell_total; obs+=nrow){
	    			 t_rmat1[newi]=rmat1[obs];
	    			 newi++;
	    	 	}
	    	 //multiply with lower triangular matrix
	         	multmat(matLcorr, t_rmat1, t_rmatcorr, ncol, ncol, 1);
	         newi=0;
	         //put result into grp row of output matrix
	    	 	for(obs=grp; obs<cell_total; obs+=nrow){
	    			 rmatcorr[obs] = t_rmatcorr[newi];
	    			 newi++;
	    	 	}
	      }
	    //examine matrix of correlated numbers (each row should contain correlated elements);
	    	 if (printme) printf("%d correlated observations for %d groups\n", ncol, nrow);
	    	 if (printme) printmat(rmatcorr, nrow, ncol);

		//Calculate the F-statistics and preserve them in array
			double fstat = 0;
			//for uncorrelated observations
				fstat = F_stat(rmat1, nrow, ncol, 0);
				if (printme) printf("F-stat with uncorrelated observations is %lf \n",fstat);
				fstats_unc[s]=fstat;
			//for correlated observations
				fstat = F_stat(rmatcorr, nrow, ncol, 0);
				if (printme) printf("F-stat with correlated observations is %lf \n",fstat);
				fstats_cor[s]=fstat;

  } /*end simulation loop*/
}

//This function sort sorts an array of n doubles in place (taken from Philz webpage)
	// z is a pointer to the array of doubles to be sorted
	// n is the number of elements in z
void sort(double *z,long n) {
  long nodd,i,j,j1,j2,k;
  double v;
  nodd = n - 1;
  for(i=0;i<nodd;i++)
     {if(z[i] > z[i+1])
	{v = z[i + 1];
	 if(z[0] >= v)k = -1;
	 else
	   {j1 = 0;
	    j2 = i;
	    k = i / 2;
	    while(z[k] != v)
	      {k = (j1 + j2) / 2;
	       if(k == j1)break;
               else
		 {if(z[k] > v)j2 = k;
		  else if(z[k] < v)j1 = k; }
              }
            }
         for(j=i;j>k;j--)z[j + 1] = z[j];
         z[k + 1] = v;
        }
     }
}

//This function can be used for scaler multiplication of a matrix.  matrix contains nelements number of elements
void smultmat(double *a, double scale, int nelements) {
   int i;
   for (i=0; i<nelements; i++) {
		a[i] = scale * a[i];
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
void GS_regress(double *matY, double *matX, int nrowX, int ncolX, int outtest){
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
			double Betas[ncolX]; //this will be a px1 matrix
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


/******************************************************************************/
/******************************************************************************/


int main(int argc, char *argv[]){

	//exercise control over which parts of program are run;
	  int run_Q1 = 1;
	  int run_Q2 = 1;
	  int run_Q3 = 1;
	  int run_Q4 = 1;

	  int nrow, ncol, cell_total, nsim;
	  int i, j, k, n, g;
	  double mycorr;

    //initialize values for random number generation
      uint64_t modulus = 1;
      modulus <<= 31;
      modulus*=2;
      const int a = 8003;
      uint64_t seed = 3498725311U;


if (run_Q1) {
printf("\n*************** QUESTION 1 ****************\n\n");

//Test the Cholesky Decomposition (remember that j*nrow+i represents the i,j element)
  int dimA = 3;
  double matA [] = { 2, -1, 0, -1, 2, -1, 0, -1, 2}; //symmetric positive definite matrix
  double matL [] = { 0, 0, 0, 0, 0, 0, 0, 0, 0}; //this is the lower triangular matrix that we want to create. stuff with zeros for now.
  double matU [] = { 0, 0, 0, 0, 0, 0, 0, 0, 0};
  cholesky(dimA, matA, matL, matU, 1);
  //verify the results are good
    double matcheck1[dimA*dimA];
    multmat(matL, matU, matcheck1, dimA, dimA, dimA);
    printf("This should be our original matrix A:\n");
    printmat(matcheck1, dimA, dimA);

} /*end of run_Q1*/



if (run_Q2) {
printf("\n*************** QUESTION 2 ****************\n\n");

	//Simulation
	    int printme = 0; //zero if you wish to suppress extra output
		nrow = 5; //number of groups
		ncol = 6; //number of observations per group
		mycorr = .7; //base level of correlation between observations belonging to the same group
		nsim = 10000; //number of simulations to run
		double fstats_cor[nsim];
		double fstats_unc[nsim];
	    simulation(modulus, a, seed, nrow, ncol, mycorr, nsim, fstats_cor, fstats_unc, 0);
	//Analyze simulation results
    	//find f_stat values that correspond to 90th, 95th, and 99th percentiles
     	double pct90 = (.90*nsim) - ((90*nsim) % 100) - 1;
     	double pct95 = (.95*nsim) - ((95*nsim) % 100) - 1;
        	double pct99 = (.99*nsim) - ((99*nsim) % 100) - 1; /*this one doesn't work*/
       		pct99 = .99*nsim - 1;
       	printf("observations to look at: %d, %d, %d \n", (int)pct90, (int)pct95, (int)pct99);
        //uncorrelated observations
	    	sort(fstats_unc, nsim);
    		if (printme) printf("Sorted F-stats from simulation with uncorreleated obs:\n");
    		if (printme) printmat(fstats_unc, 1, nsim);
    		printf("Uncorrelated obs percentiles: 90th=%lf, 95th=%lf, 99th=%lf\n",fstats_unc[(int)pct90],fstats_unc[(int)pct95],fstats_unc[(int)pct99]);
    	//correlated observations
	    	sort(fstats_cor, nsim);
    		if (printme) printf("Sorted F-stats from simulation with correleated obs:\n");
			if (printme) printmat(fstats_cor, 1, nsim);
    		printf("Correlated obs percentiles: 90th=%lf, 95th=%lf, 99th=%lf\n",fstats_cor[(int)pct90],fstats_cor[(int)pct95],fstats_cor[(int)pct99]);

} /*end run_Q2*/


if (run_Q3) {
printf("\n\n*************** QUESTION 3 ****************\n\n");

	//create matrices for testing
	  int ncolX3 = 4; //this is the value for p
	  int nrowX3 = 5; //this is the value for n
	  int gscells3 = ncolX3*nrowX3;
	  double matX3[] = { 1, 5, 0, 1, 2, 4, 6, 2, 13, 5, 14, 10, 7, 0, 33, 17, 8, 4, 9, 10};
	  printf("Matrix X3, GS orthogonalization original\n");
      printmat(matX3,nrowX3,ncolX3);
	  double matQ3[gscells3];
	  double matR3[gscells3];
	  zeromat (matQ3, gscells3);
	  zeromat (matR3, gscells3);

	//run GS orthogonalization and check results by verifying that X=QR
    	GS_ortho(matX3, matQ3, matR3, nrowX3, ncolX3, 0);
    	printf("Matrix Q3 (should be orthonormal)\n");
    	printmat(matQ3,nrowX3,ncolX3);
    	printf("Matrix R3 (should be upper triangular)\n");
    	printmat(matR3,ncolX3,ncolX3);
    	double matX3_remake[gscells3];
    	multmat(matQ3, matR3, matX3_remake, nrowX3, ncolX3, ncolX3);
    	printf("Matrix X3, recreated\n");
    	printmat(matX3_remake,nrowX3,ncolX3);


} /*end run_Q3*/

if (run_Q4) {
printf("\n\n*************** QUESTION 4 ****************\n");


 //Run regression using matrix established here in C code
 printf("\n\n***** Regression 1: data entered in C code *****\n\n");
	//create matrices for testing
	  int ncolX = 4; //this is the value for p
	  int nrowX = 5; //this is the value for n
	  int gscells = ncolX*nrowX;
	  double matX[] = { 2, 4, 6, 2, 13, 5, 14, 10, 7, 0, 33, 17, 8, 4, 9, 10, 1, 5, 0, 1 };
	  printf("Matrix X, GS regression original\n");
	  printmat(matX,nrowX,ncolX);
	  double matY[] = { 4, 9, 5, 2, 13 };
	  printf("Matrix Y for GS regression\n");
      printmat(matY,nrowX,1);
    //Run Gram-Schmidt Regression on matrices created, above
	  GS_regress(matY, matX, nrowX, ncolX, 0);

  //Run regression using data read from text file
  printf("\n\n***** Regression 2: data read from text file *****\n\n");
	//read data from text file
		int i,j;
		int nrowX2;
		int ncolX2;
		infile=fopen("X_data.txt", "r");//The first 2 elements of this data file are the number of rows and columns
		if (!infile) {
			printf("Cannot Open File");
			return 1;
		}
		fscanf(infile, "%d", &nrowX2);
		fscanf(infile, "%d", &ncolX2);

        double matX2[nrowX2*ncolX2];
		for (i=0; i<nrowX2; i++) {
			for (j=0; j<ncolX2 ; j++){
				fscanf(infile, "%lf", matX2+i+j*nrowX2);
			}
	    }
		fclose(infile);
		printf("Matrix X2 (from file), GS regression original\n");
		printmat(matX2,nrowX2,ncolX2);

		infile=fopen("Y_data.txt", "r");
		if (!infile) {
			printf("Cannot Open File");
			return(1);
		}
		double matY2[nrowX2];
		for (i=0; i<nrowX2; i++) fscanf(infile, "%lf", matY2+i);
	    fclose(infile);
	    printf("Matrix Y2 (from file) for GS regression\n");
	    printmat(matY2, nrowX2, 1);

	//Run Gram-Schmidt Regression on data read from file, above
	  GS_regress(matY2, matX2, nrowX2, ncolX2, 0);

} /*end run_Q4*/


}
