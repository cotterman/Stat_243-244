#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define LOC(a,b)  ((a) * p + (b))
#include "Stat244_hw1_header1.h" //this basically means to copy and paste code from HW4_header.h into the top of this code document
#include "timer.h"
FILE *infile;

/*For the LAPACK routine in this code to work, you must log onto an SCF machine (which have the LAPACK library already installed) and run the following:
	gcc -g Stat244_hw1_main2.c -o hw1exe -llapack
	./hw1exe -llpack > hw1_c_out.txt
(you can then copy the output file onto your home PC) */


/******************************************************************************/
/******************************************************************************/


//obsnum from data will become leader of clustercount + 1
void make_leader(int obsj, int *clustercount, double *leaders, double *matD, int obstot, int vars){
	int p;
	for(p=0;p<vars;p++){
 		leaders[*clustercount + p*obstot] = matD[obsj + p*obstot]; //row clustercount in leaders matrix should get row obsnum from matD
	}
	matD[obsj + (vars)*obstot] = *clustercount; // last column of matD gets value = clustercount, indexed starting with zero
	(*clustercount)++; //parentheses to get order of operations right, so that we increment value at clustercount rather than the address
}


//obtain cluster number of leader that minimizes distance to observation j, and what that distance is
void min_dist(int obsj, double *mydist, int *mycluster, double *leaders, double *matD, int clustercount, int obstot, int vars) {
	int p, k;
	for (k=0; k<clustercount; k++){
		double t_mydist = 0;
		for (p=1; p<vars; p++){
			double tmp = leaders[k + p*obstot] - matD[obsj + p*obstot];
			tmp = tmp*tmp;
			t_mydist += tmp;
		}
		if (k==0) { //if only one equal sign then k will get assigned the value of 0, which would be wrong
			*mydist = t_mydist;
			*mycluster = k; //I will label my clusters starting with 0
		}
		else if (t_mydist < *mydist) {
			*mydist = t_mydist;
			*mycluster = k;
		}
	}
}

extern dsyevd_(char* jobz, char* uplo, int* n, double* a, int* lda, double* w,
               double* work, int* lwork, int* iwork, int* liwork, int* info); //I am using a fortran function from LAPACK here

void dsyev(double *A, int n, double *E){
  char jobz, uplo;
  int lda, lwork, info, i;

  jobz = 'V';
  uplo = 'U';
  lda = n;

  lwork = 3*n-1;
  double work[lwork];

  dsyev_(&jobz, &uplo, &n, A, &lda, E, work, &lwork, &info);
}


int main(int argc, char *argv[]){

	//exercise control over which parts of program are run;
	  int run_Q1 = 1;
	  int run_Q2 = 1;


if (run_Q1) {
printf("\n\n*************** QUESTION 1 (Calculating Eigenvalues & Eigenvectors) ****************\n\n");

	double sum_triangle = 0;
	double eps1 = .0000001; //will use in stopping rule
	int i, j, t, timeme;
	int mytimes = 100000; //this is the number of times we will run each routine (in order to time them)

	//create matrix for which we want to know the eigenvalues and eigenvectors.  Should be full-rank and symmetric.
	  int dimA = 5;
	  int cellsA = dimA*dimA; //number of columns and rows that matrix A has
	  double matA[cellsA];
	  zeromat (matA, dimA*dimA);
	  double triA[] = {7, 15, 11, 1, 2, 4, 6, 7, 13, 5};
	  double diagA[] = {14, 10, 16, 8, 9};
	  //fill the diagonal of matrix A
	  	for(i=0;i<dimA;i++){
				matA[i*dimA+i] = diagA[i];
	  	}
	  //fill in the lower triangle of matrix A
	  	int n = 0;
	  	for (i=0; i<dimA; i++) { //cols
	  		   		for(j=0; j<dimA; j++) { //rows
	  			 		if (i!=j && i>j) {
							matA[i+dimA*j] = triA[n];
							matA[j+dimA*i] = triA[n];
							n++;
						}
	  				}
	  	}
	  printf("Matrix A, original (should be full-rank and symmetric)\n");
      printmat(matA,dimA,dimA);


	/////the power method

	TIME0
	  double matT[cellsA];
	  double matR[cellsA];
	  double matC[cellsA];
	  double matX[cellsA];
	  for (timeme=0; timeme<mytimes; timeme++){
     	 for (i=0; i<cellsA; i++) {
	 	 		  matX[i] = matA[i];
	 	 }
     	 while (1) {
			zeromat(matT, cellsA);
			zeromat(matR, cellsA);
			multmat(matX, matA, matT, dimA, dimA, dimA); //multiply A * X (=T)
			gs(matT,matR,dimA,dimA); //find QR decomposition of T using Gram-Schmidt.  Resulting orthogonal matrix will be contained in new T.
			//see if stopping rule says to stop
			sum_triangle = 0;
			for (i=0; i<dimA; i++) {
	 	  		for(j=0; j<dimA; j++) {
			 		if (i != j ) //note: fabs is required if you want to return a double (abs will round to the nearest integer)
			 		    sum_triangle = sum_triangle + fabs(matR[j*dimA+i]); //add the values that are not on the diagonal
				}
			}
			for (i=0; i<cellsA; i++){
				matX[i]=matT[i];
			}
			if(sum_triangle < eps1)
				break;
	 	 }
  	   }
  	printf("Time to run the Power Method %i times", timeme);
	TIME1(" ")
	printf("\n");

	//check results
    	printf("Matrix R (should be diagonal matrix)\n");
    	printmat(matR,dimA,dimA);
    	zeromat (matC, cellsA);
    	xtransx(matX, matC, dimA, dimA); // X'X should be the identity matrix if X is orthonormal
    	printf("Matrix X'X (should be the identity matrix)\n");
    	printmat(matC, dimA, dimA);
    //report results
        printf("Matrix X (eigenvectors appear as columns)\n");
    	printmat(matX,dimA, dimA);
    	printf("The eigenvalues are as follows:\n");
			for (i=0; i<dimA; i++){
				printf("%lf  ",matR[i*dimA+i]);
			}
			printf("\n \n");


	////////now use the standard routine from LAPACK
	TIME0
	double E[dimA];
	for (timeme=0; timeme<mytimes; timeme++){
		for (i=0; i<cellsA; i++) {
			matX[i] = matA[i];
	 	}
    	dsyev(matX, dimA, E);
		}
	printf("Time to run the LAPACK DSYEV Routine %i times", timeme);
	TIME1(" ")
	printf("\n");

	//report results
		printf("Matrix from LAPACK (eigenvectors appear as rows)\n");
		printmat(matX,dimA, dimA);
		printf("The eigenvalues from LAPACK are as follows:\n");
					for (i=0; i<dimA; i++){
						printf("%lf  ",E[i]);
					}
		printf("\n");

} /*end run_Q1*/

if (run_Q2) {
printf("\n\n*************** QUESTION 2 (Cluster Analysis) ****************\n");

	int i, j, p, n, k;
	double myd, mytemp;
	double colmean = 0;
	double colvar = 0;
	int printme = 1;
	int printextra = 0;

	//read in the crime dataset
		//inputs obtained from background and data files
		int obstot = 16; //number of observations
		int vars = 8; //number of variables
		int cellsD = obstot * vars;
		int cellsDP = cellsD + obstot;
		double mycol[obstot]; //temp place to store each column
		infile=fopen("crime", "r"); //1st col is city code, remaining columns contain crime stats
			if (!infile) {
				printf("Cannot Open File");
				return 1;
			}
            double matD[cellsD + obstot]; //matD will hold data, plus extra column to indicate cluster membership
			for (i=0; i<obstot; i++) {
				for (n=0; n<vars; n++) {
					fscanf(infile, "%lf", &myd); //fscanf reads row-by-row
					matD[i+obstot*n] = myd;
				}
		    }
			fclose(infile);
			printf("Crime Data Matrix (from file)\n");
			printmat(matD,obstot,vars);

	//normalize the variables by subtracting their mean and dividing by their standard deviation
	  for (p=1; p<vars; p++) {
	 	 colmean = 0;
	 	 colvar = 0;
	 	 for (j=0; j<obstot; j++){
		 	 mycol[j] = matD[j+obstot*p]; //fill array with variable p values
		 }
  	 	 var_mean(mycol, obstot, &colmean, &colvar);
  	 	 for (j=0; j<obstot; j++){
		 	matD[j+obstot*p] = (matD[j+obstot*p] - colmean) / sqrt(colvar); //fill array with variable p values
		 }
  	  }
  	  printf("Normlized Crime Data Matrix\n");
	  printmat(matD,obstot,vars);

	double threshold; //start a new cluster if more than this distance away from existing leaders
	for(threshold=1;threshold<31;threshold++){

	//reset matrices and other values
		for (i=0; i<obstot; i++) {
				matD[i+vars*obstot] = -999; //this is the column that will eventually contain the cluster number
		    }
		double leaders[cellsD]; //to hold cluster leaders (and their covariate values)
		  memset(leaders, 0, sizeof(leaders)); //fills leaders with zeros
		double avgclust[cellsD]; //to hold the average x-variable values for each of the clusters
		  memset(avgclust, 0, sizeof(avgclust));
		double errtot = 0; //to hold the within cluster sum of squares
		int clustercount = 0; //current total number of clusters created
		double mydist = 0; //to hold distance to closest cluster leader
		int mycluster = 0; //to hold cluster number of leader that minimizes distance

	//The leader algorithm
		make_leader(0, &clustercount, leaders, matD, obstot, vars); //designate first observation as the leader of a cluster
		for (j=1; j<obstot; j++) {
			// Determine min distance between observation j and the leaders of all existing clusters (return min distance and corresponding cluster number)
			    min_dist(j, &mydist, &mycluster, leaders, matD, clustercount, obstot, vars);
			// If the min distance is less than or equal to t, then place observation j in the cluster with a leader closest to it.
			    if (mydist <= threshold) matD[j + vars*obstot] = mycluster; // last column of matD gets value = mycluster for observation obsnum
			// If the the minimum of these distances is greater than t, then create a new cluster with observation j as its leader.
			    else make_leader(j, &clustercount, leaders, matD, obstot, vars);
		}

		//Find the mean values of the x variables for each cluster
		for (k=0; k<clustercount;k++){ //for each cluster
			for (j=0;j<obstot;j++){ //loop through observations (rows) in matD
				if (matD[j+vars*obstot] == k) { //just for the observations belonging to cluster k
					avgclust[k] = avgclust[k] + 1; //first column of avgclust will contain the number of observations belonging to that row's cluster
					for(p=1; p<vars; p++){
						avgclust[k+p*obstot] += matD[j+p*obstot]; //running sum each of the x variables belonging to cluster k
					}
				}
			}
			for(p=1; p<vars; p++){
				if (avgclust[k]!=0) avgclust[k+p*vars] = avgclust[k+p*vars] / avgclust[k]; //divide running sum to get mean of x variables belonging to cluster k
			}
		}
		//Calculate the within cluster sum of squares
			errtot = 0;
			for (j=0;j<obstot;j++){ //loop through observations (rows) in matD
				for (p=1; p<vars; p++){
					k = matD[j+vars*obstot]; //the row number in avgclust that matches the cluster number of observation j
					mytemp = matD[j + p*obstot] - avgclust[k + p*obstot] ;
					mytemp = mytemp*mytemp;
					errtot += mytemp;
				}
			}

	//Ouput results
		printf("*** RESULTS USING TRESHOLD %fl ***\n", threshold);
		if (printextra==1){
			printf("Crime Data with Cluster Group No.\n");
		  	printmat(matD,obstot,vars+1);
			printf("Cluster Leaders \n");
		  	printmat(leaders,obstot,vars);
		  	printf("Cluster Means \n");
		  	printmat(avgclust,obstot,vars);
		}
		if (printme==1){
	  	printf("Total Number of clusters: %i \n", clustercount);
	  	printf("Total within cluster sum of squares: %fl \n", errtot);
	  	printf("Members of each cluster: \n");
		  	for (k=0; k<clustercount;k++){ //for each cluster
		  		printf("Cluster %i -> ", k);
		  		for (j=0;j<obstot;j++){ //loop through observations (rows) in matD
					if (matD[j+vars*obstot] == k) printf("%fl ",matD[j]);
				}
		  		printf("\n");
			}
			printf("\n");
		}

	} //end of threshold loop


} /*end run_Q2*/


}
