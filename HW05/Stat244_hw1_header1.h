//this header file will be called by the program Stat244_hw1_main1.c

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

void zeromat (double *matA, int ncells) {
	int i = 0;
	for(i=0;i<ncells;i++){
		matA[i] = 0;
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

//This function performs a Gram - Schmidt orthogonalization of an n x p matrix, x.
	//The orthogonal (q) part of the decomposition is returned in x, and the lower triangular part in r.
	//Both x and r must be allocated by the calling program, and r should be filled with zeroes.
	//The matrix x is assumed to be stored by columns.
	//This function was provided by Phil Spector

void gs(x,r,n,p)
 double *x,*r;
 long n,p;
{
  long i,j,k;
  double t;
  double *xnow;

  for(j=0;j<p;j++)
     {t = 0.;
      xnow = x + j;
      for(i=0;i<n;i++,xnow += p)t += *xnow * *xnow;
      t = sqrt(t);
      r[LOC(j,j)] = t;
      xnow = x + j;
      for(i=0;i<n;i++,xnow += p)*xnow /= t;
      for(k=j+1;k<p;k++)
	 {t = 0;
	  xnow = x;
	  for(i=0;i<n;i++,xnow += p)t += xnow[j] * xnow[k];
          r[LOC(j,k)] = t;
	  xnow = x;
          for(i=0;i<n;i++,xnow +=p)xnow[k] -= t * xnow[j]; }
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


