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

void zeromat (double *matA, int ncells) {
	int i = 0;
	for(i=0;i<ncells;i++){
		matA[i] = 0;
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

