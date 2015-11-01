

//This function will sweep the columns of a pxp square matrix called "matrix"
//it is an "in place" algorithm in the sense that the original matrix will be written over with the sweep result
//the columns to be swept have indices contained in the vector called "indices"
	//note that column 1 is the 1st column (i.e., there is no column 0 by this convention) -- this is more intuitive for typical R users
//nsweep contains the number of columns to be swept and should correspond to the number of elements in "indices"

void beatonsweep(double *matrix, int p, int *indices, int nsweep)
{
	int i, j, k, c;
	double d, b;

	if(p < nsweep){
		printf("Number of columns to sweep exceeds matrix dimension.\n");
		exit(1);
	}
	for(c=0; c<nsweep; c++){
		if(p < indices[c] || indices[c] < 1) {
			printf("Matrix is %d x %d.  Cannot sweep column %d.\n", p, p, indices[c]);
			exit(1);
		}
		k = indices[c];
		d = matrix[(k-1)*p + (k-1)]; //d contains element (k,k) of matrix
		if(d < sqrt(0.00000001)){
			printf("Column %d is not linearly independent of previously swept columns.\n", k);
			exit(1);
		}
		for(i=0; i<p; i++){
			matrix[i*p + (k-1)] /= d; //element i,k
		}
		for(i=0; i<p; i++) {
			if(i!=(k-1)) {
				b = matrix[(k-1)*p + i];
				for(j=0; j<p; j++) {
					matrix[j*p + i] = matrix[j*p + i] - (b * matrix[j*p + (k-1)]);
				}
				matrix[(k-1)*p + i] = -b/d;
			}
		}
		matrix[(k-1)*p + (k-1)] = 1/d;
	}
}
