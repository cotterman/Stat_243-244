

# This file contains the bsweep() function for R.
# The function takes only 2 input arguments, but we need to calculate and pass all of the arguments that our C function requires.  
# We use the as.integer() function to ensure that the integer data types are consistent with those required for our C program 
# We use the storage.mode() function to preserve the matrix data structure while ensuring that the data types are consistent with those required for our C program.

bsweep = function(mat, vec) {
	nr = nrow(mat)
	nc = ncol(mat)
	nrv = nrow(vec)
	ncv = ncol(vec)

	if( nr == nc && (nrv==1 || ncv==1) ) {
		if(!is.loaded('bsweepr')) dyn.load('C:/Users/Carolyn/Documents/R_myfiles/win-library/bsweep.dll')
		storage.mode(mat) = 'double'
		storage.mode(vec) = 'integer'
		z = .C( 'beatonsweepr', result = mat, as.integer(nc), vec, as.integer(max(c(nrv, ncv))) )
		return(z$result)
	} else return(0)
}
		