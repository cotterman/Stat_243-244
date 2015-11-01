/*
This is the wrapper for beatonsweep().
Since I am using the .C method for calling C from R, I need this wrapper.
The.C method requires that all of the arguments passed to the function be pointers.
But if we include this wrapper, then we don't need to modify the beatonsweep() function.

To create the shared object that R needs for accessing the function, I must then do the following:
1) install rtools.exe (see http://www.murdoch-sutherland.com/Rtools/) if not already done so
2) search "cmd" from the start menu and change evironmental variables (if necessary) for this command prompt
	You must make sure you can access the directories containing R and all of the other relevant tools.
	(see http://www.stat.columbia.edu/~gelman/stuff_for_blog/AlanRPackageTutorial.pdf for help)
3) place beatonsweep.c and beatonsweepr.c in C:\Users\Carolyn
4) type the following into the command prompt
	R CMD SHLIB beatonsweepr.c beatonsweep.c -o bsweep.dll
On a windows machine, the resulting shared object that is created will be called bsweep.dll


*/

void beatonsweep(double *, int, int *, int);

	beatonsweepr(double *matrix, int *p, int *indices, int *nsweep) {
	beatonsweep(matrix, *p, indices, *nsweep);

}
