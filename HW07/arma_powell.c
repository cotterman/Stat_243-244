/* the following routines were adapted from those in Numerical Recipes in C,
	The Art of Scientific Computing, Press, et. al., Cambridge.  */

#include <math.h>
#define ITMAX 200

static float sqrarg;
#define SQR(a)  (sqrarg=(a),sqrarg*sqrarg)


//void arma_amoeba(double sdata[], int nobs, double **p, double y[], int ar_p, int ma_q, double ftol, int *nfunk)

//void powell(p,xi,n,ftol,iter,fret,func) -- original function prototype
	//double *p,**xi,ftol,*fret;
	//long n,*iter;
void arma_powell(double sdata[], double p[], int nobs, double **xi, int ar_p, int ma_q, double ftol, int *iter, double *fret)
	//p is the initial starting point ([1,...n])
	//n is the number of parameters
	//xi is the n x n pointer-array matrix whose columns are the set of directions n_i (usually unit vectors)
	//ftol is the fractional tolerance in the function value such that failure to decrease by more than this amount on one iteration signals doneness.
	//OUTPUT: p is set to the best point found
		//fret is the returned function value at p
		//iter is the number of iterations taken
{
	long i,ibig,j;
	double t,fptt,fp,del;
	double *pt,*ptt,*xit,*vector();
	void arma_linmin(),nrerror(),free_vector();
	int n = ar_p + ma_q; //number of parameters
	int z;
	double alpha[ar_p];
    double beta[ma_q];

	pt = vector(n);
	ptt = vector(n);
	xit= vector(n);
	//*fret = (*func)(p);
	for(z=0;z<(ar_p+ma_q);z++){
		if (z<ar_p) alpha[z] = p[z];
		else beta[z-ar_p] = p[z];
	}
	doarma(sdata, alpha, beta, &nobs, &ar_p, &ma_q, fret);
	for (j=0;j<n;j++)
		pt[j] = p[j];
	for(*iter=1;;(*iter)++){
		fp = *fret;
		ibig = 0;
		del = 0.0;
		for (i=0;i<n;i++){
			for(j=0;j<n;j++)
				xit[j] = xi[j][i];
			fptt = *fret;
			//linmin(p,xit,n,fret,func);
			arma_linmin(sdata, nobs, ar_p, ma_q, p,xit,n,fret);
			if(fabs(fptt-*fret) > del){
				del=fabs(fptt-*fret);
				ibig=i;
			}
		}
		if(2.0*fabs(fp-*fret) <= ftol*(fabs(fp)+fabs(*fret))){
			free_vector(xit);
			free_vector(ptt);
			free_vector(pt);
			return;
			}
		if(*iter == ITMAX)
			nrerror("Too many iterations in routine POWELL");
		for(j=0;j<n;j++){
			ptt[j] = 2.0 * p[j] - pt[j];
			xit[j] = p[j] - pt[j];
			pt[j] = p[j];
		}
		//fptt = (*func)(ptt);
		for(z=0;z<(ar_p+ma_q);z++){
			if (z<ar_p) alpha[z] = ptt[z];
			else beta[z-ar_p] = ptt[z];
		}
		doarma(sdata, alpha, beta, &nobs, &ar_p, &ma_q, &fptt);
		if (fptt < fp){
			t = 2.0 * (fp-2.*(*fret)+fptt) * SQR(fp-(*fret)-del)
				  - del * SQR(fp-fptt);
			if(t < 0.0){
				//linmin(p,xit,n,fret);
				arma_linmin(sdata, nobs, ar_p, ma_q, p,xit,n,fret);
			for(j=0;j<n;j++)xi[j][ibig] = xit[j];
			}
		}
	}
}

