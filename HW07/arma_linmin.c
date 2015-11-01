#include "arma_headers.h"

#define TOL 2.0e-4

int ncom = 0;
double *pcom = 0,*xicom = 0;
//double (*nrfunc)();

void arma_linmin(double sdata[], int nobs, int ar_p, int ma_q, double *p,double *xi,int n,double *fret)
//double *p,*xi,*fret;
//long n;
{
	long j;
	double xx,xmin,fx,fb,fa,bx,ax,tol=1.0e-4;
	double arma_brent(),f1dim(),*vector();
	void arma_mnbrak(),free_vector();
	int z;
	double alpha[ar_p];
    double beta[ma_q];

	ncom=n;
	pcom = vector(n);
	xicom = vector(n);
	//nrfunc = func;
	for (j=0;j<n;j++){
		pcom[j] = p[j];
		xicom[j] = xi[j];
		}

	ax = 0.;
	xx = 1.0;
	bx = 2.0;
	arma_mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim,sdata,nobs,ar_p,ma_q);
	*fret=arma_brent(ax,xx,bx,f1dim,TOL,&xmin,sdata,nobs,ar_p,ma_q);
	for(j=0;j<n;j++){
		xi[j] *= xmin;
		p[j] += xi[j];
	}
	free_vector(xicom);
	free_vector(pcom);
}

double f1dim(double x, double sdata[], int nobs, int ar_p, int ma_q)
{
	int j;
	double f,*xt,*vector();
	void free_vector();
	int z;
	double alpha[ar_p];
    double beta[ma_q];

	xt = vector(ncom);
	for(j=0;j<ncom;j++)xt[j] = pcom[j] + x * xicom[j];
	//f = (*nrfunc)(xt);
	for(z=0;z<(ar_p+ma_q);z++){
		if (z<ar_p) alpha[z] = xt[z];
		else beta[z-ar_p] = xt[z];
	}
	doarma(sdata, alpha, beta, &nobs, &ar_p, &ma_q, &f);
	free_vector(xt);
	return f;
}

