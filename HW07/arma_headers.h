
//Including this header file will give me more specific error messages if I accidentally pass the incorrect parameters to any of my functions
	//Else, my program will stop with a "segmentation fault" message

double arma_brent(double ax, double bx, double cx, double (*f)(), double tol, double *xmin, double sdata[],int nobs,int ar_p, int ma_q);
void arma_linmin(double sdata[], int nobs, int ar_p, int ma_q, double *p,double *xi,int n,double *fret);
void arma_powell(double sdata[], double p[], int nobs, double **xi, int ar_p, int ma_q, double ftol, int *iter, double *fret);
double amotry(double **p, double y[], double psum[], long ndim, long ihi, double fac, double beta[], double alpha[], int nobs, int ar_p, int ma_q, double sdata[]);
void arma_mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc, double (*func)(),double sdata[],int nobs,int ar_p, int ma_q);
double f1dim(double x, double sdata[], int nobs, int ar_p, int ma_q);
