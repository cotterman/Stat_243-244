#include <math.h>

#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY  1.0e-20
#define MAX(a,b) ((a) > (b)  ? (a) : (b))
#define SIGN(a,b)  ((b)  > 0.0 ? fabs(a) : -fabs(a))
#define SHFT(a,b,c,d)  (a)=(b);(b)=(c);(c)=(d);


//void mnbrak(ax,bx,cx,fa,fb,fc,func)
//double *ax, *bx, *cx, *fa, *fb, *fc;
//double (*func)();
void arma_mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc, double (*func)(),double sdata[],int nobs,int ar_p, int ma_q)
{
	double ulim,u,r,q,fu,dum;

	//*fa=(*func)(*ax);
	*fa=(*func)(*ax,sdata,nobs,ar_p,ma_q);
	//*fb=(*func)(*bx);
	*fb=(*func)(*bx,sdata,nobs,ar_p,ma_q);
	if(*fb > *fa){
		SHFT(dum,*ax,*bx,dum)
		SHFT(dum,*fb,*fa,dum)
	}
	*cx = *bx + GOLD * (*bx-*ax);
	//*fc=(*func)(*cx);
	*fc=(*func)(*cx,sdata,nobs,ar_p,ma_q);
	while(*fb > *fc){
		r=(*bx-*ax)*(*fb-*fc);
		q=(*bx-*cx)*(*fb-*fa);
		u = *bx-((*bx-*cx)*q-(*bx-*ax)*r) /
		    (2.*SIGN(MAX(fabs(q-r),TINY),q-r));
		ulim = *bx + GLIMIT * (*cx-*bx);
		if((*bx-u) * (u-*cx) > 0.0){
			//fu=(*func)(u);
			fu=(*func)(u,sdata,nobs,ar_p,ma_q);
			if(fu < *fc){
				*ax = *bx;
				*bx = u;
				*fa = *fb;
				*fb = fu;
				return;
   	 			}
			else if(fu > *fb){
				*cx=u;
				*fc=fu;
				return;
				}
			u = *cx + GOLD * (*cx-*bx);
			//fu = (*func)(u);
			fu = (*func)(u,sdata,nobs,ar_p,ma_q);
			}
		else if((*cx-u)*(u-ulim) > 0.0){
			//fu=(*func)(u);
			fu=(*func)(u,sdata,nobs,ar_p,ma_q);
			if(fu < *fc){
				SHFT(*bx,*cx,u,*cx + GOLD * (*cx - *bx));
				//SHFT(*fb,*fc,fu,(*func)(u));
				SHFT(*fb,*fc,fu,(*func)(u,sdata,nobs,ar_p,ma_q));
			}
		}
		else if((u-ulim)*(ulim-*cx) >= 0.0){
			u=ulim;
			//fu=(*func)(u);
			fu=(*func)(u,sdata,nobs,ar_p,ma_q);
		}
		else{
			u = *cx+GOLD*(*cx-*bx);
			//fu=(*func)(u);
			fu=(*func)(u,sdata,nobs,ar_p,ma_q);
			}
		SHFT(*ax,*bx,*cx,u)
		SHFT(*fa,*fb,*fc,fu)
		}
	return;
}
