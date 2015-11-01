/*************************************************************
 function to sweep the columns of a matrix 

Arguments:
  a    double *    pointer to the matrix (must be square, 
		   usually symmetric)
  n       long     number of columns in the matrix a.
  tol   double     tolerance value (1.e-8 is a reasonable 
		   choice)
  sw     long *    pointer to indices (base 0) of columns to
		   be swept
  ns      long     number of indices (columns to be swept) in
		   sw.
  ierr   long *    pointer to return code:
		      0 -> no problem
		      i -> zero pivot for column i-1
***************************************************************/
void sweep(double *a,long n,double tol,long *sw,long ns,long *ierr)

{ int i,j,k,k1;
  double d,b,v;
  double *atmp,*wtmp;

  double fabs();

  *ierr = 0;
  for(k1=0;k1<ns;k1++) 
     {k = sw[k1];
      if(fabs(d = a[k * n + k]) < tol) /* zero out kth row */
	{atmp = a + k;
	 for(i=0;i<n;i++,atmp += n)*atmp = 0.;
	 atmp = a + k * n;
	 for(i=0;i<n;i++,atmp++)*atmp = 0.;
	 *ierr = k + 1;
	 return; }

      wtmp = a + k * n;
      for(i=0;i<n;i++)
	 {if(i == k)continue;
	  atmp = a + i * n;
	  for(j=0;j<n;j++)
	     {if(j == k)continue;
	      atmp[j] -= atmp[k] * wtmp[j] / d; }
         }


      atmp = a + k;
      for(i=0;i<n;i++,atmp += n)
	 {if(i == k)continue;
	  *atmp = -(*atmp) / d; }

      atmp = a + k * n;
      for(i=0;i<n;i++,atmp++)
	 {if(i == k)continue;
	  *atmp = *atmp / d; }

      a[k * n + k] = 1 / d;
     }
}
