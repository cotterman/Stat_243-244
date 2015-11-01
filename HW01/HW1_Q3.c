#include <stdio.h>

void cdfchi_(int*, double*, double*, double*, double*, int*, double*);

int main(int argc, char *argv[]){

  double Q, P, X, bound, df;
  int status, which;
  char repeat = 'y';

  while (repeat=='y') {
    printf("Please enter the probability level: ");
    scanf("%lf",&Q);
    P = 1-Q;
    X = 9999;
    bound = 9999;
    status=-1;
    which=2; /*indicates that I am providing prob level and DFs and that I want to obtain chi2 value*/
    printf("\nPlease enter the degrees of freedom: ");
    scanf("%lf",&df);

    cdfchi_(&which, &P, &Q, &X, &df, &status, &bound);

    printf("\nThe chi-square value for probability level %6.5lf with %12.lf dfs is %7.5lf\n", Q, df, X);
    printf("\nWould you like to do that again? (y or n)");
    scanf("\n%c",&repeat);
  }

  printf("\nWell alrighty then.  Hope you got what you wanted.  Goodbye!\n");

  return 0;
}
