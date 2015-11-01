#include <stdio.h>

double sum(double, double); /*this is a function protype that works*/
//double sum(); /*A function declaration done this way will not work*/
double mysum;

main() {
  mysum = sum(3,4);
  printf("my sum is %lf\n", mysum);

}

/*declaration in C is not necessarily a prototype*/
/*without proper prototype, C doesn't know what I am passing*/
