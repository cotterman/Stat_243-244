#include <stdio.h>
//#include "../options.h" /*this is a file that I created for debugging help*/

main() {

int i_prior = 1; /*long integer*/
int i_new = 2;
short s_prior = 1; /*short integer*/
short s_new = 1;
unsigned int ui_prior = 1;
unsigned int ui_new = 1;
unsigned short us_prior = 1;
unsigned short us_new = 1;
float f_prior = 1; /*single precision floating point number*/
float f_new = 2;
float f_old = 1;
double d_prior = 1; /*double precision floating point number*/
double d_new = 2;
double d_old = 1;
double etry = 1;
double esum = 2;


/*find the largest long integer that can  be handled by my computer*/
while(i_new>0) {
  i_prior*=2;
  i_new*=2;
  //printf("i_prior is %d\n y is %d\n",i_prior, i_new); /*d is good for printing signed integers*/
}
int myval;
myval = (i_prior-1)*2 + 1;
printf("largest long integer is %d\n", myval);

/*find the largest short integer that can be handled by my computer*/
while(s_new>0) {
  s_new = s_prior * 2;
  //printf("s_prior is %d\n s_new is %d\n",s_prior, s_new);
  if (s_new>0) s_prior = s_new;
}
short myshort;
myshort = (s_prior-1)*2 + 1;
printf("largest short integer is %d\n", myshort);

/*find the largest unsigned [long] integer that can be handled by my computer*/
while(ui_new>0) {
  ui_new = ui_prior * 2;
  //printf("ui_prior is %u\n ui_new is %u\n",ui_prior, ui_new); /*u is good for unsigned ints*/
  if (ui_new>0) ui_prior = ui_new;
}
unsigned int myunsignedi;
myunsignedi = (ui_prior-1)*2 + 1;
printf("largest unsigned [long] integer is %u\n", myunsignedi);

/*find the largest unsigned short [integer] that can be handled by my computer*/
while(us_new>0) {
  us_new = us_prior * 2;
  //printf("us_prior is %u\n us_new is %u\n",us_prior, us_new); /*u is good for unsigned ints*/
  if (us_new>0) us_prior = us_new;
}
unsigned short myunsigneds;
myunsigneds = (us_prior-1)*2 + 1;
printf("largest unsigned short [integer] is %u\n", myunsigneds);

/*find the largest [single precision] floating point that can be handled by my computer*/
while(f_new>0) {
  f_new = f_prior * 2;
  //printf("f_prior is %f\n f_new is %f\n",f_prior, f_new); /*f is good for floats*/
  if (f_new==f_prior) break;
  f_old = f_prior;
  f_prior = f_new;
}
/*f_old contains the max exponent component of a float*/
printf("largest [single precision] exponent component is %f\n", f_old);
float max_exp = f_old;
float mantissa = 0;
float m_add = 1;
float fm_prior = 0;
float fm_new = 0;
float fm_old = 0;
while(mantissa>=0) {
  m_add = m_add/2;
  mantissa = mantissa + m_add;
  fm_new = max_exp*(1+mantissa);
  //printf("fm_prior is %f\n fm_new is %f\n",fm_prior, fm_new); /*f is good for floats*/
  if (fm_new==fm_prior) break;
  fm_old = fm_prior;
  fm_prior = fm_new;
}
/*fm_old contains the max value of a float -- test it*/
float toohigh = fm_old + .5;
printf("largest [single precision] float is %f\n", fm_old);


/*find the largest double [precision floating point] that can be handled by my computer*/
while(d_new>0) {
  d_new = d_prior * 2;
  //printf("d_prior is %f\n d_new is %f\n",d_prior, d_new); /*f is good for floats*/
  if (d_new==d_prior) break;
  d_old = d_prior;
  d_prior = d_new;
}
/*d_old contains the max exponent component of a double [precision floating point]*/
printf("largest double [precision floating point] exponent component is %f\n", d_old);
double dmax_exp = d_old;
double dmantissa = 0;
double dm_add = 1;
double dm_prior = 0;
double dm_new = 0;
double dm_old = 0;
while(dmantissa>=0) {
  dm_add = dm_add/2;
  dmantissa = dmantissa + dm_add;
  dm_new = dmax_exp*(1+dmantissa);
  //printf("dm_prior is %f\n dm_new is %f\n",dm_prior, dm_new); /*f is good for floats*/
  if (dm_new==dm_prior) break;
  dm_old = dm_prior;
  dm_prior = dm_new;
}
/*dm_old contains the max value of a float -- test it*/
double dtoohigh = dm_old + .5;
printf("largest double [precision floating point] is %f\n", dm_old);

/*find the value of machine epsilon*/
while(esum>1) {
  esum = 1 + etry;
  //printf("etry is %.32lf\n",etry); /*d is good for double-precision floats*/
  if (esum>1) etry = etry/2;
}
printf("machine epsilon is %.32lf\n", etry);


}




