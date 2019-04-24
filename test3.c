#include "eMath.h"

void fun_f1(double* ans,double x){
  *ans = x*x + x;
}

int main(void){
  double d;
  integrate_function(&d,fun_f1,5,7,100);
  printf("%lf\n",d);
}
