#include "eMath.h"

void fun_f(double *o,double y,double x){
  *o =x*y;
}
int main(void){
  eMatrix data;
  eVector x;
  linspace_mesh(&x,0,5,100);
  ode1_solver(&data,fun_f,1,x);
  eMatrix d1;
  init_matrix(&d1,2,x.length);
  trans_matrix(&data,data);
  d1.data[0] = x;
  d1.data[1] = data.data[0];
  trans_matrix(&d1,d1);
  save_matrix("myoutput.dat",d1);
}
