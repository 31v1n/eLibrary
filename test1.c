#include "eMath.h"

int main(void){
  eVector a,v1,v2;
  char v[100] = "1 8 2 3 9";
  char t1[100]="1 2 4 11",t2[100]=" 3 5 12";
  str_2_vector(&a,v,5);
  str_2_vector(&v1,t1,4);
  str_2_vector(&v2,t2,3);
  print_vector(v1);
  print_vector(v2);
  merge_data(&a,v1,v2);
  print_vector(a);
}
