#include"eNeural.h"

int main(void){
  eANN n;
  int i[4] = {2,3,5,1};
  create_net(&n,i,3);
  eMatrix in_data,out_data,aux;
  char in_str[100] = "0 0;0 1;1 0;1 1";
  char out_str[100] = "1;0;0;1";
  str_2_matrix(&in_data,in_str,4,2);
  str_2_matrix(&out_data,out_str,4,1);
  train_net(&n,in_data,out_data,1.0,3000);
  forward_net(&n,in_data);
  print_matrix(n.net[n.n_layers-1].o);
}
